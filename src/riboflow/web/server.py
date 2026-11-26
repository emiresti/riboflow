"""Zero-dependency HTTP/JSON server for the RiboFlow interactive web dashboard."""

from __future__ import annotations

import json
import mimetypes
import webbrowser
from dataclasses import asdict
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Dict, Tuple

from ..harmonizer import compare_harmonization
from ..mutation import apply_synonymous_mutations, scan_critical_silent_mutations
from ..organisms import ORGANISMS, get_organism
from ..samples import SAMPLES
from ..sequence import analyze_sequence, parse_fasta_or_raw
from ..simulator import simulate_ribosome_traffic

STATIC_DIR = Path(__file__).parent / "static"


def _json_response(handler: BaseHTTPRequestHandler, status: int, payload: Any) -> None:
    body = json.dumps(payload).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("Cache-Control", "no-store")
    handler.end_headers()
    handler.wfile.write(body)


def handle_api_post(path: str, payload: Dict[str, Any]) -> Tuple[int, Dict[str, Any]]:
    """Pure request router for JSON API endpoints (easily unit-testable)."""
    try:
        if path == "/api/simulate":
            raw_seq = str(payload.get("sequence", ""))
            _, dna = parse_fasta_or_raw(raw_seq)
            org = get_organism(str(payload.get("organism", "ecoli")))
            init_rate = payload.get("initiation_rate")
            alpha = float(init_rate) if init_rate is not None else None
            duration = min(300.0, max(20.0, float(payload.get("duration", 120.0))))
            seed = int(payload.get("seed", 42))

            analysis = analyze_sequence(dna, org)
            sim = simulate_ribosome_traffic(
                dna,
                org,
                initiation_rate=alpha,
                duration_seconds=duration,
                capture_frames=160,
                seed=seed,
            )
            return HTTPStatus.OK, {
                "analysis": asdict(analysis),
                "simulation": asdict(sim),
            }

        if path == "/api/mutate":
            raw_seq = str(payload.get("sequence", ""))
            _, dna = parse_fasta_or_raw(raw_seq)
            org = get_organism(str(payload.get("organism", "ecoli")))
            raw_edits = payload.get("edits", {})
            edits = {int(k): str(v) for k, v in raw_edits.items()}
            init_rate = payload.get("initiation_rate")
            alpha = float(init_rate) if init_rate is not None else None
            seed = int(payload.get("seed", 42))

            comp = apply_synonymous_mutations(
                dna,
                org,
                edits=edits,
                initiation_rate=alpha,
                duration_seconds=100.0,
                seed=seed,
            )
            return HTTPStatus.OK, asdict(comp)

        if path == "/api/scan-mutations":
            raw_seq = str(payload.get("sequence", ""))
            _, dna = parse_fasta_or_raw(raw_seq)
            org = get_organism(str(payload.get("organism", "ecoli")))
            init_rate = payload.get("initiation_rate")
            alpha = float(init_rate) if init_rate is not None else None
            top_k = min(12, max(1, int(payload.get("top_k", 6))))
            seed = int(payload.get("seed", 42))

            hotspots = scan_critical_silent_mutations(
                dna,
                org,
                initiation_rate=alpha,
                top_k=top_k,
                seed=seed,
            )
            return HTTPStatus.OK, {"hotspots": [asdict(h) for h in hotspots]}

        if path == "/api/harmonize":
            raw_seq = str(payload.get("sequence", ""))
            _, dna = parse_fasta_or_raw(raw_seq)
            src = get_organism(str(payload.get("source", "human")))
            host = get_organism(str(payload.get("host", "ecoli")))
            init_rate = payload.get("initiation_rate")
            alpha = float(init_rate) if init_rate is not None else None
            seed = int(payload.get("seed", 42))

            report = compare_harmonization(
                dna,
                source=src,
                host=host,
                initiation_rate=alpha,
                duration_seconds=100.0,
                seed=seed,
            )
            return HTTPStatus.OK, asdict(report)

        return HTTPStatus.NOT_FOUND, {"error": f"Unknown API endpoint: {path}"}
    except ValueError as exc:
        return HTTPStatus.BAD_REQUEST, {"error": str(exc)}
    except Exception as exc:  # pragma: no cover - defensive guard
        return HTTPStatus.INTERNAL_SERVER_ERROR, {"error": str(exc)}


def get_metadata_payload() -> Dict[str, Any]:
    """Return built-in organisms and curated gene presets."""
    organisms_list = [
        {
            "key": o.key,
            "name": o.name,
            "common_name": o.common_name,
            "ribosome_type": o.ribosome_type,
            "mean_elongation_rate": o.mean_elongation_rate,
            "default_initiation_rate": o.default_initiation_rate,
            "ribosome_footprint": o.ribosome_footprint,
        }
        for o in ORGANISMS.values()
    ]
    samples_list = [asdict(s) for s in SAMPLES.values()]
    return {"organisms": organisms_list, "samples": samples_list}


class RiboFlowRequestHandler(BaseHTTPRequestHandler):
    def log_message(self, format: str, *args: Any) -> None:
        return

    def do_GET(self) -> None:
        clean_path = self.path.split("?", 1)[0]
        if clean_path == "/api/meta":
            _json_response(self, HTTPStatus.OK, get_metadata_payload())
            return

        rel = "index.html" if clean_path in ("/", "") else clean_path.lstrip("/")
        target = (STATIC_DIR / rel).resolve()
        if not str(target).startswith(str(STATIC_DIR.resolve())) or not target.is_file():
            self.send_error(HTTPStatus.NOT_FOUND, "File not found")
            return

        mime, _ = mimetypes.guess_type(str(target))
        content = target.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", f"{mime or 'application/octet-stream'}; charset=utf-8")
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(content)

    def do_POST(self) -> None:
        clean_path = self.path.split("?", 1)[0]
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length > 0 else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError:
            _json_response(self, HTTPStatus.BAD_REQUEST, {"error": "Invalid JSON body."})
            return

        status, response = handle_api_post(clean_path, payload)
        _json_response(self, status, response)


def run_server(host: str = "127.0.0.1", port: int = 8765, open_browser: bool = False) -> None:
    server = ThreadingHTTPServer((host, port), RiboFlowRequestHandler)
    url = f"http://{host}:{port}"
    print(f"RiboFlow Web Laboratory running at {url}")
    print("Press Ctrl+C to stop.")
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down RiboFlow server.")
    finally:
        server.server_close()
