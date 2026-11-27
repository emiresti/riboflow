"""Command-line interface for RiboFlow v2.0."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional

from . import __version__
from .folding import fold_5prime_window, unfold_5prime_synonymously
from .harmonizer import compare_harmonization
from .mutation import apply_synonymous_mutations, scan_critical_silent_mutations
from .organisms import ORGANISMS, get_organism
from .samples import SAMPLES
from .sequence import analyze_sequence, parse_fasta_or_raw
from .simulator import simulate_ribosome_traffic, sweep_phase_diagram
from .synbio import design_synbio_sequence


def _load_sequence(args: argparse.Namespace) -> tuple[str, str]:
    if getattr(args, "sample", None):
        sample_id = args.sample
        if sample_id not in SAMPLES:
            valid = ", ".join(SAMPLES.keys())
            raise ValueError(f"Unknown sample '{sample_id}'. Valid samples: {valid}")
        s = SAMPLES[sample_id]
        return s.name, s.sequence

    if getattr(args, "fasta", None):
        content = Path(args.fasta).read_text(encoding="utf-8")
        return parse_fasta_or_raw(content)

    if getattr(args, "sequence", None):
        return parse_fasta_or_raw(args.sequence)

    raise ValueError("Provide one of --sample, --fasta, or --sequence.")


def _add_sequence_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--sample",
        choices=sorted(SAMPLES.keys()),
        help="Built-in biological sample gene.",
    )
    group.add_argument("--fasta", help="Path to a FASTA or raw DNA/mRNA file.")
    group.add_argument("--sequence", help="Raw coding DNA/mRNA sequence string.")


def _cmd_analyze(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    org = get_organism(args.organism)
    res = analyze_sequence(seq, org)

    if args.json:
        print(json.dumps(asdict(res), indent=2))
        return 0

    print(f"\nRiboFlow Sequence Analysis — {name}")
    print("=" * 64)
    print(f"Organism            : {org.name} ({org.ribosome_type})")
    print(f"Length              : {res.codon_count} codons ({len(res.sequence)} nt)")
    print(f"CAI                 : {res.cai:.4f}")
    print(f"GC / GC3            : {res.gc_percent:.1f}% / {res.gc3_percent:.1f}%")
    print(f"Mean Elongation     : {res.mean_elongation_rate:.2f} codons/s")
    print(f"%MinMax Range       : {res.min_minmax:+.1f} to {res.max_minmax:+.1f}")
    print(f"Pause Codons (<= -15): {len(res.pause_sites)} site(s)")
    return 0


def _cmd_fold(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    if args.melt:
        res = unfold_5prime_synonymously(seq)
        if args.json:
            print(json.dumps(asdict(res), indent=2))
            return 0
        print(f"\nRiboFlow 5' mRNA Hairpin Unfolding — {name}")
        print("=" * 68)
        print(f"Before : dG = {res.before.mfe_kcal_mol:+.2f} kcal/mol | Access = {res.before.initiation_accessibility * 100:.1f}% | {res.before.severity}")
        print(f"         {res.before.dot_bracket}")
        print(f"After  : dG = {res.after.mfe_kcal_mol:+.2f} kcal/mol | Access = {res.after.initiation_accessibility * 100:.1f}% | {res.after.severity}")
        print(f"         {res.after.dot_bracket}")
        print(f"Gain   : ΔMFE = {res.delta_mfe:+.2f} kcal/mol | Access {res.accessibility_gain_percent:+.1f}% | Codons changed: {res.codons_changed} | Protein 100% same: {res.protein_preserved}")
        return 0

    rep = fold_5prime_window(seq, window_nt=args.window)
    if args.json:
        print(json.dumps(asdict(rep), indent=2))
        return 0
    print(f"\nRiboFlow 5' mRNA Secondary Structure — {name}")
    print("=" * 68)
    print(f"5' Window ({rep.window_nt} nt) : {rep.rna_sequence}")
    print(f"Dot-Bracket      : {rep.dot_bracket}")
    print(f"MFE (ΔG_5')      : {rep.mfe_kcal_mol:+.2f} kcal/mol ({rep.severity})")
    print(f"AUG-Zone Paired  : {rep.start_region_paired * 100:.1f}%")
    print(f"Init Access      : {rep.initiation_accessibility * 100:.1f}%")
    return 0


def _cmd_simulate(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    org = get_organism(args.organism)
    sim = simulate_ribosome_traffic(
        seq,
        org,
        initiation_rate=args.initiation_rate,
        duration_seconds=args.duration,
        capture_frames=0,
        couple_5prime_folding=not args.no_5prime_coupling,
        enable_tunnel_drag=not args.no_tunnel_drag,
        seed=args.seed,
    )

    if args.json:
        print(json.dumps(asdict(sim), indent=2))
        return 0

    print(f"\nRiboFlow v2.0 l-TASEP Ribosome Traffic Simulation — {name}")
    print("=" * 66)
    print(f"Organism            : {org.name} (footprint = {sim.ribosome_footprint} codons)")
    print(f"Nominal -> Eff alpha: {sim.initiation_rate:.3f} -> {sim.effective_initiation_rate:.3f} /s (5' ΔG = {sim.folding_5prime.mfe_kcal_mol:+.1f} kcal/mol)")
    print(f"Peak Tunnel Charge  : +{sim.peak_tunnel_charge:.1f} (exit-tunnel drag: {'ON' if sim.enable_tunnel_drag else 'OFF'})")
    print(f"Traffic Phase       : {sim.phase_regime}")
    print(f"Protein Flux (J)    : {sim.protein_flux_per_min:.2f} proteins/min")
    print(f"Mean Ribosomes      : {sim.mean_ribosomes_on_mrna:.2f} on mRNA")
    print(f"Disome Collision %  : {sim.collision_fraction * 100:.2f}% ({sim.collision_events} events)")
    print(f"Max Queue Length    : {sim.max_queue_length} ribosomes")
    print(f"NGD / RQC Risk      : {sim.ngd_risk_score:.1f} / 100")
    print(f"Peak Bottleneck     : codon #{sim.bottleneck_codon + 1}")
    return 0


def _cmd_phase(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    org = get_organism(args.organism)
    rep = sweep_phase_diagram(seq, org, seed=args.seed)
    if args.json:
        print(json.dumps(asdict(rep), indent=2))
        return 0
    print(f"\nRiboFlow TASEP Phase Transition Sweep — {name} ({org.name})")
    print(f"Critical Jamming Threshold alpha_c: {rep.critical_alpha or 'None (jam-free)'} | Max Flux: {rep.max_flux_per_min:.2f}/min")
    print("=" * 74)
    print(f"{'Alpha':>7} {'Eff a':>7} {'Flux/min':>10} {'Density':>9} {'Coll %':>9} {'NGD':>7}  Phase")
    print("-" * 74)
    for pt in rep.points:
        print(
            f"{pt.alpha:>7.2f} {pt.effective_alpha:>7.2f} {pt.protein_flux_per_min:>10.2f} "
            f"{pt.mean_density:>9.3f} {pt.collision_fraction * 100:>8.1f}% {pt.ngd_risk_score:>7.1f}  {pt.phase_regime}"
        )
    return 0


def _cmd_synbio(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    src = get_organism(args.source)
    host = get_organism(args.host)
    rep = design_synbio_sequence(
        seq,
        source=src,
        host=host,
        preserve_pauses=not args.no_pauses,
        unfold_5prime=not args.no_unfold_5prime,
        uridine_deplete=args.vaccine_u_deplete,
        scrub_restriction_sites=not args.no_scrub_sites,
        avoid_homopolymers=not args.no_homopolymer_filter,
        initiation_rate=args.initiation_rate,
        seed=args.seed,
    )
    if args.json:
        print(json.dumps(asdict(rep), indent=2))
        return 0

    print(f"\nRiboFlow v2.0 Multi-Objective SynBio & mRNA Designer — {name}")
    print(f"Transfer: {src.name} -> {host.name} | Protein 100% Preserved: {rep.protein_preserved} | Changed: {rep.percent_changed}%")
    print("=" * 76)
    print(f"5' Hairpin MFE (ΔG) : {rep.before_qc.folding_5prime.mfe_kcal_mol:+.2f} -> {rep.after_qc.folding_5prime.mfe_kcal_mol:+.2f} kcal/mol")
    print(f"5' Init Access      : {rep.before_qc.folding_5prime.initiation_accessibility * 100:.1f}% -> {rep.after_qc.folding_5prime.initiation_accessibility * 100:.1f}%")
    print(f"Folding Rhythm (r)  : {rep.original_rhythm_r:+.3f} -> {rep.designed_rhythm_r:+.3f}")
    print(f"Uridine / UU Pairs  : {rep.before_qc.uridine_percent:.1f}% ({rep.before_qc.uu_dinucleotides} UU) -> {rep.after_qc.uridine_percent:.1f}% ({rep.after_qc.uu_dinucleotides} UU)")
    print(f"Cloning Sites Left  : {len(rep.before_qc.restriction_hits)} -> {len(rep.after_qc.restriction_hits)}")
    print(f"Homopolymer >=5 nt  : {rep.before_qc.homopolymer_violations} -> {rep.after_qc.homopolymer_violations}")
    print(f"Protein Flux (J)    : {rep.before_sim.protein_flux_per_min:.2f} -> {rep.after_sim.protein_flux_per_min:.2f} /min")
    print(f"NGD Collision Risk  : {rep.before_sim.ngd_risk_score:.1f} -> {rep.after_sim.ngd_risk_score:.1f} / 100")

    if args.output_fasta:
        out = Path(args.output_fasta)
        out.write_text(f">{name}_riboflow_v2_synbio\n{rep.designed_sequence}\n", encoding="utf-8")
        print(f"\nSaved SynBio FASTA to: {out}")
    return 0


def _cmd_harmonize(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    src = get_organism(args.source)
    host = get_organism(args.host)
    rep = compare_harmonization(
        seq,
        source=src,
        host=host,
        initiation_rate=args.initiation_rate,
        duration_seconds=args.duration,
        seed=args.seed,
    )

    if args.json:
        print(json.dumps(asdict(rep), indent=2))
        return 0

    print(f"\nRiboFlow Pause-Preserving Codon Harmonization — {name}")
    print(f"Transfer: {src.name} -> {host.name} | Protein 100% preserved: {rep.protein_preserved}")
    print("=" * 78)
    print(
        f"{'Variant':<24} {'CAI':>6} {'Rhythm r':>9} {'Flux/min':>9} {'Coll %':>8} {'NGD Risk':>9}"
    )
    print("-" * 78)
    for v in (rep.wildtype_in_host, rep.naive_optimized, rep.harmonized):
        print(
            f"{v.label:<24} {v.cai_in_host:>6.3f} {v.rhythm_correlation:>+9.3f} "
            f"{v.simulation.protein_flux_per_min:>9.2f} "
            f"{v.simulation.collision_fraction * 100:>7.2f}% "
            f"{v.simulation.ngd_risk_score:>8.1f}"
        )
    if args.output_fasta:
        out = Path(args.output_fasta)
        out.write_text(
            f">{name}_harmonized_{src.key}_to_{host.key}\n{rep.harmonized.sequence}\n",
            encoding="utf-8",
        )
        print(f"\nSaved harmonized FASTA to: {out}")
    return 0


def _cmd_mutate(args: argparse.Namespace) -> int:
    name, seq = _load_sequence(args)
    org = get_organism(args.organism)

    if args.scan:
        hotspots = scan_critical_silent_mutations(
            seq,
            org,
            initiation_rate=args.initiation_rate,
            top_k=args.top_k,
            seed=args.seed,
        )
        if args.json:
            print(json.dumps([asdict(h) for h in hotspots], indent=2))
            return 0

        print(f"\nRiboFlow Critical Synonymous ('Silent') Mutation Hotspots — {name} ({org.name})")
        print("=" * 78)
        print(
            f"{'Codon #':>8} {'AA':>5} {'Swap':>10} {'Rate Drop':>10} {'dFlux':>9} {'dNGD Risk':>10}"
        )
        print("-" * 78)
        for h in hotspots:
            print(
                f"{h.index + 1:>8} {h.aa_three:>5} {h.from_codon + '->' + h.to_codon:>10} "
                f"{-h.rate_drop_percent:>+9.1f}% {h.delta_flux_per_min:>+9.2f} "
                f"{h.delta_ngd_risk:>+10.1f}"
            )
        return 0

    if not args.edit:
        raise ValueError("Provide --scan or at least one --edit INDEX:CODON (1-based index).")

    edits: Dict[int, str] = {}
    for spec in args.edit:
        if ":" not in spec:
            raise ValueError(f"Invalid edit '{spec}'. Expected format: 1_BASED_INDEX:CODON (e.g. 42:AGG)")
        idx_s, codon = spec.split(":", 1)
        edits[int(idx_s) - 1] = codon.strip().upper()

    comp = apply_synonymous_mutations(
        seq,
        org,
        edits=edits,
        initiation_rate=args.initiation_rate,
        duration_seconds=args.duration,
        seed=args.seed,
    )
    if args.json:
        print(json.dumps(asdict(comp), indent=2))
        return 0

    print(f"\nRiboFlow Synonymous Mutation Impact — {name} ({org.name})")
    print(f"Protein 100% Identical: {comp.protein_identical}")
    print("=" * 64)
    for e in comp.edits:
        print(
            f"  Codon #{e.index + 1} ({e.aa_three}): {e.from_codon} -> {e.to_codon} "
            f"({e.from_rate:.2f} -> {e.to_rate:.2f} codons/s)"
        )
    print("-" * 64)
    print(
        f"Protein Flux  : {comp.wildtype_sim.protein_flux_per_min:.2f} -> "
        f"{comp.mutant_sim.protein_flux_per_min:.2f} /min ({comp.delta_flux_percent:+.1f}%)"
    )
    print(
        f"Collision %   : {comp.wildtype_sim.collision_fraction * 100:.2f}% -> "
        f"{comp.mutant_sim.collision_fraction * 100:.2f}%"
    )
    print(
        f"NGD Risk      : {comp.wildtype_sim.ngd_risk_score:.1f} -> "
        f"{comp.mutant_sim.ngd_risk_score:.1f} ({comp.delta_ngd_risk:+.1f})"
    )
    return 0


def _cmd_serve(args: argparse.Namespace) -> int:
    from .web.server import run_server

    run_server(host=args.host, port=args.port, open_browser=args.open)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="riboflow",
        description="RiboFlow v2.0: Biophysical ribosome traffic, 5' mRNA folding & SynBio/mRNA vaccine designer.",
    )
    parser.add_argument("--version", action="version", version=f"riboflow {__version__}")

    sub = parser.add_subparsers(dest="command", required=True)
    org_choices = sorted(ORGANISMS.keys())

    # analyze
    p_an = sub.add_parser("analyze", help="Compute CAI, GC3, RSCU, and %%MinMax pause profile.")
    _add_sequence_args(p_an)
    p_an.add_argument("-o", "--organism", choices=org_choices, default="ecoli")
    p_an.add_argument("--json", action="store_true")
    p_an.set_defaults(func=_cmd_analyze)

    # fold
    p_fold = sub.add_parser("fold", help="Compute 5' mRNA hairpin MFE (dG) or unfold synonymously.")
    _add_sequence_args(p_fold)
    p_fold.add_argument("-w", "--window", type=int, default=45)
    p_fold.add_argument("--melt", action="store_true", help="Melt inhibitory 5' hairpin synonymously.")
    p_fold.add_argument("--json", action="store_true")
    p_fold.set_defaults(func=_cmd_fold)

    # simulate
    p_sim = sub.add_parser("simulate", help="Run l-TASEP simulation with 5' hairpin & exit-tunnel physics.")
    _add_sequence_args(p_sim)
    p_sim.add_argument("-o", "--organism", choices=org_choices, default="ecoli")
    p_sim.add_argument("-a", "--initiation-rate", type=float, default=None)
    p_sim.add_argument("-d", "--duration", type=float, default=120.0)
    p_sim.add_argument("--no-5prime-coupling", action="store_true")
    p_sim.add_argument("--no-tunnel-drag", action="store_true")
    p_sim.add_argument("--seed", type=int, default=42)
    p_sim.add_argument("--json", action="store_true")
    p_sim.set_defaults(func=_cmd_simulate)

    # phase
    p_ph = sub.add_parser("phase", help="Sweep initiation rate alpha to map TASEP phase transition.")
    _add_sequence_args(p_ph)
    p_ph.add_argument("-o", "--organism", choices=org_choices, default="ecoli")
    p_ph.add_argument("--seed", type=int, default=42)
    p_ph.add_argument("--json", action="store_true")
    p_ph.set_defaults(func=_cmd_phase)

    # synbio
    p_sb = sub.add_parser("synbio", help="Multi-objective SynBio & mRNA vaccine (U-depletion + site scrubbing) designer.")
    _add_sequence_args(p_sb)
    p_sb.add_argument("-s", "--source", choices=org_choices, default="human")
    p_sb.add_argument("-t", "--host", choices=org_choices, default="ecoli")
    p_sb.add_argument("--vaccine-u-deplete", action="store_true", help="Enable mRNA vaccine uridine depletion.")
    p_sb.add_argument("--no-pauses", action="store_true")
    p_sb.add_argument("--no-unfold-5prime", action="store_true")
    p_sb.add_argument("--no-scrub-sites", action="store_true")
    p_sb.add_argument("--no-homopolymer-filter", action="store_true")
    p_sb.add_argument("-a", "--initiation-rate", type=float, default=None)
    p_sb.add_argument("--seed", type=int, default=42)
    p_sb.add_argument("--output-fasta", default=None)
    p_sb.add_argument("--json", action="store_true")
    p_sb.set_defaults(func=_cmd_synbio)

    # harmonize
    p_harm = sub.add_parser("harmonize", help="Compare native, naive CAI, and pause-preserving harmonization.")
    _add_sequence_args(p_harm)
    p_harm.add_argument("-s", "--source", choices=org_choices, default="human")
    p_harm.add_argument("-t", "--host", choices=org_choices, default="ecoli")
    p_harm.add_argument("-a", "--initiation-rate", type=float, default=None)
    p_harm.add_argument("-d", "--duration", type=float, default=100.0)
    p_harm.add_argument("--seed", type=int, default=42)
    p_harm.add_argument("--output-fasta", default=None)
    p_harm.add_argument("--json", action="store_true")
    p_harm.set_defaults(func=_cmd_harmonize)

    # mutate
    p_mut = sub.add_parser("mutate", help="Test synonymous ('silent') mutations or scan critical hotspots.")
    _add_sequence_args(p_mut)
    p_mut.add_argument("-o", "--organism", choices=org_choices, default="ecoli")
    p_mut.add_argument("--edit", action="append", help="1-based silent edit INDEX:CODON (repeatable).")
    p_mut.add_argument("--scan", action="store_true", help="Scan for highest-risk single-nt silent mutations.")
    p_mut.add_argument("--top-k", type=int, default=6)
    p_mut.add_argument("-a", "--initiation-rate", type=float, default=None)
    p_mut.add_argument("-d", "--duration", type=float, default=100.0)
    p_mut.add_argument("--seed", type=int, default=42)
    p_mut.add_argument("--json", action="store_true")
    p_mut.set_defaults(func=_cmd_mutate)

    # serve
    p_srv = sub.add_parser("serve", help="Launch the interactive RiboFlow v2.0 web laboratory.")
    p_srv.add_argument("--host", default="127.0.0.1")
    p_srv.add_argument("--port", type=int, default=8765)
    p_srv.add_argument("--open", action="store_true", help="Open web browser automatically.")
    p_srv.set_defaults(func=_cmd_serve)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
