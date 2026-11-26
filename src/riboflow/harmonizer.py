"""Pause-preserving codon harmonization vs. naive CAI maximization.

Naive codon optimization replaces every codon with the host organism's most
frequent synonymous codon. While this maximizes CAI, it flattens evolutionarily
conserved translational pause sites required for co-translational protein folding,
often leading to misfolding and inclusion bodies.

Codon harmonization (Angov et al., 2008; Clarke & Clark %MinMax rhythm matching)
maps each source codon's relative usage rank to the closest equivalent rank in
the target host, preserving folding pauses while clamping extreme dead-stop
codons that trigger ribosome collisions and No-Go Decay.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

from .organisms import GENETIC_CODE, SYNONYMOUS_CODONS, OrganismProfile
from .sequence import (
    analyze_sequence,
    compute_cai,
    compute_percent_minmax,
    translate_codons,
    validate_cds,
)
from .simulator import TrafficSimulationResult, simulate_ribosome_traffic


@dataclass(frozen=True)
class VariantSummary:
    label: str
    sequence: str
    cai_in_host: float
    gc_percent: float
    rhythm_correlation: float     # Pearson r vs native %MinMax profile in source organism
    codons_changed: int
    percent_changed: float
    minmax_profile: List[float]
    simulation: TrafficSimulationResult


@dataclass(frozen=True)
class HarmonizationReport:
    source_organism: str
    host_organism: str
    protein_sequence: str
    protein_preserved: bool
    native_source_minmax: List[float]
    wildtype_in_host: VariantSummary
    naive_optimized: VariantSummary
    harmonized: VariantSummary


def _pearson_r(xs: List[float], ys: List[float]) -> float:
    n = min(len(xs), len(ys))
    if n < 2:
        return 0.0
    mx = sum(xs[:n]) / n
    my = sum(ys[:n]) / n
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    dx = math.sqrt(sum((xs[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ys[i] - my) ** 2 for i in range(n)))
    if dx < 1e-9 or dy < 1e-9:
        return 0.0
    return round(num / (dx * dy), 4)


def naive_optimize_codons(codons: List[str], host: OrganismProfile) -> List[str]:
    """Replace each sense codon with the highest-frequency host codon."""
    freq = host.codon_freq_per_1000
    out: List[str] = []
    for idx, c in enumerate(codons):
        aa = GENETIC_CODE[c]
        if idx == 0 or aa == "*":
            out.append(c)
            continue
        best = max(SYNONYMOUS_CODONS[aa], key=lambda syn: freq[syn])
        out.append(best)
    return out


def harmonize_codons(
    codons: List[str],
    source: OrganismProfile,
    host: OrganismProfile,
    min_safe_freq: float = 4.5,
) -> List[str]:
    """Map source codon usage rhythm to host while avoiding ultra-rare dead-stops.

    For each sense codon:
      1. Look up its normalized intra-synonymous usage rank in `source` (0=rarest, 1=most common).
      2. Filter host synonymous codons to avoid ultra-rare dead-stops (freq < min_safe_freq)
         whenever safer synonymous choices exist, preventing frameshift/NGD stalls.
      3. Select the host synonymous codon whose normalized rank best matches the source rank.
    """
    src_ranks = source.normalized_usage_ranks()
    host_ranks = host.normalized_usage_ranks()
    host_freq = host.codon_freq_per_1000

    out: List[str] = []
    for idx, c in enumerate(codons):
        aa = GENETIC_CODE[c]
        if idx == 0 or aa == "*":
            out.append(c)
            continue

        candidates = SYNONYMOUS_CODONS[aa]
        if len(candidates) == 1:
            out.append(c)
            continue

        safe_candidates = [syn for syn in candidates if host_freq[syn] >= min_safe_freq]
        pool = safe_candidates if safe_candidates else candidates

        target_rank = src_ranks[c]
        best = min(
            pool,
            key=lambda syn: (abs(host_ranks[syn] - target_rank), -host_freq[syn]),
        )
        out.append(best)
    return out


def compare_harmonization(
    sequence: str,
    source: OrganismProfile,
    host: OrganismProfile,
    initiation_rate: Optional[float] = None,
    duration_seconds: float = 100.0,
    seed: int = 42,
) -> HarmonizationReport:
    """Compare Wild-Type, Naive CAI Optimization, and Pause-Preserving Harmonization."""
    wt_codons = validate_cds(sequence)
    protein = translate_codons(wt_codons)

    native_source_minmax = compute_percent_minmax(wt_codons, source)

    naive_codons = naive_optimize_codons(wt_codons, host)
    harm_codons = harmonize_codons(wt_codons, source, host)

    protein_preserved = (
        translate_codons(naive_codons) == protein
        and translate_codons(harm_codons) == protein
    )

    def build_variant(label: str, cand_codons: List[str]) -> VariantSummary:
        cand_seq = "".join(cand_codons)
        analysis = analyze_sequence(cand_seq, host)
        changed = sum(1 for a, b in zip(wt_codons, cand_codons) if a != b)
        pct_changed = round((changed / len(wt_codons)) * 100.0, 1)
        mm = compute_percent_minmax(cand_codons, host)
        sim = simulate_ribosome_traffic(
            cand_seq,
            host,
            initiation_rate=initiation_rate,
            duration_seconds=duration_seconds,
            capture_frames=80,
            seed=seed,
        )
        return VariantSummary(
            label=label,
            sequence=cand_seq,
            cai_in_host=compute_cai(cand_codons, host),
            gc_percent=analysis.gc_percent,
            rhythm_correlation=_pearson_r(native_source_minmax, mm),
            codons_changed=changed,
            percent_changed=pct_changed,
            minmax_profile=mm,
            simulation=sim,
        )

    return HarmonizationReport(
        source_organism=source.key,
        host_organism=host.key,
        protein_sequence=protein,
        protein_preserved=protein_preserved,
        native_source_minmax=native_source_minmax,
        wildtype_in_host=build_variant("Wild-Type (Unmodified)", wt_codons),
        naive_optimized=build_variant("Naive CAI Optimized", naive_codons),
        harmonized=build_variant("RiboFlow Harmonized", harm_codons),
    )
