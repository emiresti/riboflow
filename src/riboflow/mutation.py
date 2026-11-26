"""Synonymous ('silent') mutation kinetic laboratory and hotspot scanner.

Demonstrates how single-nucleotide synonymous substitutions—which leave the
encoded amino acid sequence 100% unchanged—can dramatically alter local
elongation speed, create ribosome bottlenecks, and trigger disome collisions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .organisms import AMINO_ACID_NAMES, GENETIC_CODE, SYNONYMOUS_CODONS, OrganismProfile
from .sequence import compute_cai, compute_percent_minmax, translate_codons, validate_cds
from .simulator import TrafficSimulationResult, simulate_ribosome_traffic


@dataclass(frozen=True)
class SilentEdit:
    index: int
    amino_acid: str
    aa_three: str
    from_codon: str
    to_codon: str
    from_rate: float
    to_rate: float
    rate_ratio: float


@dataclass(frozen=True)
class MutationComparison:
    organism: str
    edits: List[SilentEdit]
    protein_identical: bool
    wildtype_sequence: str
    mutant_sequence: str
    wildtype_cai: float
    mutant_cai: float
    delta_flux_per_min: float
    delta_flux_percent: float
    delta_collision_fraction: float
    delta_ngd_risk: float
    wildtype_minmax: List[float]
    mutant_minmax: List[float]
    wildtype_sim: TrafficSimulationResult
    mutant_sim: TrafficSimulationResult


@dataclass(frozen=True)
class CriticalSilentHotspot:
    index: int
    amino_acid: str
    aa_three: str
    from_codon: str
    to_codon: str
    nucleotide_changes: int
    rate_drop_percent: float
    delta_flux_per_min: float
    delta_ngd_risk: float
    mutant_ngd_risk: float
    mutant_collision_fraction: float


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def apply_synonymous_mutations(
    sequence: str,
    organism: OrganismProfile,
    edits: Dict[int, str],
    initiation_rate: Optional[float] = None,
    duration_seconds: float = 120.0,
    seed: int = 42,
) -> MutationComparison:
    """Apply synonymous codon edits and compare WT vs. Silent Mutant kinetics."""
    wt_codons = validate_cds(sequence)
    mut_codons = list(wt_codons)
    rates = organism.elongation_rates()

    applied_edits: List[SilentEdit] = []
    for idx, new_codon_raw in sorted(edits.items()):
        if idx <= 0 or idx >= len(wt_codons) - 1:
            raise ValueError(
                f"Codon index {idx} cannot be mutated (start and terminal stop codons are locked)."
            )
        new_codon = new_codon_raw.strip().upper().replace("U", "T")
        if new_codon not in GENETIC_CODE:
            raise ValueError(f"Invalid codon '{new_codon_raw}'.")

        old_codon = wt_codons[idx]
        old_aa = GENETIC_CODE[old_codon]
        new_aa = GENETIC_CODE[new_codon]
        if old_aa != new_aa:
            raise ValueError(
                f"Edit at codon {idx + 1} ({old_codon}->{new_codon}) is NOT synonymous: "
                f"{old_aa} != {new_aa}."
            )

        mut_codons[idx] = new_codon
        if new_codon != old_codon:
            r_old = rates[old_codon]
            r_new = rates[new_codon]
            applied_edits.append(
                SilentEdit(
                    index=idx,
                    amino_acid=old_aa,
                    aa_three=AMINO_ACID_NAMES[old_aa],
                    from_codon=old_codon,
                    to_codon=new_codon,
                    from_rate=r_old,
                    to_rate=r_new,
                    rate_ratio=round(r_new / r_old, 3) if r_old > 0 else 1.0,
                )
            )

    wt_seq = "".join(wt_codons)
    mut_seq = "".join(mut_codons)
    protein_identical = translate_codons(wt_codons) == translate_codons(mut_codons)

    wt_sim = simulate_ribosome_traffic(
        wt_seq,
        organism,
        initiation_rate=initiation_rate,
        duration_seconds=duration_seconds,
        capture_frames=120,
        seed=seed,
    )
    mut_sim = simulate_ribosome_traffic(
        mut_seq,
        organism,
        initiation_rate=initiation_rate,
        duration_seconds=duration_seconds,
        capture_frames=120,
        seed=seed,
    )

    delta_flux = round(mut_sim.protein_flux_per_min - wt_sim.protein_flux_per_min, 2)
    delta_flux_pct = (
        round((delta_flux / wt_sim.protein_flux_per_min) * 100.0, 1)
        if wt_sim.protein_flux_per_min > 0
        else 0.0
    )

    return MutationComparison(
        organism=organism.key,
        edits=applied_edits,
        protein_identical=protein_identical,
        wildtype_sequence=wt_seq,
        mutant_sequence=mut_seq,
        wildtype_cai=compute_cai(wt_codons, organism),
        mutant_cai=compute_cai(mut_codons, organism),
        delta_flux_per_min=delta_flux,
        delta_flux_percent=delta_flux_pct,
        delta_collision_fraction=round(
            mut_sim.collision_fraction - wt_sim.collision_fraction, 4
        ),
        delta_ngd_risk=round(mut_sim.ngd_risk_score - wt_sim.ngd_risk_score, 1),
        wildtype_minmax=compute_percent_minmax(wt_codons, organism),
        mutant_minmax=compute_percent_minmax(mut_codons, organism),
        wildtype_sim=wt_sim,
        mutant_sim=mut_sim,
    )


def scan_critical_silent_mutations(
    sequence: str,
    organism: OrganismProfile,
    initiation_rate: Optional[float] = None,
    top_k: int = 6,
    single_nt_only: bool = True,
    seed: int = 42,
) -> List[CriticalSilentHotspot]:
    """Discover synonymous mutations that cause the largest ribosome collision spikes.

    Uses a two-stage screen:
      1. Fast kinetic screen: rank degenerate codons by largest decoding-rate drop.
      2. Full TASEP verification on the top candidate sites to quantify exact
         disome collision / NGD risk increase.
    """
    wt_codons = validate_cds(sequence)
    rates = organism.elongation_rates()

    wt_sim = simulate_ribosome_traffic(
        "".join(wt_codons),
        organism,
        initiation_rate=initiation_rate,
        duration_seconds=80.0,
        burn_in_seconds=20.0,
        capture_frames=0,
        seed=seed,
    )

    # Stage 1: rank candidate silent swaps by relative rate drop and downstream position
    candidates = []
    for idx in range(1, len(wt_codons) - 1):
        old_c = wt_codons[idx]
        aa = GENETIC_CODE[old_c]
        syns = SYNONYMOUS_CODONS[aa]
        if len(syns) <= 1:
            continue

        valid_syns = [
            s
            for s in syns
            if s != old_c and (not single_nt_only or _hamming(old_c, s) == 1)
        ]
        if not valid_syns:
            continue

        slowest = min(valid_syns, key=lambda s: rates[s])
        r_old, r_new = rates[old_c], rates[slowest]
        if r_new >= r_old * 0.85:
            continue

        # Favor codons downstream of the first footprint (idx >= 10) where queues can form
        queue_weight = 1.25 if idx >= organism.ribosome_footprint else 0.75
        drop_ratio = (r_old - r_new) / r_old
        heuristic = drop_ratio * (1.0 / max(r_new, 0.2)) * queue_weight
        candidates.append((heuristic, idx, old_c, slowest, aa, r_old, r_new))

    candidates.sort(key=lambda item: item[0], reverse=True)
    shortlist = candidates[: max(top_k * 3, 14)]

    # Stage 2: run TASEP on shortlisted silent mutants
    results: List[CriticalSilentHotspot] = []
    for _, idx, old_c, new_c, aa, r_old, r_new in shortlist:
        mut_codons = list(wt_codons)
        mut_codons[idx] = new_c
        mut_sim = simulate_ribosome_traffic(
            "".join(mut_codons),
            organism,
            initiation_rate=initiation_rate,
            duration_seconds=80.0,
            burn_in_seconds=20.0,
            capture_frames=0,
            seed=seed,
        )
        results.append(
            CriticalSilentHotspot(
                index=idx,
                amino_acid=aa,
                aa_three=AMINO_ACID_NAMES[aa],
                from_codon=old_c,
                to_codon=new_c,
                nucleotide_changes=_hamming(old_c, new_c),
                rate_drop_percent=round(((r_old - r_new) / r_old) * 100.0, 1),
                delta_flux_per_min=round(
                    mut_sim.protein_flux_per_min - wt_sim.protein_flux_per_min, 2
                ),
                delta_ngd_risk=round(mut_sim.ngd_risk_score - wt_sim.ngd_risk_score, 1),
                mutant_ngd_risk=mut_sim.ngd_risk_score,
                mutant_collision_fraction=mut_sim.collision_fraction,
            )
        )

    results.sort(
        key=lambda r: (r.delta_ngd_risk, -r.delta_flux_per_min, r.rate_drop_percent),
        reverse=True,
    )
    return results[:top_k]
