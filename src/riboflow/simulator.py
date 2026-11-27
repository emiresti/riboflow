"""Exact continuous-time l-TASEP stochastic ribosome traffic simulator (v2.0).

Models ribosomes of footprint length `l` (default 10 codons = ~30 nt) translating
an mRNA transcript of length L codons using a Gillespie / kinetic Monte Carlo
Totally Asymmetric Simple Exclusion Process (l-TASEP), now coupled with:

  1. 5' mRNA Secondary Structure Accessibility: stable hairpins near AUG
     attenuate initiation rate via alpha_eff = alpha * accessibility(ΔG_5').
  2. Ribosome Exit-Tunnel Electrostatic Drag (Charneski & Hurst, PLOS Biol 2013):
     clusters of positively charged nascent residues (K, R, H) inside the
     negatively charged ribosomal exit tunnel slow down local elongation.
  3. Statistical Physics Phase Diagram Sweeper: sweeps initiation rate alpha
     to reveal the Low-Density -> Maximum-Current -> Traffic-Jam phase transition.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import List, Optional, Tuple

from .folding import FivePrimeFoldingReport, fold_5prime_window
from .organisms import GENETIC_CODE, OrganismProfile
from .sequence import validate_cds


@dataclass(frozen=True)
class FrameRibosome:
    id: int
    pos: int
    colliding: bool
    tunnel_stalled: bool = False


@dataclass(frozen=True)
class SimulationFrame:
    time: float
    completed: int
    ribosomes: List[FrameRibosome]


@dataclass(frozen=True)
class TrafficSimulationResult:
    organism: str
    codon_count: int
    initiation_rate: float            # nominal alpha requested
    effective_initiation_rate: float  # alpha_eff after 5' hairpin coupling
    termination_rate: float
    ribosome_footprint: int
    duration_seconds: float
    burn_in_seconds: float
    couple_5prime_folding: bool
    enable_tunnel_drag: bool
    folding_5prime: FivePrimeFoldingReport
    completed_proteins: int
    protein_flux_per_min: float
    mean_ribosomes_on_mrna: float
    mean_density: float
    collision_fraction: float
    collision_events: int
    max_queue_length: int
    ngd_risk_score: float
    phase_regime: str
    bottleneck_codon: int
    peak_tunnel_charge: float
    codon_occupancy: List[float]
    collision_profile: List[float]
    codon_rates: List[float]          # effective decoding rates after tunnel drag
    tunnel_charge_profile: List[float]
    tunnel_drag_factors: List[float]
    frames: List[SimulationFrame]


@dataclass(frozen=True)
class PhasePoint:
    alpha: float
    effective_alpha: float
    protein_flux_per_min: float
    mean_density: float
    collision_fraction: float
    ngd_risk_score: float
    phase_regime: str


@dataclass(frozen=True)
class PhaseDiagramReport:
    organism: str
    critical_alpha: Optional[float]  # alpha_c where collision_fraction crosses 10%
    max_flux_per_min: float
    points: List[PhasePoint]


def compute_exit_tunnel_profile(
    codons: List[str],
    tunnel_window: int = 12,
) -> Tuple[List[float], List[float]]:
    """Compute nascent-peptide positive charge inside the ribosomal exit tunnel.

    Returns (charge_profile, drag_factors) per codon:
      - Positive residues K (+1.0), R (+1.0), H (+0.35) create electrostatic drag
        against the negatively charged rRNA tunnel wall.
      - Negative residues D (-0.4), E (-0.4) partially screen positive charge.
      - When net positive charge in the preceding 12-aa tunnel window exceeds +3.0,
        elongation rate is attenuated smoothly down to 0.28x.
    """
    aa_charge = {"K": 1.0, "R": 1.0, "H": 0.35, "D": -0.4, "E": -0.4}
    aas = [GENETIC_CODE[c] for c in codons]
    charges: List[float] = []
    drags: List[float] = []

    for i in range(len(aas)):
        lo = max(0, i - tunnel_window)
        window = aas[lo:i]  # nascent chain already extruded into exit tunnel
        net_pos = sum(aa_charge.get(a, 0.0) for a in window)
        net_pos = round(max(0.0, net_pos), 2)
        charges.append(net_pos)

        if net_pos <= 3.0:
            drags.append(1.0)
        else:
            excess = net_pos - 3.0
            factor = max(0.26, math.exp(-0.24 * excess))
            drags.append(round(factor, 3))

    return charges, drags


def _classify_phase(
    mean_density: float,
    collision_fraction: float,
    effective_alpha: float,
    min_rate: float,
) -> str:
    if collision_fraction >= 0.14 or mean_density >= 0.45:
        return "High-Density (Traffic Jam)"
    if collision_fraction >= 0.04 or effective_alpha >= min_rate * 0.45:
        return "Maximum-Current (Near Capacity)"
    return "Low-Density (Initiation-Limited)"


def _compute_ngd_risk(
    collision_fraction: float,
    max_queue_length: int,
    peak_collision_density: float,
) -> float:
    score = (
        collision_fraction * 145.0
        + max(0, max_queue_length - 1) * 9.5
        + peak_collision_density * 65.0
    )
    return round(min(100.0, max(0.0, score)), 1)


def simulate_ribosome_traffic(
    sequence: str,
    organism: OrganismProfile,
    initiation_rate: Optional[float] = None,
    termination_rate: Optional[float] = None,
    duration_seconds: float = 120.0,
    burn_in_seconds: float = 30.0,
    capture_frames: int = 160,
    couple_5prime_folding: bool = True,
    enable_tunnel_drag: bool = True,
    seed: int = 42,
) -> TrafficSimulationResult:
    """Run exact continuous-time l-TASEP ribosome traffic simulation."""
    codons = validate_cds(sequence)
    dna = "".join(codons)
    L = len(codons)
    footprint = organism.ribosome_footprint

    nominal_alpha = (
        initiation_rate
        if initiation_rate is not None
        else organism.default_initiation_rate
    )
    if nominal_alpha <= 0:
        raise ValueError("initiation_rate must be positive.")

    folding_rep = fold_5prime_window(dna, window_nt=45)
    alpha_eff = (
        nominal_alpha * folding_rep.initiation_accessibility
        if couple_5prime_folding
        else nominal_alpha
    )

    rate_map = organism.elongation_rates()
    base_rates = [rate_map[c] for c in codons]
    tunnel_charges, tunnel_drags = compute_exit_tunnel_profile(codons)

    if enable_tunnel_drag:
        rates = [round(r * d, 3) for r, d in zip(base_rates, tunnel_drags)]
    else:
        rates = base_rates
        tunnel_drags = [1.0] * L

    beta = (
        termination_rate
        if termination_rate is not None
        else max(rates[-1], organism.mean_elongation_rate * 0.8)
    )
    if beta <= 0:
        raise ValueError("termination_rate must be positive.")

    rng = random.Random(seed)
    positions: List[int] = []
    ribo_ids: List[int] = []
    next_ribo_id = 1

    occ_time = [0.0] * L
    coll_time = [0.0] * L

    total_time = burn_in_seconds + duration_seconds
    t = 0.0
    completed = 0
    total_ribo_time = 0.0
    total_coll_ribo_time = 0.0
    collision_events = 0
    max_queue_length = 0

    frames: List[SimulationFrame] = []
    frame_interval = (
        duration_seconds / capture_frames if capture_frames > 0 else float("inf")
    )
    next_frame_t = burn_in_seconds

    def record_frame(current_t: float) -> None:
        n = len(positions)
        snap: List[FrameRibosome] = []
        for idx in range(n):
            pos = positions[idx]
            is_coll = (idx + 1 < n) and (positions[idx + 1] - pos == footprint)
            is_tunnel = enable_tunnel_drag and tunnel_drags[pos] < 0.78
            snap.append(
                FrameRibosome(
                    id=ribo_ids[idx],
                    pos=pos,
                    colliding=is_coll,
                    tunnel_stalled=is_tunnel,
                )
            )
        frames.append(
            SimulationFrame(
                time=round(current_t - burn_in_seconds, 2),
                completed=completed,
                ribosomes=snap,
            )
        )

    while t < total_time:
        n = len(positions)
        can_initiate = (n == 0) or (positions[0] >= footprint)
        total_rate = alpha_eff if can_initiate else 0.0

        move_rates: List[float] = [0.0] * n
        for idx in range(n):
            pos = positions[idx]
            if pos == L - 1:
                r = beta
            elif idx + 1 < n and positions[idx + 1] - pos <= footprint:
                r = 0.0
            else:
                r = rates[pos]
            move_rates[idx] = r
            total_rate += r

        if total_rate <= 0.0:
            break

        u = rng.random()
        dt = -math.log(max(u, 1e-12)) / total_rate
        t_next = min(t + dt, total_time)

        if t_next > burn_in_seconds:
            seg_start = max(t, burn_in_seconds)
            seg_dt = t_next - seg_start
            if seg_dt > 0.0 and n > 0:
                total_ribo_time += n * seg_dt
                current_queue = 0
                for idx in range(n):
                    pos = positions[idx]
                    occ_time[pos] += seg_dt
                    is_coll = (idx + 1 < n) and (positions[idx + 1] - pos == footprint)
                    if is_coll:
                        coll_time[pos] += seg_dt
                        total_coll_ribo_time += seg_dt
                        current_queue += 1
                        if current_queue + 1 > max_queue_length:
                            max_queue_length = current_queue + 1
                    else:
                        current_queue = 0

            while next_frame_t <= t_next and len(frames) < capture_frames:
                record_frame(next_frame_t)
                next_frame_t += frame_interval

        t = t_next
        if t >= total_time:
            break

        threshold = rng.random() * total_rate
        if can_initiate:
            if threshold < alpha_eff:
                positions.insert(0, 0)
                ribo_ids.insert(0, next_ribo_id)
                next_ribo_id += 1
                if t >= burn_in_seconds and len(positions) > 1 and positions[1] == footprint:
                    collision_events += 1
                continue
            threshold -= alpha_eff

        chosen = n - 1
        acc = 0.0
        for idx, r in enumerate(move_rates):
            if r <= 0.0:
                continue
            acc += r
            if threshold <= acc:
                chosen = idx
                break

        if positions[chosen] == L - 1:
            positions.pop(chosen)
            ribo_ids.pop(chosen)
            if t >= burn_in_seconds:
                completed += 1
        else:
            positions[chosen] += 1
            if (
                t >= burn_in_seconds
                and chosen + 1 < len(positions)
                and positions[chosen + 1] - positions[chosen] == footprint
            ):
                collision_events += 1

    meas = max(duration_seconds, 1e-9)
    codon_occupancy = [round(x / meas, 4) for x in occ_time]
    collision_profile = [round(x / meas, 4) for x in coll_time]

    mean_ribos = total_ribo_time / meas
    max_capacity = max(1.0, L / footprint)
    mean_density = min(1.0, mean_ribos / max_capacity)
    collision_fraction = (
        total_coll_ribo_time / total_ribo_time if total_ribo_time > 0 else 0.0
    )
    protein_flux_per_min = round((completed / meas) * 60.0, 2)

    peak_coll_density = max(collision_profile) if collision_profile else 0.0
    ngd_risk = _compute_ngd_risk(collision_fraction, max_queue_length, peak_coll_density)
    phase = _classify_phase(mean_density, collision_fraction, alpha_eff, min(rates))

    if peak_coll_density > 0.0:
        bottleneck = max(range(L), key=lambda i: (collision_profile[i], codon_occupancy[i]))
    else:
        bottleneck = max(range(L), key=lambda i: codon_occupancy[i])

    return TrafficSimulationResult(
        organism=organism.key,
        codon_count=L,
        initiation_rate=round(nominal_alpha, 3),
        effective_initiation_rate=round(alpha_eff, 3),
        termination_rate=round(beta, 3),
        ribosome_footprint=footprint,
        duration_seconds=duration_seconds,
        burn_in_seconds=burn_in_seconds,
        couple_5prime_folding=couple_5prime_folding,
        enable_tunnel_drag=enable_tunnel_drag,
        folding_5prime=folding_rep,
        completed_proteins=completed,
        protein_flux_per_min=protein_flux_per_min,
        mean_ribosomes_on_mrna=round(mean_ribos, 2),
        mean_density=round(mean_density, 3),
        collision_fraction=round(collision_fraction, 4),
        collision_events=collision_events,
        max_queue_length=max_queue_length,
        ngd_risk_score=ngd_risk,
        phase_regime=phase,
        bottleneck_codon=bottleneck,
        peak_tunnel_charge=max(tunnel_charges) if tunnel_charges else 0.0,
        codon_occupancy=codon_occupancy,
        collision_profile=collision_profile,
        codon_rates=rates,
        tunnel_charge_profile=tunnel_charges,
        tunnel_drag_factors=tunnel_drags,
        frames=frames,
    )


def sweep_phase_diagram(
    sequence: str,
    organism: OrganismProfile,
    alpha_values: Optional[List[float]] = None,
    couple_5prime_folding: bool = True,
    enable_tunnel_drag: bool = True,
    duration_seconds: float = 65.0,
    seed: int = 42,
) -> PhaseDiagramReport:
    """Sweep initiation rate alpha to map the TASEP statistical-physics phase diagram."""
    alphas = alpha_values or [0.08, 0.16, 0.28, 0.42, 0.60, 0.85, 1.10, 1.40, 1.75, 2.20]
    points: List[PhasePoint] = []
    critical_alpha: Optional[float] = None

    for a in alphas:
        sim = simulate_ribosome_traffic(
            sequence,
            organism,
            initiation_rate=a,
            duration_seconds=duration_seconds,
            burn_in_seconds=20.0,
            capture_frames=0,
            couple_5prime_folding=couple_5prime_folding,
            enable_tunnel_drag=enable_tunnel_drag,
            seed=seed,
        )
        if critical_alpha is None and sim.collision_fraction >= 0.10:
            critical_alpha = round(a, 3)
        points.append(
            PhasePoint(
                alpha=round(a, 3),
                effective_alpha=sim.effective_initiation_rate,
                protein_flux_per_min=sim.protein_flux_per_min,
                mean_density=sim.mean_density,
                collision_fraction=sim.collision_fraction,
                ngd_risk_score=sim.ngd_risk_score,
                phase_regime=sim.phase_regime,
            )
        )

    max_flux = max((p.protein_flux_per_min for p in points), default=0.0)
    return PhaseDiagramReport(
        organism=organism.key,
        critical_alpha=critical_alpha,
        max_flux_per_min=max_flux,
        points=points,
    )
