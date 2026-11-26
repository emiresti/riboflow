"""Exact continuous-time l-TASEP stochastic ribosome traffic simulator.

Models ribosomes of footprint length `l` (default 10 codons = ~30 nt) translating
an mRNA transcript of length L codons using a Gillespie / kinetic Monte Carlo
Totally Asymmetric Simple Exclusion Process (l-TASEP):

  - Initiation at codon 0 with rate alpha when codons [0 .. l-1] are unobstructed.
  - Codon-specific elongation from i -> i+1 with rate k_i when i + l is vacant.
  - Termination at codon L-1 with rate beta.

Tracks simulated Ribo-seq (per-codon A-site occupancy), Disome-seq (steric
collision density), protein production flux J (proteins/min), and No-Go Decay
(NGD) / RQC risk.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Dict, List, Optional

from .organisms import OrganismProfile
from .sequence import validate_cds


@dataclass(frozen=True)
class FrameRibosome:
    id: int
    pos: int
    colliding: bool


@dataclass(frozen=True)
class SimulationFrame:
    time: float
    completed: int
    ribosomes: List[FrameRibosome]


@dataclass(frozen=True)
class TrafficSimulationResult:
    organism: str
    codon_count: int
    initiation_rate: float
    termination_rate: float
    ribosome_footprint: int
    duration_seconds: float
    burn_in_seconds: float
    completed_proteins: int
    protein_flux_per_min: float
    mean_ribosomes_on_mrna: float
    mean_density: float
    collision_fraction: float       # fraction of ribosome-time spent in disome collision
    collision_events: int           # number of newly formed disome collision contacts
    max_queue_length: int           # max contiguous collided ribosome chain observed
    ngd_risk_score: float           # 0..100 No-Go Decay / RQC activation risk
    phase_regime: str               # Low-Density / Maximum-Current / High-Density
    bottleneck_codon: int           # 0-based index of peak occupancy/stall
    codon_occupancy: List[float]    # simulated Ribo-seq density per codon
    collision_profile: List[float]  # simulated Disome-seq density per codon
    codon_rates: List[float]        # decoding rate k_i (codons/s) per codon
    frames: List[SimulationFrame]   # animation snapshots for live visualizer


def _classify_phase(
    mean_density: float,
    collision_fraction: float,
    initiation_rate: float,
    min_rate: float,
) -> str:
    if collision_fraction >= 0.14 or mean_density >= 0.45:
        return "High-Density (Traffic Jam)"
    if collision_fraction >= 0.04 or initiation_rate >= min_rate * 0.45:
        return "Maximum-Current (Near Capacity)"
    return "Low-Density (Initiation-Limited)"


def _compute_ngd_risk(
    collision_fraction: float,
    max_queue_length: int,
    peak_collision_density: float,
) -> float:
    """Estimate No-Go Decay (NGD) / RQC risk score on a 0-100 scale."""
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
    seed: int = 42,
) -> TrafficSimulationResult:
    """Run exact continuous-time l-TASEP ribosome traffic simulation."""
    codons = validate_cds(sequence)
    L = len(codons)
    footprint = organism.ribosome_footprint

    alpha = (
        initiation_rate
        if initiation_rate is not None
        else organism.default_initiation_rate
    )
    if alpha <= 0:
        raise ValueError("initiation_rate must be positive.")

    rate_map = organism.elongation_rates()
    rates = [rate_map[c] for c in codons]
    beta = (
        termination_rate
        if termination_rate is not None
        else max(rates[-1], organism.mean_elongation_rate * 0.8)
    )
    if beta <= 0:
        raise ValueError("termination_rate must be positive.")

    rng = random.Random(seed)

    # Ribosome state kept sorted ascending by codon index: positions[0] < ... < positions[-1]
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

    # Animation frame sampling during measurement phase
    frames: List[SimulationFrame] = []
    frame_interval = (
        duration_seconds / capture_frames if capture_frames > 0 else float("inf")
    )
    next_frame_t = burn_in_seconds

    def record_frame(current_t: float) -> None:
        n = len(positions)
        snap: List[FrameRibosome] = []
        for idx in range(n):
            is_coll = (idx + 1 < n) and (positions[idx + 1] - positions[idx] == footprint)
            snap.append(FrameRibosome(id=ribo_ids[idx], pos=positions[idx], colliding=is_coll))
        frames.append(
            SimulationFrame(
                time=round(current_t - burn_in_seconds, 2),
                completed=completed,
                ribosomes=snap,
            )
        )

    while t < total_time:
        n = len(positions)

        # Build propensity list
        # Event type: -1 = initiation, i in [0..n-1] = step/terminate ribosome i
        can_initiate = (n == 0) or (positions[0] >= footprint)
        total_rate = alpha if can_initiate else 0.0

        move_rates: List[float] = [0.0] * n
        for idx in range(n):
            pos = positions[idx]
            if pos == L - 1:
                r = beta
            elif idx + 1 < n and positions[idx + 1] - pos <= footprint:
                r = 0.0  # sterically blocked (disome collision!)
            else:
                r = rates[pos]
            move_rates[idx] = r
            total_rate += r

        if total_rate <= 0.0:
            break

        # Sample waiting time from Exp(total_rate)
        u = rng.random()
        dt = -math.log(max(u, 1e-12)) / total_rate
        t_next = min(t + dt, total_time)

        # Accumulate measurement statistics for interval [max(t, burn_in), t_next]
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

        # Choose which event fires
        threshold = rng.random() * total_rate
        if can_initiate:
            if threshold < alpha:
                positions.insert(0, 0)
                ribo_ids.insert(0, next_ribo_id)
                next_ribo_id += 1
                if t >= burn_in_seconds and len(positions) > 1 and positions[1] == footprint:
                    collision_events += 1
                continue
            threshold -= alpha

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

    # Normalize metrics over measurement window
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
    phase = _classify_phase(mean_density, collision_fraction, alpha, min(rates))

    # Bottleneck codon: prefer highest collision density, fallback to highest occupancy
    if peak_coll_density > 0.0:
        bottleneck = max(range(L), key=lambda i: (collision_profile[i], codon_occupancy[i]))
    else:
        bottleneck = max(range(L), key=lambda i: codon_occupancy[i])

    return TrafficSimulationResult(
        organism=organism.key,
        codon_count=L,
        initiation_rate=round(alpha, 3),
        termination_rate=round(beta, 3),
        ribosome_footprint=footprint,
        duration_seconds=duration_seconds,
        burn_in_seconds=burn_in_seconds,
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
        codon_occupancy=codon_occupancy,
        collision_profile=collision_profile,
        codon_rates=rates,
        frames=frames,
    )
