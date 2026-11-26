"""Coding sequence (CDS) validation, translation, and codon-usage metrics."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .organisms import (
    AMINO_ACID_NAMES,
    GENETIC_CODE,
    SYNONYMOUS_CODONS,
    OrganismProfile,
)


@dataclass(frozen=True)
class CodonMetrics:
    index: int
    codon: str
    amino_acid: str
    aa_three: str
    frequency_per_1000: float
    relative_adaptiveness: float
    rscu: float
    elongation_rate: float
    minmax_score: float  # %MinMax value in [-100, +100]
    synonymous_options: List[str]


@dataclass(frozen=True)
class SequenceAnalysis:
    organism: str
    sequence: str
    protein: str
    codon_count: int
    gc_percent: float
    gc3_percent: float
    cai: float
    mean_elongation_rate: float
    min_minmax: float
    max_minmax: float
    pause_sites: List[int]  # 0-based codon indices in strong %MinMax troughs
    codons: List[CodonMetrics]


def parse_fasta_or_raw(text: str) -> Tuple[str, str]:
    """Extract (header, cleaned DNA sequence) from FASTA or raw DNA/mRNA input."""
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("Input sequence is empty.")

    header = "Custom Sequence"
    seq_parts: List[str] = []
    for line in lines:
        if line.startswith(">"):
            if not seq_parts:
                header = line[1:].strip() or header
            else:
                # Only read first FASTA record
                break
        else:
            seq_parts.append(line)

    raw = "".join(seq_parts).upper().replace(" ", "").replace("\t", "")
    # Normalize RNA (U) to DNA (T)
    dna = raw.replace("U", "T")
    return header, dna


def validate_cds(sequence: str) -> List[str]:
    """Validate a coding DNA sequence and return its list of codons."""
    dna = sequence.strip().upper().replace("U", "T").replace(" ", "").replace("\n", "")
    if not dna:
        raise ValueError("Sequence cannot be empty.")

    invalid_chars = sorted(set(dna) - {"A", "T", "G", "C"})
    if invalid_chars:
        raise ValueError(
            f"Sequence contains non-canonical nucleotide(s): {', '.join(invalid_chars)}. "
            "Only A, T/U, G, C are supported."
        )

    if len(dna) % 3 != 0:
        raise ValueError(
            f"Coding sequence length ({len(dna)} nt) must be a multiple of 3."
        )

    codons = [dna[i : i + 3] for i in range(0, len(dna), 3)]
    if len(codons) < 15:
        raise ValueError(
            f"Sequence is too short ({len(codons)} codons). Provide at least 15 codons "
            "so ribosome footprint (10 codons) kinetics can be modeled."
        )

    if codons[0] != "ATG":
        raise ValueError(
            f"Coding sequence must begin with start codon ATG (found '{codons[0]}')."
        )

    internal_stops = [
        idx for idx, c in enumerate(codons[:-1]) if GENETIC_CODE[c] == "*"
    ]
    if internal_stops:
        first = internal_stops[0]
        raise ValueError(
            f"Premature internal stop codon '{codons[first]}' found at codon position {first + 1}."
        )

    return codons


def translate_codons(codons: List[str]) -> str:
    """Translate codons to single-letter amino acid sequence (excluding terminal stop)."""
    aas = [GENETIC_CODE[c] for c in codons]
    if aas and aas[-1] == "*":
        aas = aas[:-1]
    return "".join(aas)


def compute_cai(codons: List[str], organism: OrganismProfile) -> float:
    """Compute Sharp & Li Codon Adaptation Index (CAI) over sense codons."""
    w = organism.relative_adaptiveness()
    sense = [c for c in codons if GENETIC_CODE[c] != "*"]
    if not sense:
        return 0.0
    log_sum = sum(math.log(max(w[c], 1e-6)) for c in sense)
    return round(math.exp(log_sum / len(sense)), 4)


def compute_percent_minmax(
    codons: List[str],
    organism: OrganismProfile,
    window: int = 11,
) -> List[float]:
    """Compute Clarke & Clark (2008) %MinMax sliding-window codon usage rhythm.

    Positive values (+1 to +100) indicate clusters of frequent (fast) codons;
    negative values (-1 to -100) indicate clusters of rare (slow/pause) codons.
    """
    n = len(codons)
    half = window // 2
    freq_table = organism.codon_freq_per_1000

    # Precompute min, max, avg synonymous frequency per amino acid
    aa_stats: Dict[str, Tuple[float, float, float]] = {}
    for aa, syn in SYNONYMOUS_CODONS.items():
        freqs = [freq_table[c] for c in syn]
        aa_stats[aa] = (min(freqs), max(freqs), sum(freqs) / len(freqs))

    scores: List[float] = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        win = [c for c in codons[lo:hi] if GENETIC_CODE[c] != "*"]
        if not win:
            scores.append(0.0)
            continue

        actual = sum(freq_table[c] for c in win) / len(win)
        f_min = sum(aa_stats[GENETIC_CODE[c]][0] for c in win) / len(win)
        f_max = sum(aa_stats[GENETIC_CODE[c]][1] for c in win) / len(win)
        f_avg = sum(aa_stats[GENETIC_CODE[c]][2] for c in win) / len(win)

        if actual >= f_avg:
            denom = f_max - f_avg
            val = ((actual - f_avg) / denom) * 100.0 if denom > 1e-9 else 0.0
        else:
            denom = f_avg - f_min
            val = -((f_avg - actual) / denom) * 100.0 if denom > 1e-9 else 0.0

        scores.append(round(max(-100.0, min(100.0, val)), 2))

    return scores


def analyze_sequence(
    sequence: str,
    organism: OrganismProfile,
    minmax_window: int = 11,
) -> SequenceAnalysis:
    """Perform full codon-level sequence and folding-rhythm analysis."""
    codons = validate_cds(sequence)
    dna = "".join(codons)
    protein = translate_codons(codons)

    w = organism.relative_adaptiveness()
    rscu = organism.rscu_table()
    rates = organism.elongation_rates()
    minmax = compute_percent_minmax(codons, organism, window=minmax_window)

    gc_count = dna.count("G") + dna.count("C")
    gc_percent = round((gc_count / len(dna)) * 100.0, 2)

    sense_codons = [c for c in codons if GENETIC_CODE[c] != "*"]
    gc3_count = sum(1 for c in sense_codons if c[2] in ("G", "C"))
    gc3_percent = round((gc3_count / len(sense_codons)) * 100.0, 2) if sense_codons else 0.0

    cai = compute_cai(codons, organism)
    mean_rate = round(
        sum(rates[c] for c in sense_codons) / len(sense_codons), 3
    )

    # Identify local pause troughs (%MinMax <= -15)
    pause_sites = [i for i, score in enumerate(minmax) if score <= -15.0]

    codon_metrics: List[CodonMetrics] = []
    for idx, c in enumerate(codons):
        aa = GENETIC_CODE[c]
        codon_metrics.append(
            CodonMetrics(
                index=idx,
                codon=c,
                amino_acid=aa,
                aa_three=AMINO_ACID_NAMES[aa],
                frequency_per_1000=organism.codon_freq_per_1000[c],
                relative_adaptiveness=round(w[c], 3),
                rscu=round(rscu[c], 3),
                elongation_rate=rates[c],
                minmax_score=minmax[idx],
                synonymous_options=list(SYNONYMOUS_CODONS[aa]),
            )
        )

    return SequenceAnalysis(
        organism=organism.key,
        sequence=dna,
        protein=protein,
        codon_count=len(codons),
        gc_percent=gc_percent,
        gc3_percent=gc3_percent,
        cai=cai,
        mean_elongation_rate=mean_rate,
        min_minmax=min(minmax),
        max_minmax=max(minmax),
        pause_sites=pause_sites,
        codons=codon_metrics,
    )
