"""Organism codon usage tables and translation elongation kinetics.

Codon usage frequencies (per 1,000 codons) are derived from the Kazusa Codon
Usage Database for four model organisms widely used in molecular biology and
recombinant protein expression:
  - Homo sapiens (human)
  - Escherichia coli K-12 (ecoli)
  - Saccharomyces cerevisiae (yeast)
  - Mus musculus (mouse)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

# Standard genetic code (DNA alphabet)
GENETIC_CODE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AMINO_ACID_NAMES: Dict[str, str] = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "Stop",
}

SYNONYMOUS_CODONS: Dict[str, List[str]] = {}
for _codon, _aa in GENETIC_CODE.items():
    SYNONYMOUS_CODONS.setdefault(_aa, []).append(_codon)


@dataclass(frozen=True)
class OrganismProfile:
    """Codon usage and physiological translation parameters for an organism."""

    key: str
    name: str
    common_name: str
    ribosome_type: str
    mean_elongation_rate: float  # codons / second (physiological mean)
    ribosome_footprint: int      # codons occluded per ribosome (~30 nt = 10 codons)
    default_initiation_rate: float  # initiation events / second
    codon_freq_per_1000: Dict[str, float]

    def relative_adaptiveness(self) -> Dict[str, float]:
        """Sharp & Li relative adaptiveness w_ij = f_ij / max(f_i) for CAI."""
        w: Dict[str, float] = {}
        for aa, codons in SYNONYMOUS_CODONS.items():
            if aa == "*":
                for c in codons:
                    w[c] = 1.0
                continue
            max_f = max(self.codon_freq_per_1000[c] for c in codons)
            for c in codons:
                w[c] = self.codon_freq_per_1000[c] / max_f if max_f > 0 else 1.0
        return w

    def rscu_table(self) -> Dict[str, float]:
        """Relative Synonymous Codon Usage (RSCU) for every sense codon."""
        rscu: Dict[str, float] = {}
        for aa, codons in SYNONYMOUS_CODONS.items():
            if aa == "*":
                for c in codons:
                    rscu[c] = 1.0
                continue
            total = sum(self.codon_freq_per_1000[c] for c in codons)
            n = len(codons)
            for c in codons:
                rscu[c] = (self.codon_freq_per_1000[c] * n / total) if total > 0 else 1.0
        return rscu

    def normalized_usage_ranks(self) -> Dict[str, float]:
        """Normalized intra-synonymous frequency in [0, 1] for codon harmonization.

        For a single-codon amino acid (M, W), returns 1.0.
        For degenerate amino acids, returns (f - f_min) / (f_max - f_min) so
        the rarest synonymous codon maps to 0.0 and the most frequent maps to 1.0.
        """
        ranks: Dict[str, float] = {}
        for aa, codons in SYNONYMOUS_CODONS.items():
            if aa == "*" or len(codons) == 1:
                for c in codons:
                    ranks[c] = 1.0
                continue
            freqs = [self.codon_freq_per_1000[c] for c in codons]
            f_min, f_max = min(freqs), max(freqs)
            span = f_max - f_min
            for c in codons:
                ranks[c] = (self.codon_freq_per_1000[c] - f_min) / span if span > 0 else 1.0
        return ranks

    def elongation_rates(self) -> Dict[str, float]:
        """Estimated per-codon decoding rate (codons/s) grounded in codon frequency.

        Combines global codon abundance (proxy for cognate tRNA availability)
        and intra-synonymous relative adaptiveness, scaled around the organism's
        physiological mean elongation rate (e.g. ~5.6 aa/s in human, ~15 aa/s in E. coli).
        Rare codons (e.g., AGG/AGA/CTA in E. coli) naturally become kinetic bottlenecks.
        """
        w = self.relative_adaptiveness()
        sense_freqs = [
            self.codon_freq_per_1000[c]
            for c, aa in GENETIC_CODE.items()
            if aa != "*"
        ]
        median_freq = sorted(sense_freqs)[len(sense_freqs) // 2]

        raw_scores: Dict[str, float] = {}
        for codon, aa in GENETIC_CODE.items():
            if aa == "*":
                continue
            freq = max(self.codon_freq_per_1000[codon], 0.4)
            # Blend global tRNA supply proxy (freq/median) and synonymous adaptiveness w
            supply_ratio = min(freq / median_freq, 2.2)
            score = 0.65 * supply_ratio + 0.35 * w[codon]
            # Smooth lower bound to avoid zero rate while preserving strong stalls
            raw_scores[codon] = max(0.08, score)

        mean_raw = sum(raw_scores.values()) / len(raw_scores)
        scale = self.mean_elongation_rate / mean_raw

        rates = {c: round(v * scale, 3) for c, v in raw_scores.items()}
        for stop_codon in SYNONYMOUS_CODONS["*"]:
            rates[stop_codon] = round(self.mean_elongation_rate * 0.8, 3)
        return rates


# Kazusa-derived codon usage frequencies per 1,000 codons
_HUMAN_FREQ: Dict[str, float] = {
    "TTT": 17.6, "TTC": 20.3, "TTA": 7.7,  "TTG": 12.9,
    "CTT": 13.2, "CTC": 19.6, "CTA": 7.2,  "CTG": 39.6,
    "ATT": 16.0, "ATC": 20.8, "ATA": 7.5,  "ATG": 22.0,
    "GTT": 11.0, "GTC": 14.5, "GTA": 7.1,  "GTG": 28.1,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
    "TAT": 12.2, "TAC": 15.3, "TAA": 1.0,  "TAG": 0.8,
    "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
    "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
    "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
    "TGT": 10.6, "TGC": 12.6, "TGA": 1.6,  "TGG": 13.2,
    "CGT": 4.5,  "CGC": 10.4, "CGA": 6.2,  "CGG": 11.4,
    "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
}

_ECOLI_FREQ: Dict[str, float] = {
    "TTT": 22.4, "TTC": 16.6, "TTA": 13.9, "TTG": 13.7,
    "CTT": 11.0, "CTC": 11.1, "CTA": 3.9,  "CTG": 52.6,
    "ATT": 30.3, "ATC": 25.1, "ATA": 4.4,  "ATG": 27.8,
    "GTT": 18.3, "GTC": 15.3, "GTA": 10.9, "GTG": 26.2,
    "TCT": 8.5,  "TCC": 8.6,  "TCA": 7.2,  "TCG": 8.9,
    "CCT": 7.0,  "CCC": 5.5,  "CCA": 8.4,  "CCG": 23.2,
    "ACT": 9.0,  "ACC": 23.4, "ACA": 7.1,  "ACG": 14.4,
    "GCT": 15.3, "GCC": 25.5, "GCA": 20.2, "GCG": 33.6,
    "TAT": 16.2, "TAC": 12.2, "TAA": 2.0,  "TAG": 0.3,
    "CAT": 12.9, "CAC": 9.7,  "CAA": 15.3, "CAG": 28.8,
    "AAT": 17.7, "AAC": 21.7, "AAA": 33.6, "AAG": 10.3,
    "GAT": 32.1, "GAC": 19.1, "GAA": 39.4, "GAG": 17.8,
    "TGT": 5.2,  "TGC": 6.5,  "TGA": 1.0,  "TGG": 15.2,
    "CGT": 20.9, "CGC": 22.0, "CGA": 3.6,  "CGG": 5.4,
    "AGT": 8.8,  "AGC": 16.1, "AGA": 2.1,  "AGG": 1.2,
    "GGT": 24.7, "GGC": 29.6, "GGA": 8.0,  "GGG": 11.1,
}

_YEAST_FREQ: Dict[str, float] = {
    "TTT": 26.1, "TTC": 18.4, "TTA": 26.2, "TTG": 27.2,
    "CTT": 12.3, "CTC": 5.4,  "CTA": 13.4, "CTG": 10.5,
    "ATT": 30.1, "ATC": 17.2, "ATA": 17.8, "ATG": 20.9,
    "GTT": 22.1, "GTC": 11.8, "GTA": 11.8, "GTG": 10.8,
    "TCT": 23.5, "TCC": 14.2, "TCA": 18.7, "TCG": 8.6,
    "CCT": 13.5, "CCC": 6.8,  "CCA": 18.3, "CCG": 5.3,
    "ACT": 20.3, "ACC": 12.7, "ACA": 17.8, "ACG": 8.0,
    "GCT": 21.2, "GCC": 12.6, "GCA": 16.2, "GCG": 6.2,
    "TAT": 18.8, "TAC": 14.8, "TAA": 1.1,  "TAG": 0.5,
    "CAT": 13.6, "CAC": 7.8,  "CAA": 27.3, "CAG": 12.1,
    "AAT": 35.7, "AAC": 24.8, "AAA": 41.9, "AAG": 30.8,
    "GAT": 37.6, "GAC": 20.2, "GAA": 45.6, "GAG": 19.2,
    "TGT": 8.1,  "TGC": 4.8,  "TGA": 0.7,  "TGG": 10.4,
    "CGT": 6.4,  "CGC": 2.6,  "CGA": 3.0,  "CGG": 1.7,
    "AGT": 14.2, "AGC": 9.8,  "AGA": 21.3, "AGG": 9.2,
    "GGT": 23.9, "GGC": 9.8,  "GGA": 10.9, "GGG": 6.0,
}

_MOUSE_FREQ: Dict[str, float] = {
    "TTT": 17.2, "TTC": 21.8, "TTA": 6.7,  "TTG": 13.4,
    "CTT": 13.4, "CTC": 20.2, "CTA": 8.1,  "CTG": 39.5,
    "ATT": 15.4, "ATC": 22.5, "ATA": 7.4,  "ATG": 22.8,
    "GTT": 10.7, "GTC": 15.4, "GTA": 7.4,  "GTG": 28.4,
    "TCT": 16.2, "TCC": 18.1, "TCA": 11.8, "TCG": 4.2,
    "CCT": 18.4, "CCC": 18.2, "CCA": 17.3, "CCG": 6.2,
    "ACT": 13.7, "ACC": 19.0, "ACA": 16.0, "ACG": 5.6,
    "GCT": 20.0, "GCC": 26.0, "GCA": 15.8, "GCG": 6.4,
    "TAT": 12.2, "TAC": 16.1, "TAA": 1.0,  "TAG": 0.8,
    "CAT": 10.6, "CAC": 15.3, "CAA": 12.0, "CAG": 34.1,
    "AAT": 15.6, "AAC": 20.3, "AAA": 21.9, "AAG": 33.6,
    "GAT": 21.0, "GAC": 26.0, "GAA": 27.0, "GAG": 39.4,
    "TGT": 11.4, "TGC": 12.3, "TGA": 1.6,  "TGG": 12.5,
    "CGT": 4.7,  "CGC": 9.4,  "CGA": 6.6,  "CGG": 10.2,
    "AGT": 12.7, "AGC": 19.7, "AGA": 12.1, "AGG": 12.2,
    "GGT": 11.4, "GGC": 21.2, "GGA": 16.8, "GGG": 15.2,
}

ORGANISMS: Dict[str, OrganismProfile] = {
    "human": OrganismProfile(
        key="human",
        name="Homo sapiens",
        common_name="Human",
        ribosome_type="80S Eukaryotic",
        mean_elongation_rate=5.6,
        ribosome_footprint=10,
        default_initiation_rate=0.32,
        codon_freq_per_1000=_HUMAN_FREQ,
    ),
    "ecoli": OrganismProfile(
        key="ecoli",
        name="Escherichia coli K-12",
        common_name="E. coli",
        ribosome_type="70S Prokaryotic",
        mean_elongation_rate=15.0,
        ribosome_footprint=10,
        default_initiation_rate=0.85,
        codon_freq_per_1000=_ECOLI_FREQ,
    ),
    "yeast": OrganismProfile(
        key="yeast",
        name="Saccharomyces cerevisiae",
        common_name="Baker's Yeast",
        ribosome_type="80S Eukaryotic",
        mean_elongation_rate=9.5,
        ribosome_footprint=10,
        default_initiation_rate=0.55,
        codon_freq_per_1000=_YEAST_FREQ,
    ),
    "mouse": OrganismProfile(
        key="mouse",
        name="Mus musculus",
        common_name="Mouse",
        ribosome_type="80S Eukaryotic",
        mean_elongation_rate=5.8,
        ribosome_footprint=10,
        default_initiation_rate=0.33,
        codon_freq_per_1000=_MOUSE_FREQ,
    ),
}


def get_organism(key: str) -> OrganismProfile:
    """Return an OrganismProfile by key (case-insensitive)."""
    normalized = key.strip().lower()
    if normalized not in ORGANISMS:
        valid = ", ".join(sorted(ORGANISMS.keys()))
        raise ValueError(f"Unknown organism '{key}'. Valid options: {valid}")
    return ORGANISMS[normalized]
