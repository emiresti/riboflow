"""5' mRNA secondary structure (Turner-calibrated MFE hairpin) & initiation accessibility engine.

Models the biophysical finding (Kudla et al., Science 2009; Goodman et al.,
Science 2013) that tight mRNA hairpins within the first ~45 nucleotides around
the AUG start codon physically occlude ribosome 30S/43S recruitment and dominate
translation initiation efficiency.

Uses a zero-dependency dynamic programming Minimum Free Energy (MFE) RNA folder
calibrated to Turner nearest-neighbor helix stacking energies, hairpin loop
initiation entropy, and internal/multi-branch loop penalties, then couples ΔG_5'
via a Boltzmann accessibility factor directly into the l-TASEP initiation rate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .organisms import GENETIC_CODE, SYNONYMOUS_CODONS
from .sequence import translate_codons, validate_cds

# Canonical RNA base pairs (DNA alphabet T represents U)
_VALID_PAIRS = {
    ("G", "C"), ("C", "G"),
    ("A", "T"), ("T", "A"),
    ("G", "T"), ("T", "G"),
}

# Turner-inspired nearest-neighbor helix stacking free energies (kcal/mol at 37°C)
# Key: (5'-outer, 3'-outer, 5'-inner, 3'-inner)
_STACK_ENERGY: Dict[Tuple[str, str, str, str], float] = {
    ("G", "C", "C", "G"): -3.4,
    ("C", "G", "G", "C"): -2.4,
    ("G", "C", "G", "C"): -3.3,
    ("C", "G", "C", "G"): -3.3,
    ("G", "C", "A", "T"): -2.4,
    ("C", "G", "T", "A"): -2.1,
    ("G", "C", "T", "A"): -2.2,
    ("C", "G", "A", "T"): -2.1,
    ("A", "T", "G", "C"): -2.1,
    ("T", "A", "C", "G"): -2.4,
    ("A", "T", "C", "G"): -2.1,
    ("T", "A", "G", "C"): -2.2,
    ("A", "T", "T", "A"): -1.1,
    ("T", "A", "A", "T"): -1.3,
    ("A", "T", "A", "T"): -0.9,
    ("T", "A", "T", "A"): -0.9,
    ("G", "C", "G", "T"): -2.1,
    ("C", "G", "T", "G"): -1.5,
    ("G", "T", "C", "G"): -2.5,
    ("T", "G", "G", "C"): -1.5,
    ("G", "T", "A", "T"): -1.3,
    ("T", "G", "T", "A"): -1.0,
    ("A", "T", "G", "T"): -1.0,
    ("T", "A", "T", "G"): -1.3,
    ("G", "T", "T", "G"): -0.5,
    ("T", "G", "G", "T"): -0.5,
}

MIN_HAIRPIN_LOOP = 3
UNSTACKED_LOOP_PENALTY = 2.6  # internal/bulge/multiloop initiation penalty (kcal/mol)


@dataclass(frozen=True)
class BasePairArc:
    i: int
    j: int
    pair: str  # e.g. "G-C", "A-U", "G-U"


@dataclass(frozen=True)
class FivePrimeFoldingReport:
    window_nt: int
    rna_sequence: str
    dot_bracket: str
    mfe_kcal_mol: float
    paired_fraction: float
    start_region_paired: float
    initiation_accessibility: float
    severity: str
    arcs: List[BasePairArc]


@dataclass(frozen=True)
class HairpinMeltResult:
    original_sequence: str
    unfolded_sequence: str
    protein_preserved: bool
    codons_changed: int
    before: FivePrimeFoldingReport
    after: FivePrimeFoldingReport
    delta_mfe: float
    accessibility_gain_percent: float


def _hairpin_initiation(loop_len: int) -> float:
    """Entropic initiation cost (kcal/mol) for closing a hairpin loop of size loop_len."""
    if loop_len < MIN_HAIRPIN_LOOP:
        return 25.0
    if loop_len == 3:
        return 5.4
    if loop_len == 4:
        return 4.6
    if loop_len <= 6:
        return 4.9
    return 4.9 + 1.05 * math.log(loop_len / 5.0)


def _boltzmann_accessibility(mfe_kcal_mol: float, start_region_paired: float) -> float:
    """Map 5' MFE (kcal/mol) and AUG-proximal pairing to initiation accessibility [0.15, 1.0].

    Calibrated to Kudla et al. (Science, 2009):
      - ΔG >= -6 kcal/mol  -> ~0.90 - 1.00 (open, optimal ribosome access)
      - ΔG ~ -11 kcal/mol  -> ~0.55 - 0.70 (moderate 5' hairpin)
      - ΔG <= -16 kcal/mol -> ~0.18 - 0.35 (severe 5' translation blockade)
    """
    effective_dg = mfe_kcal_mol - (start_region_paired * 2.5)
    x = -(effective_dg + 12.2) / 2.6
    raw = 1.0 / (1.0 + math.exp(max(-20.0, min(20.0, x))))
    return round(min(1.0, max(0.15, raw)), 3)


def fold_5prime_window(sequence: str, window_nt: int = 45) -> FivePrimeFoldingReport:
    """Compute Turner-calibrated Minimum Free Energy (MFE) structure of the 5' mRNA window."""
    dna = sequence.strip().upper().replace("U", "T")
    sub = dna[: min(len(dna), max(18, window_nt))]
    n = len(sub)
    INF = 1e6

    # V[i][j]: MFE of sub[i..j] given that (i, j) forms a base pair
    V = [[INF] * n for _ in range(n)]
    V_tr: List[List[Tuple[str, int, int] | None]] = [[None] * n for _ in range(n)]

    # W[i][j]: optimal MFE of sub[i..j] (unrestricted; 0.0 = completely unfolded)
    W = [[0.0] * n for _ in range(n)]
    W_tr: List[List[Tuple[str, int, int] | None]] = [[None] * n for _ in range(n)]

    for span in range(MIN_HAIRPIN_LOOP + 1, n):
        for i in range(n - span):
            j = i + span

            # 1) Compute V[i][j] when (i, j) can base-pair
            if (sub[i], sub[j]) in _VALID_PAIRS:
                # Case A: Terminal hairpin loop closed by (i, j)
                best_v = _hairpin_initiation(j - i - 1) - (0.8 if (sub[i], sub[j]) in {("G", "C"), ("C", "G")} else 0.0)
                best_vt: Tuple[str, int, int] = ("H", i, j)

                # Case B: Helix stacking directly onto (i+1, j-1)
                if V[i + 1][j - 1] < INF:
                    stack_e = _STACK_ENERGY.get(
                        (sub[i], sub[j], sub[i + 1], sub[j - 1]), -1.2
                    )
                    cand = V[i + 1][j - 1] + stack_e
                    if cand < best_v:
                        best_v = cand
                        best_vt = ("S", i + 1, j - 1)

                # Case C: Internal loop / bulge / multiloop closed by (i, j)
                if W[i + 1][j - 1] < 0.0:
                    cand = W[i + 1][j - 1] + UNSTACKED_LOOP_PENALTY - (
                        0.6 if (sub[i], sub[j]) in {("G", "C"), ("C", "G")} else 0.0
                    )
                    if cand < best_v:
                        best_v = cand
                        best_vt = ("M", i + 1, j - 1)

                V[i][j] = best_v
                V_tr[i][j] = best_vt

            # 2) Compute W[i][j]
            best_w = 0.0
            best_wt: Tuple[str, int, int] | None = None

            if W[i + 1][j] < best_w:
                best_w = W[i + 1][j]
                best_wt = ("L", i + 1, j)
            if W[i][j - 1] < best_w:
                best_w = W[i][j - 1]
                best_wt = ("R", i, j - 1)
            if V[i][j] < best_w:
                best_w = V[i][j]
                best_wt = ("V", i, j)

            for k in range(i + 1, j):
                cand = W[i][k] + W[k + 1][j]
                if cand < best_w - 1e-9:
                    best_w = cand
                    best_wt = ("B", k, j)

            W[i][j] = best_w
            W_tr[i][j] = best_wt

    # Traceback on (table, i, j)
    pairs: List[Tuple[int, int]] = []
    stack_tb: List[Tuple[str, int, int]] = [("W", 0, n - 1)]
    while stack_tb:
        tbl, i, j = stack_tb.pop()
        if i >= j:
            continue
        if tbl == "W":
            tr = W_tr[i][j]
            if tr is None:
                continue
            kind, a, b = tr
            if kind in ("L", "R"):
                stack_tb.append(("W", a, b))
            elif kind == "V":
                stack_tb.append(("V", a, b))
            elif kind == "B":
                stack_tb.append(("W", i, a))
                stack_tb.append(("W", a + 1, b))
        else:
            tr = V_tr[i][j]
            if tr is None:
                continue
            pairs.append((i, j))
            kind, a, b = tr
            if kind == "S":
                stack_tb.append(("V", a, b))
            elif kind == "M":
                stack_tb.append(("W", a, b))

    pairs.sort()
    chars = ["."] * n
    arcs: List[BasePairArc] = []
    paired_positions = set()

    for i, j in pairs:
        chars[i] = "("
        chars[j] = ")"
        paired_positions.add(i)
        paired_positions.add(j)
        b1 = "U" if sub[i] == "T" else sub[i]
        b2 = "U" if sub[j] == "T" else sub[j]
        arcs.append(BasePairArc(i=i, j=j, pair=f"{b1}-{b2}"))

    mfe = round(min(0.0, W[0][n - 1]), 2)
    paired_frac = round(len(paired_positions) / n, 3) if n else 0.0
    start_span = min(15, n)
    start_paired = round(
        sum(1 for pos in range(start_span) if pos in paired_positions) / start_span, 3
    )
    acc = _boltzmann_accessibility(mfe, start_paired)

    if acc >= 0.80:
        severity = "Open (Optimal Initiation)"
    elif acc >= 0.50:
        severity = "Moderate 5' Hairpin"
    else:
        severity = "Severe 5' Blockade"

    return FivePrimeFoldingReport(
        window_nt=n,
        rna_sequence=sub.replace("T", "U"),
        dot_bracket="".join(chars),
        mfe_kcal_mol=mfe,
        paired_fraction=paired_frac,
        start_region_paired=start_paired,
        initiation_accessibility=acc,
        severity=severity,
        arcs=arcs,
    )


def unfold_5prime_synonymously(
    sequence: str,
    window_codons: int = 15,
    passes: int = 2,
) -> HairpinMeltResult:
    """Melt inhibitory 5' mRNA hairpins using strictly synonymous swaps in codons 1..14."""
    codons = validate_cds(sequence)
    wt_protein = translate_codons(codons)
    before = fold_5prime_window("".join(codons), window_nt=window_codons * 3)

    best_codons = list(codons)
    limit = min(len(best_codons) - 1, window_codons)
    cur_report = before

    for _ in range(passes):
        improved = False
        for idx in range(1, limit):
            orig = best_codons[idx]
            aa = GENETIC_CODE[orig]
            syns = SYNONYMOUS_CODONS[aa]
            if len(syns) <= 1:
                continue

            best_choice = orig
            best_score = (
                cur_report.initiation_accessibility,
                cur_report.mfe_kcal_mol,
                -cur_report.start_region_paired,
            )

            for cand in syns:
                if cand == orig:
                    continue
                best_codons[idx] = cand
                rep = fold_5prime_window("".join(best_codons), window_nt=window_codons * 3)
                sc = (rep.initiation_accessibility, rep.mfe_kcal_mol, -rep.start_region_paired)
                if sc > best_score:
                    best_score = sc
                    best_choice = cand
                    cur_report = rep
                    improved = True
            best_codons[idx] = best_choice
        if not improved:
            break

    unfolded_seq = "".join(best_codons)
    after = fold_5prime_window(unfolded_seq, window_nt=window_codons * 3)
    changed = sum(1 for a, b in zip(codons, best_codons) if a != b)
    gain_pct = (
        round(
            ((after.initiation_accessibility - before.initiation_accessibility) / before.initiation_accessibility)
            * 100.0,
            1,
        )
        if before.initiation_accessibility > 0
        else 0.0
    )

    return HairpinMeltResult(
        original_sequence="".join(codons),
        unfolded_sequence=unfolded_seq,
        protein_preserved=(translate_codons(best_codons) == wt_protein),
        codons_changed=changed,
        before=before,
        after=after,
        delta_mfe=round(after.mfe_kcal_mol - before.mfe_kcal_mol, 2),
        accessibility_gain_percent=gain_pct,
    )
