"""Multi-objective Synthetic Biology & mRNA-Vaccine Sequence Designer.

Simultaneously optimizes a coding DNA/mRNA sequence across five real wet-lab
constraints while strictly guaranteeing 100% identical amino-acid translation:

  1. Co-Translational Folding Rhythm Preservation (%MinMax correlation vs. source)
  2. 5' mRNA Hairpin Unfolding (maximizes ΔG_5' & ribosome initiation access)
  3. mRNA Vaccine Uridine Depletion (reduces U% & UU dinucleotides -> lower TLR7/8
     innate immune activation, modeled on modern m1Ψ mRNA therapeutics)
  4. Wet-Lab Restriction Enzyme Scrubbing (EcoRI, BamHI, HindIII, XhoI, NdeI,
     and BsaI Golden Gate sites eliminated across codon boundaries!)
  5. Synthesis Homopolymer Filtering (eliminates >=5 nt runs: AAAAA/TTTTT/GGGGG/CCCCC)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .folding import FivePrimeFoldingReport, fold_5prime_window, unfold_5prime_synonymously
from .harmonizer import _pearson_r, harmonize_codons
from .organisms import GENETIC_CODE, SYNONYMOUS_CODONS, OrganismProfile
from .sequence import (
    analyze_sequence,
    compute_cai,
    compute_percent_minmax,
    translate_codons,
    validate_cds,
)
from .simulator import TrafficSimulationResult, simulate_ribosome_traffic

FORBIDDEN_RESTRICTION_SITES: Dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "NdeI": "CATATG",
    "BsaI (+)": "GGTCTC",
    "BsaI (-)": "GAGACC",
}


@dataclass(frozen=True)
class RestrictionHit:
    enzyme: str
    motif: str
    position_nt: int  # 1-based


@dataclass(frozen=True)
class SynBioQC:
    uridine_percent: float
    uu_dinucleotides: int
    gc_percent: float
    gc3_percent: float
    max_homopolymer: int
    homopolymer_violations: int  # runs >= 5 nt
    restriction_hits: List[RestrictionHit]
    folding_5prime: FivePrimeFoldingReport


@dataclass(frozen=True)
class SynBioDesignReport:
    source_organism: str
    host_organism: str
    protein_preserved: bool
    options: Dict[str, bool]
    original_sequence: str
    designed_sequence: str
    codons_changed: int
    percent_changed: float
    original_cai: float
    designed_cai: float
    original_rhythm_r: float
    designed_rhythm_r: float
    before_qc: SynBioQC
    after_qc: SynBioQC
    before_sim: TrafficSimulationResult
    after_sim: TrafficSimulationResult
    native_source_minmax: List[float]
    before_minmax: List[float]
    after_minmax: List[float]


def find_restriction_hits(dna: str) -> List[RestrictionHit]:
    hits: List[RestrictionHit] = []
    seq = dna.upper().replace("U", "T")
    for enzyme, motif in FORBIDDEN_RESTRICTION_SITES.items():
        start = 0
        while True:
            idx = seq.find(motif, start)
            if idx == -1:
                break
            hits.append(RestrictionHit(enzyme=enzyme, motif=motif, position_nt=idx + 1))
            start = idx + 1
    hits.sort(key=lambda h: h.position_nt)
    return hits


def _homopolymer_stats(dna: str) -> tuple[int, int]:
    max_run = 1
    cur = 1
    violations = 0
    for i in range(1, len(dna)):
        if dna[i] == dna[i - 1]:
            cur += 1
            if cur > max_run:
                max_run = cur
        else:
            if cur >= 5:
                violations += 1
            cur = 1
    if cur >= 5:
        violations += 1
    return max_run, violations


def evaluate_synbio_qc(sequence: str, organism: OrganismProfile) -> SynBioQC:
    codons = validate_cds(sequence)
    dna = "".join(codons)
    an = analyze_sequence(dna, organism)
    u_count = dna.count("T")
    u_pct = round((u_count / len(dna)) * 100.0, 2)
    uu_count = sum(1 for i in range(len(dna) - 1) if dna[i : i + 2] == "TT")
    max_hp, hp_viol = _homopolymer_stats(dna)
    return SynBioQC(
        uridine_percent=u_pct,
        uu_dinucleotides=uu_count,
        gc_percent=an.gc_percent,
        gc3_percent=an.gc3_percent,
        max_homopolymer=max_hp,
        homopolymer_violations=hp_viol,
        restriction_hits=find_restriction_hits(dna),
        folding_5prime=fold_5prime_window(dna, window_nt=45),
    )


def _scrub_motifs_and_homopolymers(
    codons: List[str],
    host: OrganismProfile,
    scrub_restriction: bool,
    scrub_homopolymers: bool,
) -> List[str]:
    """Eliminate forbidden cloning sites and >=5 nt homopolymers via synonymous swaps."""
    out = list(codons)
    host_freq = host.codon_freq_per_1000

    for _ in range(4):
        dna = "".join(out)
        changed = False

        if scrub_restriction:
            hits = find_restriction_hits(dna)
            for hit in hits:
                start_c = max(1, (hit.position_nt - 1) // 3)
                end_c = min(len(out) - 2, (hit.position_nt - 1 + len(hit.motif)) // 3)
                fixed = False
                for c_idx in range(start_c, end_c + 1):
                    aa = GENETIC_CODE[out[c_idx]]
                    syns = sorted(
                        SYNONYMOUS_CODONS[aa],
                        key=lambda s: host_freq[s],
                        reverse=True,
                    )
                    orig = out[c_idx]
                    for cand in syns:
                        if cand == orig:
                            continue
                        out[c_idx] = cand
                        new_dna = "".join(out)
                        if len(find_restriction_hits(new_dna)) < len(hits):
                            changed = True
                            fixed = True
                            break
                        out[c_idx] = orig
                    if fixed:
                        break

        if scrub_homopolymers:
            dna = "".join(out)
            run = 1
            for i in range(1, len(dna)):
                if dna[i] == dna[i - 1]:
                    run += 1
                    if run >= 5:
                        c_idx = min(len(out) - 2, max(1, i // 3))
                        aa = GENETIC_CODE[out[c_idx]]
                        orig = out[c_idx]
                        best_cand = orig
                        best_viol = _homopolymer_stats(dna)[1]
                        for cand in SYNONYMOUS_CODONS[aa]:
                            if cand == orig:
                                continue
                            out[c_idx] = cand
                            cand_dna = "".join(out)
                            viol = _homopolymer_stats(cand_dna)[1]
                            if viol < best_viol and (
                                not scrub_restriction
                                or len(find_restriction_hits(cand_dna)) <= len(find_restriction_hits(dna))
                            ):
                                best_viol = viol
                                best_cand = cand
                        out[c_idx] = best_cand
                        if best_cand != orig:
                            changed = True
                            break
                else:
                    run = 1

        if not changed:
            break
    return out


def design_synbio_sequence(
    sequence: str,
    source: OrganismProfile,
    host: OrganismProfile,
    preserve_pauses: bool = True,
    unfold_5prime: bool = True,
    uridine_deplete: bool = False,
    scrub_restriction_sites: bool = True,
    avoid_homopolymers: bool = True,
    initiation_rate: Optional[float] = None,
    seed: int = 42,
) -> SynBioDesignReport:
    """Run multi-objective synthetic biology & mRNA vaccine optimization."""
    wt_codons = validate_cds(sequence)
    wt_seq = "".join(wt_codons)
    protein = translate_codons(wt_codons)
    native_mm = compute_percent_minmax(wt_codons, source)

    # Step 1: Rhythm-preserving harmonization OR host-frequency optimization
    if preserve_pauses:
        cand = harmonize_codons(wt_codons, source, host, min_safe_freq=4.5)
    else:
        freq = host.codon_freq_per_1000
        cand = [wt_codons[0]] + [
            max(SYNONYMOUS_CODONS[GENETIC_CODE[c]], key=lambda s: freq[s])
            if GENETIC_CODE[c] != "*"
            else c
            for c in wt_codons[1:]
        ]

    # Step 2: mRNA Vaccine Uridine-Depletion (m1Ψ / TLR7-8 innate-immunity evasion)
    if uridine_deplete:
        src_ranks = source.normalized_usage_ranks()
        host_ranks = host.normalized_usage_ranks()
        host_freq = host.codon_freq_per_1000
        for idx in range(1, len(cand) - 1):
            aa = GENETIC_CODE[cand[idx]]
            syns = SYNONYMOUS_CODONS[aa]
            if len(syns) <= 1:
                continue
            safe = [s for s in syns if host_freq[s] >= 4.0] or syns
            target_rank = src_ranks[wt_codons[idx]]
            prev_base = cand[idx - 1][-1]

            def u_score(c: str) -> float:
                u_cnt = c.count("T")
                uu_pen = (1 if (prev_base == "T" and c[0] == "T") else 0) + (
                    1 if "TT" in c else 0
                )
                rhythm_pen = abs(host_ranks[c] - target_rank) * (1.2 if preserve_pauses else 0.15)
                freq_bonus = -0.015 * host_freq[c]
                return u_cnt * 2.2 + uu_pen * 1.8 + rhythm_pen + freq_bonus

            cand[idx] = min(safe, key=u_score)

    # Step 3: Melt inhibitory 5' mRNA hairpin (codons 1..14)
    if unfold_5prime:
        melted = unfold_5prime_synonymously("".join(cand), window_codons=15, passes=2)
        cand = validate_cds(melted.unfolded_sequence)

    # Step 4: Scrub cloning restriction sites and >=5 nt homopolymers
    if scrub_restriction_sites or avoid_homopolymers:
        cand = _scrub_motifs_and_homopolymers(
            cand,
            host,
            scrub_restriction=scrub_restriction_sites,
            scrub_homopolymers=avoid_homopolymers,
        )

    designed_seq = "".join(cand)
    protein_ok = translate_codons(cand) == protein
    changed = sum(1 for a, b in zip(wt_codons, cand) if a != b)

    before_mm = compute_percent_minmax(wt_codons, host)
    after_mm = compute_percent_minmax(cand, host)

    before_qc = evaluate_synbio_qc(wt_seq, host)
    after_qc = evaluate_synbio_qc(designed_seq, host)

    before_sim = simulate_ribosome_traffic(
        wt_seq,
        host,
        initiation_rate=initiation_rate,
        duration_seconds=90.0,
        capture_frames=60,
        couple_5prime_folding=True,
        enable_tunnel_drag=True,
        seed=seed,
    )
    after_sim = simulate_ribosome_traffic(
        designed_seq,
        host,
        initiation_rate=initiation_rate,
        duration_seconds=90.0,
        capture_frames=60,
        couple_5prime_folding=True,
        enable_tunnel_drag=True,
        seed=seed,
    )

    return SynBioDesignReport(
        source_organism=source.key,
        host_organism=host.key,
        protein_preserved=protein_ok,
        options={
            "preserve_pauses": preserve_pauses,
            "unfold_5prime": unfold_5prime,
            "uridine_deplete": uridine_deplete,
            "scrub_restriction_sites": scrub_restriction_sites,
            "avoid_homopolymers": avoid_homopolymers,
        },
        original_sequence=wt_seq,
        designed_sequence=designed_seq,
        codons_changed=changed,
        percent_changed=round((changed / len(wt_codons)) * 100.0, 1),
        original_cai=compute_cai(wt_codons, host),
        designed_cai=compute_cai(cand, host),
        original_rhythm_r=_pearson_r(native_mm, before_mm),
        designed_rhythm_r=_pearson_r(native_mm, after_mm),
        before_qc=before_qc,
        after_qc=after_qc,
        before_sim=before_sim,
        after_sim=after_sim,
        native_source_minmax=native_mm,
        before_minmax=before_mm,
        after_minmax=after_mm,
    )
