"""Curated biological coding sequences (CDS) for immediate exploration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class SampleGene:
    id: str
    name: str
    organism: str
    recommended_host: str
    summary_en: str
    summary_tr: str
    sequence: str


SAMPLES: Dict[str, SampleGene] = {
    "human_ins": SampleGene(
        id="human_ins",
        name="Human Preproinsulin (INS)",
        organism="human",
        recommended_host="ecoli",
        summary_en=(
            "Full 110-aa human preproinsulin CDS. Classic recombinant expression target "
            "for E. coli; contains Arginine (AGA/AGG) and Proline clusters that stall "
            "70S ribosomes unless harmonized."
        ),
        summary_tr=(
            "110 amino asitlik tam insan preproinsülin kodlayan dizisi. E. coli'de "
            "rekombinant üretimde klasik modeldir; harmonizasyon yapılmazsa nadir "
            "Arjinin (AGA/AGG) kodonlarında 70S ribozom sıkışması görülür."
        ),
        sequence=(
            "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGAC"
            "CCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTAC"
            "CTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGAC"
            "CTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTG"
            "GCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGC"
            "TCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
        ),
    ),
    "sars_cov2_rbd": SampleGene(
        id="sars_cov2_rbd",
        name="SARS-CoV-2 Spike RBD (mRNA Vaccine Model)",
        organism="human",
        recommended_host="human",
        summary_en=(
            "Native AT/U-rich Spike receptor-binding domain segment (34.8% Uridine). "
            "Ideal for testing mRNA Vaccine Uridine-Depletion (m1Ψ / TLR7-8 immune evasion) "
            "and 5' hairpin unfolding in the SynBio Designer."
        ),
        summary_tr=(
            "Doğal Uridin-zengin SARS-CoV-2 Spike RBD segmenti (%34.8 Uridin). "
            "mRNA Aşı Tasarımcısında Uridin azaltma (TLR7/8 bağışıklık kaçışı) ve "
            "5' saç tokası açma optimizasyonunu denemek için idealdir."
        ),
        sequence=(
            "ATGAGAGTCCAACCAACAGAATCTATTGTTAGATTTCCTAATATTACAAACTTGTGCCCT"
            "TTTGGTGAAGTTTTTAACGCCACCAGATTTGCATCTGTTTATGCTTGGAACAGGAAGAGA"
            "ATCAGCAACTGTGTTGCTGATTATTCTGTCCTATATAATTCCGCATCATTTTCCACTTTT"
            "AAGTGTTATGGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCA"
            "GATTCATTTGTAATTAGAGGTGATGAAGTCAGACAAATCGCTCCAGGGCAAACTGGAAAG"
            "ATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGA"
        ),
    ),
    "polybasic_arrest": SampleGene(
        id="polybasic_arrest",
        name="5' Hairpin + Polybasic Exit-Tunnel Arrest + Cloning Trap",
        organism="human",
        recommended_host="ecoli",
        summary_en=(
            "Synthetic stress-test gene combining three real wet-lab failure modes: "
            "a tight 5' GC hairpin blocking AUG initiation, a poly-Arg/Lys (KRKKRKR) "
            "exit-tunnel electrostatic stall tract, and forbidden EcoRI/BamHI/BsaI cloning sites."
        ),
        summary_tr=(
            "Üç gerçek ıslak-lab sorununu birleştiren sentetik stres testi: AUG'yi kilitleyen "
            "sıkı 5' GC saç tokası, çıkış tünelinde elektrostatik frenleme yapan polibazik "
            "(KRKKRKR) küme ve yasaklı EcoRI/BamHI/BsaI klonlama kesim bölgeleri!"
        ),
        sequence=(
            "ATGGCGGCCGCCGGCGGCGGCCGCCATGGCCTGAAAGCTGAATTCCTGAAAGCTGGATCC"
            "GCTGTTCTGAAAGCTGCTGAAGTTCTGGCTCTGAAAGGTCTCCTGAAAGCTGTTCTGAAA"
            "AAGCGTAAGAAACGTAAACGTAAGAAGCGCAAAGCTCTGAAAGCTGTTCTGAAAGCTGCT"
            "GAAGTTCTGGCTCTGAAAGCTGTTCTGAAAGCTGCTGAAGTTCTGGCTCTGAAAGTTTAA"
        ),
    ),
    "cftr_nbd1": SampleGene(
        id="cftr_nbd1",
        name="Human CFTR NBD1 Folding Segment",
        organism="human",
        recommended_host="ecoli",
        summary_en=(
            "Nucleotide-binding domain 1 segment around the F508/I507 folding locus. "
            "Contains conserved slow-codon regions important for co-translational "
            "folding—ideal for testing how a single synonymous mutation alters kinetics."
        ),
        summary_tr=(
            "F508/I507 katlanma bölgesi çevresindeki CFTR NBD1 segmenti. Eş-zamanlı "
            "katlanma (co-translational folding) için kritik yavaş kodon kümeleri "
            "içerir; tek bir sessiz mutasyonun etkisini test etmek için idealdir."
        ),
        sequence=(
            "ATGATCATGGAGAACATCCAGAGCGTTAAGGCTTACTGCTGGGAAGAAGCAATGGAAAAA"
            "ATGATTGAAAACTTAAGACAAACAGAACTGAAACTGACTCGGAAGGCAGCCTATGTGAGA"
            "TACTTCAATAGCTCAGCCTTCTTCTTCTCAGGGTTCTTTGTGGTGTTTTTATCTGTGCTT"
            "CCCTATGCACTAATCAAAGGAATCATCCTCCGGAAAATATTCACCACCATCTCATTCTGC"
            "ATTGTTCTGCGCATGGCGGTCACTCGGCAATTTCCCTGGGCTGTACAAACATGGTATGAC"
            "TCTCTTGGAGCAATAAACAAAATACAGGATTTCTTACAAAAGCAAGAATATAAGACATAG"
        ),
    ),
    "egfp_core": SampleGene(
        id="egfp_core",
        name="eGFP Chromophore Core (A. victoria -> Humanized)",
        organism="human",
        recommended_host="yeast",
        summary_en=(
            "Human codon-optimized GFP N-terminal core encompassing the SYG chromophore "
            "triad. Expressing this GC-rich sequence in yeast reveals strong organismal "
            "codon mismatch."
        ),
        summary_tr=(
            "İnsan için optimize edilmiş eGFP'nin SYG kromofor bölgesini içeren "
            "çekirdek dizisi. Maya (S. cerevisiae) konağına aktarıldığında belirgin "
            "kodon uyumsuzluğu ve duraklama farkları oluşturur."
        ),
        sequence=(
            "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGAC"
            "GGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTAC"
            "GGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACC"
            "CTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG"
            "CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTC"
            "TTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCTAA"
        ),
    ),
    "jam_benchmark": SampleGene(
        id="jam_benchmark",
        name="Ribosome Traffic-Jam Benchmark (Rare Bottleneck)",
        organism="human",
        recommended_host="ecoli",
        summary_en=(
            "Synthetic benchmark gene featuring a high-speed 5' ramp followed by a "
            "rare Arg/Ile/Leu cluster (AGG-AGA-ATA-CTA-CGA) in E. coli—triggers dramatic "
            "ribosome queueing and disome collisions for physics exploration."
        ),
        summary_tr=(
            "Hızlı 5' başlangıç rampasını takip eden nadir Arg/Ile/Leu kümesi "
            "(AGG-AGA-ATA-CTA-CGA) içeren sentetik test dizisi. E. coli'de belirgin "
            "ribozom kuyruklanması ve disome çarpışmalarını gözlemlemek için idealdir."
        ),
        sequence=(
            "ATGGCTCTGCTGAAAGCTGTTGCTGCTGAAAAACTGCTGGTTAAAGCTGCTGAAGTTCTG"
            "GCTCTGAAAGCTGAAGTTCTGCTGGCTGTTGAAAAAGCTGCTCTGAAAGTTGCTGCTGAA"
            "CTGAAAGCTGTTCTGAAAGCTAGGAGAATACTACGAAGGAGAATACTACGAAGGAGAATA"
            "GCTGTTCTGAAAGCTGCTGAAGTTCTGGCTCTGAAAGCTGTTCTGAAAGCTGAAGTTTAA"
        ),
    ),
}
