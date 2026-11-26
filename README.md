# RiboFlow

[![Live Demo](https://img.shields.io/badge/▶_Live_Demo-emiresti.github.io/riboflow-2dd4bf?style=for-the-badge&labelColor=0b1022)](https://emiresti.github.io/riboflow/) [![CI](https://github.com/emiresti/riboflow/actions/workflows/ci.yml/badge.svg)](https://github.com/emiresti/riboflow/actions/workflows/ci.yml)

**Co-Translational Ribosome Traffic ($\ell$-TASEP) Simulator, Synonymous Mutation Impact Lab & Pause-Preserving Codon Harmonizer**

> Why can a single *"silent"* DNA mutation—changing **zero** amino acids—disrupt protein folding or trigger mRNA degradation?
> And why does naive codon optimization frequently produce misfolded inclusion bodies when expressing human genes in *E. coli*?

**RiboFlow** bridges **molecular biology** and **computational statistical physics** in a zero-dependency Python package, CLI, and interactive web laboratory.

---

## Key Features

### 1. Exact Continuous-Time $\ell$-TASEP Ribosome Traffic Simulator
- Simulates individual ribosomes (steric footprint $\ell = 10$ codons $\approx 30\text{ nt}$) translating an mRNA strand using a kinetic Monte Carlo / Gillespie **Totally Asymmetric Simple Exclusion Process ($\ell$-TASEP)**.
- Models codon-specific decoding rates across four model organisms (*Homo sapiens*, *Escherichia coli* K-12, *Saccharomyces cerevisiae*, and *Mus musculus*).
- Quantifies **protein production flux ($J$)**, **simulated Ribo-seq occupancy**, **Disome-seq steric collision profiles**, **traffic phase regime** (*Low-Density*, *Maximum-Current*, *High-Density Jam*), and **No-Go Decay (NGD) / RQC activation risk**.

### 2. Synonymous ("Silent") Mutation Impact Laboratory
- Mutate any degenerate codon to a synonymous alternative while guaranteeing **100% amino-acid sequence identity**.
- Measure exact kinetic shifts ($\Delta J$, $\Delta$ collision time, $\Delta$ NGD risk, and Ribo-seq peak shifts).
- Includes an automated **Critical Silent Hotspot Scanner** that discovers single-nucleotide synonymous substitutions most likely to trigger ribosome pileups.

### 3. Pause-Preserving Codon Harmonizer
- Naive CAI maximization replaces every codon with the host's fastest codon—flattening conserved translational pause sites needed for co-translational protein folding.
- **RiboFlow Harmonizer** maps each source codon's intra-synonymous usage rank to the matching rank in the target host while avoiding catastrophic ultra-rare dead-stop codons.
- Benchmarks **Wild-Type vs. Naive CAI Optimized vs. RiboFlow Harmonized** side-by-side on `%MinMax` folding-rhythm correlation ($r$), protein flux, and disome collision risk.

### 4. Interactive Live Web Laboratory (TR / EN)
- Real-time HTML5 Canvas animation showing 70S/80S ribosomes moving along the mRNA with growing nascent polypeptide chains and glowing crimson disome collision pulses.
- Interactive codon strip to click-and-mutate codons live, synchronized Ribo-seq / `%MinMax` charts, and one-click FASTA export.

---

## Quick Start

RiboFlow requires **Python 3.10+** and uses **zero external runtime dependencies**.

```bash
# 1. Clone and install in editable mode
git clone https://github.com/emiresti/riboflow.git
cd riboflow
pip install -e .

# 2. Launch the interactive web laboratory
riboflow serve --open
```

Or run directly without installing:

```bash
PYTHONPATH=src python -m riboflow serve --open
```

---

## CLI Usage

```bash
# Analyze CAI, GC3, and %MinMax pause rhythm of Human Preproinsulin in E. coli
riboflow analyze --sample human_ins --organism ecoli

# Run l-TASEP ribosome traffic simulation on the rare-codon bottleneck benchmark
riboflow simulate --sample jam_benchmark --organism ecoli --initiation-rate 1.2

# Compare Wild-Type vs. Naive CAI vs. Pause-Preserving Codon Harmonization (Human -> E. coli)
riboflow harmonize --sample human_ins --source human --host ecoli --output-fasta ins_harmonized.fasta

# Scan a gene for the most dangerous single-nucleotide synonymous ("silent") mutations
riboflow mutate --sample cftr_nbd1 --organism human --scan

# Apply a specific synonymous mutation (1-based codon index)
riboflow mutate --sample human_ins --organism ecoli --edit 44:AGG
```

---

## Python API

```python
from riboflow import (
    SAMPLES,
    get_organism,
    simulate_ribosome_traffic,
    compare_harmonization,
    scan_critical_silent_mutations,
)

seq = SAMPLES["human_ins"].sequence
human = get_organism("human")
ecoli = get_organism("ecoli")

sim = simulate_ribosome_traffic(seq, ecoli, initiation_rate=1.0, seed=42)
print(f"Flux: {sim.protein_flux_per_min} proteins/min | NGD Risk: {sim.ngd_risk_score}/100")

report = compare_harmonization(seq, source=human, host=ecoli)
print("Harmonized Rhythm r:", report.harmonized.rhythm_correlation)
print("Naive CAI Rhythm r :", report.naive_optimized.rhythm_correlation)
```

---

## Running Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

---

## Scientific Background & References
- **$\ell$-TASEP Translation Physics:** Shaw, Zia & Lee (2003), *Totally asymmetric exclusion process with extended objects: a model for protein synthesis*, Phys. Rev. E.
- **%MinMax Folding Rhythm:** Clarke & Clark (2008), *Rare codons cluster*, PLoS ONE.
- **Codon Harmonization:** Angov et al. (2008), *Heterologous protein expression is enhanced by harmonizing the codon usage frequencies of the target gene with those of the expression host*, PLoS ONE.
- **Ribosome Collisions & NGD/RQC:** Simms, Yan & Zaher (2017), *Ribosome collision is critical for quality control during No-Go Decay*, Molecular Cell.

## License
MIT License — Copyright (c) 2026 Muhammed Emir Esti
