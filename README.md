# RiboFlow v2.0

[![Live Demo](https://img.shields.io/badge/▶_Live_Demo-emiresti.github.io/riboflow-2dd4bf?style=for-the-badge&labelColor=0b1022)](https://emiresti.github.io/riboflow/) [![CI](https://github.com/emiresti/riboflow/actions/workflows/ci.yml/badge.svg)](https://github.com/emiresti/riboflow/actions/workflows/ci.yml)

**Biophysical Ribosome Traffic ($\ell$-TASEP), Turner 5' mRNA Folding Thermodynamics, Exit-Tunnel Electrostatic Drag & Multi-Objective mRNA Vaccine / SynBio Designer**

> Why do identical proteins encoded by different synonymous mRNAs differ by **100× in expression**?
> Because real translation is governed simultaneously by **5' mRNA hairpin accessibility ($\Delta G_{\text{5'}}$)**, **steric ribosome traffic collisions (Disomes / No-Go Decay)**, **nascent-peptide positive charge inside the ribosomal exit tunnel**, and **co-translational folding pauses (`%MinMax`)**.

**RiboFlow v2.0** unifies all four mechanisms in a **zero-dependency Python 3.10+ package, CLI, and live WebAssembly (Pyodide) browser laboratory**.

---

## What Makes RiboFlow v2.0 Unique?

### 1. Turner-Calibrated 5' mRNA Hairpin Thermodynamics ($\Delta G_{\text{5'}}$) Coupled Directly to $\ell$-TASEP
- Runs a dynamic programming Minimum Free Energy (MFE) RNA folder with Turner nearest-neighbor helix stacking & loop initiation energies over the critical first 45 nt (`AUG` initiation zone).
- Computes an interactive **Base-Pair Arc Diagram** + dot-bracket notation and couples $\Delta G_{\text{5'}}$ directly into the $\ell$-TASEP initiation rate ($\alpha_{\text{eff}} = \alpha_0 \cdot \text{Boltzmann}(\Delta G_{\text{5'}})$).
- Includes a 1-click **Synonymous 5' Hairpin Melter** that unfolds severe start-codon blockades (e.g., `-27.4 → 0.0 kcal/mol`, raising initiation accessibility from `15% → 99%`) while keeping the protein **100% identical**.

### 2. Exact Continuous-Time $\ell$-TASEP Simulator + Ribosome Exit-Tunnel Electrostatic Drag
- Gillespie kinetic Monte Carlo simulation of ribosomes (steric footprint $\ell = 10$ codons $\approx 30\text{ nt}$) across *Homo sapiens*, *E. coli* K-12, *S. cerevisiae*, and *M. musculus*.
- Models **positively charged nascent polypeptide braking (`K`/`R` polybasic clusters)** inside the negatively charged rRNA ribosomal exit tunnel (Charneski & Hurst, *PLOS Biology* 2013).
- Computes **protein flux ($J$)**, **simulated Ribo-seq & Disome-seq profiles**, **No-Go Decay (NGD / RQC) risk**, and **Statistical Physics Phase Transition Diagrams** ($\alpha$ sweep revealing the exact critical jamming threshold $\alpha_c$).

### 3. Multi-Objective SynBio & mRNA Vaccine ($m^1\Psi$ Uridine-Depletion) Designer
Simultaneously co-optimizes **5 wet-lab constraints** while guaranteeing **100% amino-acid identity**:
1. **Pause-Preserving Codon Harmonization** (maintains native `%MinMax` co-translational folding rhythm)
2. **5' Hairpin Unfolding** (opens the AUG initiation window)
3. **mRNA Vaccine Uridine-Depletion Mode** (reduces `U%` and `UU` dinucleotides to evade TLR7/8 innate immunity—e.g. drops SARS-CoV-2 Spike RBD uridine from `34.7%` → `16.2%` and `UU` pairs from `45` → `9`!)
4. **Cloning Restriction-Site Scrubbing** (eliminates *EcoRI*, *BamHI*, *HindIII*, *XhoI*, *NdeI*, and *BsaI* Golden Gate sites across codon boundaries)
5. **Homopolymer Filtering** (removes synthesis-breaking runs $\ge 5\text{ nt}$)

### 4. Interactive Browser Lab (Live NCBI RefSeq Import + PNG & Printable Lab Reports)
- Try immediately in your browser at **https://emiresti.github.io/riboflow/** (zero installation via Pyodide WebAssembly, **TR / EN** toggle).
- Fetch any real RefSeq gene live from **NCBI E-utilities** (`INS`, `TP53`, `CFTR`, `BRCA1`, `NM_000207`), export high-res **PNG** figures for papers/presentations, or print a complete **Biophysical Lab Report**.

---

## Quick Start (Local Python / CLI)

Requires **Python 3.10+** and **zero external runtime dependencies**.

```bash
git clone https://github.com/emiresti/riboflow.git
cd riboflow
pip install -e .

# Launch the local interactive v2.0 laboratory
riboflow serve --open
```

## CLI Examples

```bash
# Fold and synonymously melt an inhibitory 5' mRNA hairpin
riboflow fold --sample polybasic_arrest --melt

# Run l-TASEP with 5' hairpin coupling & exit-tunnel charge drag
riboflow simulate --sample polybasic_arrest --organism ecoli

# Run mRNA vaccine Uridine-depletion + cloning-site scrubbing + pause preservation
riboflow synbio --sample sars_cov2_rbd --source human --host human --vaccine-u-deplete

# Sweep initiation rate alpha to map the TASEP phase transition & critical alpha_c
riboflow phase --sample jam_benchmark --organism ecoli

# Scan for the most dangerous single-nucleotide synonymous ("silent") mutations
riboflow mutate --sample cftr_nbd1 --organism ecoli --scan
```

## Running Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Scientific References
- **5' mRNA Folding & Expression:** Kudla et al. (2009), *Coding-sequence determinants of gene expression in Escherichia coli*, *Science*; Goodman et al. (2013), *Science*.
- **$\ell$-TASEP Translation Physics:** Shaw, Zia & Lee (2003), *Totally asymmetric exclusion process with extended objects*, *Phys. Rev. E*.
- **Exit-Tunnel Electrostatic Drag:** Charneski & Hurst (2013), *Positively charged residues are the major determinants of ribosomal velocity*, *PLOS Biology*.
- **%MinMax & Codon Harmonization:** Clarke & Clark (2008), *PLOS ONE*; Angov et al. (2008), *PLOS ONE*.
- **Ribosome Collisions & NGD/RQC:** Simms, Yan & Zaher (2017), *Molecular Cell*.
- **mRNA Vaccine Uridine Depletion:** Karikó et al. (2005), *Immunity*; Thess et al. (2015), *Molecular Therapy*.

## License
MIT License — Copyright (c) 2026 Muhammed Emir Esti
