# BigNemEvol

# PanNemaEvol: Nematode Macroevolutionary & Phylogenomic Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Data: Zenodo](https://img.shields.io/badge/Data-Zenodo-blue.svg)](https://zenodo.org/)

Official code repository for the macroevolutionary and phylogenomic analyses described in our manuscript evaluating the deep evolutionary fluidity and reversibility of parasitism across the phylum Nematoda.

## 📖 Abstract Summary
Nematodes exhibit remarkable macroevolutionary plasticity, adapting to both free-living and parasitic niches. Challenging classical paradigms of biological irreversibility (Dollo’s Law), our analyses reveal that obligate parasitism was not a terminal specialization, but a highly dynamic and reversible ecological state. 

By comprehensively mining ~6 million public accessions and assembling a contiguous supermatrix of full ribosomal operons (18S–ITS1–5.8S–ITS2–28S) for 3,654 unique species, we deployed stringent topological validation and hidden-state Markov models (corHMM). Our results demonstrate pervasive secondary reversions to free-living states and show that the basal emergence of parasitism completely decoupled genome size from molecular evolutionary rates.

## 🗂 Repository Structure

```text
├── 01_Data_Acquisition/     # Python scripts for NCBI/SILVA bulk mining and filtering
├── 02_Operon_Assembly/      # Bioinformatic pipeline for de novo ribosomal operon reconstruction
├── 03_Phylogenetics/        # Scripts for MAFFT alignment, ClipKit trimming, and IQ-TREE 2
├── 04_Macroevolution/       # R scripts for PCMs (PGLS, corHMM, SIMMAP, Root-Sweep)
├── data/                    # Processed metadata and curated trait matrices (.csv)
├── trees/                   # Final phylograms and timetrees (.nwk)
└── README.md                # Project documentation
```

## ⚙️ Prerequisites and Dependencies

### Bioinformatics Tools
Ensure the following tools are installed and available in your system `$PATH`:
* **VSEARCH** (Dereplication and UCHIME3 chimera removal)
* **SeqKit** (Sequence standardization)
* **MAFFT v7.505** (Multiple sequence alignment)
* **ClipKit** (Alignment trimming using `-m smart-gap`)
* **IQ-TREE 2 v2.1.3** (Maximum Likelihood phylogenetic inference)
* **FastTree 2** (Initial topological screening)

### Python Environment (>= 3.9)
* `Biopython`, `pandas`, `DendroPy`, `ETE3` (for tree processing and visualization)

### R Environment (>= 4.3)
Phylogenetic Comparative Methods (PCMs) rely on the following R packages:
* `ape`, `phytools`, `geiger`, `nlme`, `corHMM`, `diversitree`, `phylolm`, `coda`

## 🚀 Pipeline Workflow

### 1. Data Acquisition and Operon Reconstruction (`01_Data_Acquisition` & `02_Operon_Assembly`)
Our custom Python pipeline accesses the NCBI E-utilities API and the SILVA database to retrieve raw nucleotide records. Sequences undergo rigorous multi-step curation:
* Standardization and chimera removal (`--uchime3_denovo`).
* Metadata-anchored *de novo* assembly merging 18S, ITS, and 28S fragments into continuous operons.
* Patristic-distance-based topological screening (dynamic thresholding) to eliminate rogue taxa and deep misidentifications.

### 2. Phylogenetic Inference (`03_Phylogenetics`)
Alignments are processed using MAFFT and optimized via ClipKit. Phylogenetic reconstruction is performed using **IQ-TREE 2**:
* **Model:** GTR+F+G (identified via ModelFinder using BIC).
* **Support:** 500 UFBoot replicates, 1000 SH-aLRT replicates, and aBayes tests.
* Only robust nodes meeting strict thresholds (UFBoot ≥95%, SH-aLRT ≥80%, aBayes ≥0.90) are retained for downstream evolutionary analysis.

### 3. Macroevolutionary Modeling (`04_Macroevolution`)
The R scripts encompass the core statistical framework used to test Dollo's Law and analyze trait evolution:
* **Root-Sweep Sensitivity Analysis:** Dynamically evaluates ancestral state polarity across alternative deep nodes to bypass Midpoint/Outgroup topological biases.
* **Ancestral State Reconstruction (ASR):** Utilizes Equal-Rates (ER), All-Rates-Different (ARD), irreversible (Dollo), and Hidden State Markov Models (`corHMM`).
* **Stochastic Character Mapping (SIMMAP):** Quantifies the probabilistic history and frequency of lifestyle transitions (reversions).
* **Continuous Trait Evolution:** Compares Brownian Motion (BM) vs. Ornstein-Uhlenbeck (OU) models for genome size evolution.
* **Phylogenetic Regression:** Uses robust PGLS (incorporating `varIdent` structures) to evaluate the decoupling of genome size from mutational turnover.
* **Diversification Dynamics:** Calculates Pybus & Harvey’s γ statistic and compares Yule vs. Birth-Death models via `diversitree`.



## 📝 Citation
If you use the code or data in this repository, please cite our paper:
> **[Author Names]** (2026). *[Paper Title]*. Nature Genetics. DOI: [Insert DOI]

## 📜 License
This project is licensed under the MIT License
