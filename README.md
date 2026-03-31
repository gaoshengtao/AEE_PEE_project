# AEE_PEE_project

> Multi-omics analysis framework for acetate- and propionate-enriched rumen fermentation ecotypes in cattle.

![Platform](https://img.shields.io/badge/platform-R%20%7C%20Bash-blue)
![Project](https://img.shields.io/badge/project-rumen%20multi--omics-green)
![Status](https://img.shields.io/badge/status-active-success)

---

## Overview

This repository contains analysis scripts for a multi-omics study of **VFA-defined rumen fermentation ecotypes** in cattle, with a particular focus on the comparison between:

- **AEE**: *Acetate-Enriched Ecotype*
- **PEE**: *Propionate-Enriched Ecotype*

The project integrates multiple analytical layers, including:

- VFA-based rumenotype classification
- untargeted metabolomics
- host GWAS
- metagenome assembly
- MAG reconstruction and taxonomic annotation
- MAG abundance and functional profiling
- gene catalog functional analysis
- multi-layer integrative analyses

The overall goal is to identify ecological, functional, metabolic, and host-genetic differences associated with acetate- versus propionate-enriched rumen fermentation patterns.

---

## Repository structure

```text
AEE_PEE_project/
├─ 00_rumenotype_vfa/
├─ 01_metabolomics/
├─ 02_host_gwas/
├─ 03_metagenome_assembly/
├─ 04_mag_reconstruction_taxonomy/
├─ 05_mag_abundance_function/
├─ 06_gene_catalog_function/
└─ 07_integrative_analysis/
```

---

## Module description

### `00_rumenotype_vfa/`
Scripts for VFA profile processing, Aitchison/CLR-based clustering, and assignment of rumen fermentation ecotypes (AEE/PEE).

**Typical tasks**
- VFA feature construction
- compositional transformation
- clustering and evaluation
- rumenotype assignment export

### `01_metabolomics/`
Scripts for untargeted metabolomics processing and downstream analyses.

**Typical tasks**
- raw MS1/MS2 processing
- feature extraction and annotation
- OPLS-DA and differential metabolite analysis
- pathway-oriented metabolite visualization
- peptide and amino acid subset analysis

### `02_host_gwas/`
Scripts for host genomic association analyses.

**Typical tasks**
- host read QC and alignment
- BAM generation and indexing
- AP-ratio GWAS
- MAG-trait GWAS
- suggestive SNP extraction
- post-GWAS processing
- heritability estimation
- bidirectional MR analyses

### `03_metagenome_assembly/`
Scripts for metagenomic read processing and assembly workflows.

**Typical tasks**
- read trimming
- primary assembly
- secondary assembly
- scaffold-based read remapping
- contig filtering and dereplication

### `04_mag_reconstruction_taxonomy/`
Scripts for MAG quality control, dereplication, taxonomy assignment, phylogeny, and similarity analysis.

**Typical tasks**
- dRep-based dereplication
- CheckM/CheckM2 reruns
- GTDB-Tk annotation
- GTDB-based tree preparation
- PhyloPhlAn tree construction
- NR-based protein similarity analysis

### `05_mag_abundance_function/`
Scripts for MAG abundance profiling, differential abundance analysis, and MAG functional analysis.

**Typical tasks**
- CoverM-based abundance processing
- LEfSe input generation
- alpha diversity analysis
- MAG differential abundance visualization
- functional matrix correction
- dimensional reduction and clustering
- KEGG module enrichment and summary plots

### `06_gene_catalog_function/`
Scripts for gene catalog construction and function-level profiling.

**Typical tasks**
- ORF calling
- non-redundant gene catalog generation
- eggNOG annotation
- kallisto quantification
- KO/module-level summarization and downstream analyses

### `07_integrative_analysis/`
Scripts for cross-layer integrative analyses, especially those linking rumenotype, MAG clusters, metabolite branches, and host-genetic signals.

**Typical tasks**
- differential MAG distribution across clusters
- multivariate testing in reduced-dimensional space
- cross-omics interpretation
- figure-oriented integrative summaries

---

## Analytical workflow

The project follows a multi-step workflow:

1. Define rumenotypes using VFA composition profiles
2. Compare metabolite landscapes between AEE and PEE
3. Assemble metagenomes and reconstruct representative MAGs
4. Annotate MAG taxonomy and function
5. Quantify MAG abundance and identify AEE/PEE-associated MAGs
6. Perform host GWAS for AP ratio and MAG traits
7. Integrate metabolomic, microbiome, and host-genetic results

---

## Data and file conventions

### Input data

This repository is designed to work with the following types of input data:

- VFA profile tables
- untargeted metabolomics matrices
- host genomic BAM/VCF/PLINK files
- metagenomic paired-end reads
- MAG fasta files
- MAG abundance tables
- annotation resources such as GTDB, eggNOG, and KEGG-related tables

### Important note

Many scripts currently contain **hard-coded absolute paths** from the development environment.  
Before running the scripts on another machine, please update:

- input file paths
- output directories
- software executable paths
- database locations
- thread and memory parameters

---

## Environment and dependencies

This project was developed using a mixture of **R** and **Bash** workflows.

### Core languages

- R
- Bash

### Common R packages used across modules

- tidyverse
- data.table
- readxl
- openxlsx
- ggplot2
- vegan
- cluster
- limma
- DESeq2
- ropls
- patchwork
- ggtree
- ggtreeExtra
- dbscan
- uwot
- Rtsne

### Common command-line tools

- samtools
- bcftools
- plink2
- gcta
- bwa
- bowtie2
- Trimmomatic
- MEGAHIT
- metaSPAdes
- cd-hit-est
- dRep
- GTDB-Tk
- PhyloPhlAn
- Prodigal
- eggNOG-mapper
- kallisto

---

## Getting started

### 1. Clone the repository

```bash
git clone https://github.com/gaoshengtao/AEE_PEE_project.git
cd AEE_PEE_project
```

### 2. Review script paths

Before execution, open each script and modify the relevant file paths and software/database locations.

### 3. Run analyses by module

A typical order is:

```text
00_rumenotype_vfa
→ 01_metabolomics
→ 03_metagenome_assembly
→ 04_mag_reconstruction_taxonomy
→ 05_mag_abundance_function
→ 06_gene_catalog_function
→ 02_host_gwas
→ 07_integrative_analysis
```

Depending on your analysis purpose, modules can also be run independently.

---

## Recommended usage strategy

Because this repository contains both exploratory and production-style scripts, it is recommended to:

- treat each folder as an analysis module
- read script headers before execution
- standardize paths before batch running
- save intermediate outputs for reproducibility
- document any local modifications to parameters or databases

---

## Reproducibility notes

To improve portability and reproducibility, future versions of this repository may include:

- unified configuration files
- Conda environment files
- standardized input/output naming rules
- wrapper scripts for one-click execution
- example test datasets

At present, the repository mainly serves as a **project-organized analysis codebase** rather than a fully packaged pipeline.

---

## Outputs

Typical outputs generated by this project include:

- rumenotype assignment tables
- differential metabolite result tables
- OPLS-DA and volcano plots
- metagenome assembly outputs
- dereplicated MAG sets
- MAG phylogenetic trees
- LEfSe and abundance plots
- KO/module enrichment results
- AP and MAG GWAS result tables
- MR summary tables and visualization figures
- integrative summary figures for manuscript preparation

---

## Intended use

This repository is intended for:

- project-level data organization
- reproducible multi-omics analysis
- figure generation for manuscripts
- methodological reference for similar rumen microbiome studies

---

## Citation

If you use this repository or adapt parts of its workflow, please cite the corresponding study when available.

**Citation information will be added after manuscript publication.**

---

## License

This repository currently does not include an explicit license.

If you would like others to reuse the code with clear permissions, consider adding one of the following:

- MIT License
- BSD 3-Clause License
- GPL-3.0 License

---

## Contact

**Author:** Shengtao Gao  
**GitHub:** [gaoshengtao](https://github.com/gaoshengtao)

For questions, collaboration, or code-related issues, please open an issue in this repository or contact the author directly.
