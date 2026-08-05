# Functional Genomic Analyses Reveal tRNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia

This repository contains the analysis code for the project:

> [**"Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia"**](https://www.medrxiv.org/content/10.64898/2026.01.26.26344826v1)

All analyses were implemented using published software and existing R packages (CRAN and Bioconductor).

---

## Table of Contents

- [Repository Structure](#repository-structure)
- [Analysis Steps](#analysis-steps)
  - [1. Gene Mapping](#1-gene-mapping)
  - [2. Fine Mapping](#2-fine-mapping)
  - [3. Gene-Based Analysis](#3-gene-based-analysis)
  - [4. Heritability Estimation](#4-heritability-estimation)
  - [5. SMR (Summary-data-based Mendelian Randomization)](#5-smr-summary-data-based-mendelian-randomization)
  - [6. TWAS (Transcriptome-Wide Association Study)](#6-twas-transcriptome-wide-association-study)
  - [7. Two-Sample MR & Colocalization](#7-two-sample-mr--colocalization)
  - [8. RNA-Seq Analysis](#8-rna-seq-analysis)
  - [9. Pathway Analysis](#9-pathway-analysis)
  - [10. Multi-Method Integration](#10-multi-method-integration)
  - [11. ColocDB Analysis](#11-colocdb-analysis)
- [Data Sources](#data-sources)
- [Requirements](#requirements)
- [Citation](#citation)

---

## Repository Structure

```
PGC_post_GWAS-main/
├── Data_source/           # GWAS summary statistics and RMP gene list
├── 01_Gene_mapping/       # SNP annotation and visualization
├── 02_Fine_mapping/       # SNP-level fine mapping
├── 03_Gene_based/         # Gene-based association tests
├── 04_Heritability/       # SNP heritability estimation
├── 05_SMR/               # Summary-data-based Mendelian Randomization
├── 06_TWAS/              # Transcriptome-Wide Association Study
├── 07_MR_coloc/          # Two-Sample MR & Colocalization
├── 08_RNA_seq/           # RNA-seq differential expression analysis
├── 09_Pathway/           # GSEA pathway enrichment analysis
├── 10_Integration/       # Multi-method sPLS-DA integration
├── 11_ColocDB/           # ColocDB database mining and re-analysis
└── README.md
```

---

## Analysis Steps

### 1. Gene Mapping

Annotates SNPs based on Gencode v26 (hg19) and generates Manhattan/QQ plots for RMP genes.

| File | Description |
|------|-------------|
| `01_Gene_mapping/Gene_mapping.r` | SNP annotation, Manhattan plots, and QQ plots |
| `01_Gene_mapping/zoom_local.r` | LocusZoom plots for significant genomic regions |

---

### 2. Fine Mapping

Performs SNP-level fine mapping using [PolyFun](https://github.com/omerwe/polyfun).

| File | Description |
|------|-------------|
| `02_Fine_mapping/PolyFun.sh` | SNP-level fine mapping pipeline |

---

### 3. Gene-Based Analysis

Performs gene-based association tests using [MAGMA](https://cncr.nl/research/magma/) and [hMAGMA](https://github.com/thewonlab/H-MAGMA).

| File | Description |
|------|-------------|
| `03_Gene_based/MAGMA.sh` | MAGMA gene-based analysis |
| `03_Gene_based/hMAGMA.sh` | hMAGMA cell-type-specific gene-based analysis |
| `03_Gene_based/Gene_based.r` | Visualization of gene-based results |

---

### 4. Heritability Estimation

Estimates SNP heritability and aggregates SNP-level estimates.

| File | Description |
|------|-------------|
| `04_Heritability/GCTB_sbayes.sh` | SBayesC heritability estimation via [GCTB](https://yanglab.westlake.edu.cn/software/gctb/) |
| `04_Heritability/sbayes_vis.R` | Aggregation and visualization of SNP-level estimates |

---

### 5. SMR (Summary-data-based Mendelian Randomization)

Integrates GWAS and eQTL/sQTL data to identify trait-associated genes.

| File | Description |
|------|-------------|
| `05_SMR/SMR.sh` | SMR analysis pipeline |
| `05_SMR/SMR.r` | Visualization of SMR results |

**Data:** Both cis-eQTL and cis-sQTL summary data are sourced from [BrainMeta v2](https://yanglab.westlake.edu.cn/software/smr/#sQTLsummarydata).

---

### 6. TWAS (Transcriptome-Wide Association Study)

Performs gene-level association tests using predicted expression weights.

| File | Description |
|------|-------------|
| `06_TWAS/TWAS.sh` | FUSION-based TWAS pipeline |
| `06_TWAS/FUSION.R` | Visualization of TWAS results |

**Data:** Precomputed weights are from [FUSION](http://gusevlab.org/projects/fusion/).

---

### 7. Two-Sample MR & Colocalization

Investigates causal relationships between gene expression and psychiatric traits.

| File | Description |
|------|-------------|
| `07_MR_coloc/MR_and_coloc.R` | Core functions for Two-Sample MR and colocalization (wrapping [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) and [coloc](https://github.com/chr1swallace/coloc)) |
| `07_MR_coloc/Run_TSMR_coloc.R` | Application of MR and colocalization to study data |

---

### 8. RNA-Seq Analysis

Differential expression analysis of postmortem brain samples.

| File | Description |
|------|-------------|
| `08_RNA_seq/BrainSeq_SCZ.r` | RNA-seq analysis for schizophrenia samples |
| `08_RNA_seq/BrainSeq_BP.r` | RNA-seq analysis for bipolar disorder samples |

**Data:** [BrainSeq Phase 2](https://eqtl.brainseq.org/phase2/)

---

### 9. Pathway Analysis

Gene set enrichment analysis (GSEA) and cross-method summary statistics for RMP genes.

| File | Description |
|------|-------------|
| `09_Pathway/gsea_pathway.R` | GSEA pathway enrichment analysis (Main Fig. 2b-2c) and sPLS-DA input preparation |
| `09_Pathway/gsea_pathway_sup.R` | Cross-method comparison of significant RMP counts (Main Fig. 2a) |

---

### 10. Multi-Method Integration

Sparse Partial Least Squares Discriminant Analysis (sPLS-DA) integrating results from multiple analytical methods.

| File | Description |
|------|-------------|
| `10_Integration/splsda_integration.R` | sPLS-DA integration of MAGMA, hMAGMA, SMR, TWAS, and colocalization results |

---

### 11. ColocDB Analysis

Mining and re-analysis of the [ColocDB](https://ngdc.cncb.ac.cn/colocdb/) database for schizophrenia and bipolar disorder colocalization evidence.

| File | Description |
|------|-------------|
| `11_ColocDB/batch_fetch_genes.py` | Python script to batch fetch COLOC and SMR data from ColocDB API |
| `11_ColocDB/colocdb_analysis.R` | Systematic analysis of ColocDB colocalization and SMR records |
| `11_ColocDB/smr_slope.R` | SMR effect slope analysis across genes (Main Fig. 6a) |

---

## Data Sources

See [`Data_source/`](./Data_source/) for detailed information:

- **GWAS summary statistics** (`GWAS_data.csv`) — Psychiatric disorder GWAS from PGC consortia
- **RMP gene list** (`RMP_list.csv`) — 123 RNA-modifying proteins curated from MODOMICS and RNAME

---

## Requirements

### Software
- [PLINK](https://www.cog-genomics.org/plink/)
- [GCTB](https://yanglab.westlake.edu.cn/software/gctb/)
- [MAGMA](https://cncr.nl/research/magma/)
- [SMR](https://yanglab.westlake.edu.cn/software/smr/)
- [FUSION](http://gusevlab.org/projects/fusion/)

### Python
- `requests`, `pandas`

### R Packages
- **CRAN:** `dplyr`, `ggplot2`, `ggrepel`, `patchwork`, `tidyr`, `tibble`, `fgsea`, `pcaMethods`, `RColorBrewer`
- **Bioconductor:** `GenomicRanges`
- **GitHub/Other:** `geni.plots`, `TwoSampleMR`, `coloc`, `ieugwasr`, `MetBrewer`, `cowplot`, `mixOmics`, `bbplot`

---

## Citation

If you use this code, please cite our paper:

```bibtex
@article{,
  title={Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia},
  year={2026},
  doi={10.64898/2026.01.26.26344826}
}
```
