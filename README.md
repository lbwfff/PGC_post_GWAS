# Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia

This Repositories was used to include codes used in the project “Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia” (https://www.biorxiv.org/content/XXX).

Most of the code in this project is implemented using published software and existing R packages (including those from CRAN and Bioconductor).

## [Data sources](./Data_source/): 
Includes the RMP list and the sources of GWAS summary statistics utilized in this study.

## Gene mapping
[Gene_mapping.r](./Gene_mapping.r) annotated SNPs based on Gencode v26 (hg19 version) annotation and generated Manhattan plots and QQ plots for RMP genes. [zoom_local.r](./zoom_local.r) generated LocusZoom plots for genes exhibiting significant genomic associations.

## Fine mapping
[PolyFun.sh](./PolyFun.sh) performs SNP-level finemapping analysis based on PolyFun. Visualization is available at [zoom_local.r](./zoom_local.r).

## Gene-based analysis
MAGMA and hMAGMA analyses were performed using [MAGMA.sh](./MAGMA.sh) and [hMAGMA.sh](./hMAGMA.sh), respectively, with annotation files sourced from previously published projects. Visualization is performed by [Gene_based.r](./Gene_based.r).

## Gene-based analysis
The heritability estimation is implemented by [GCTB_sbayes.sh](GCTB_sbayes.sh), while [sbayes.R](sbayes.R) performs the aggregation and visualization of SNP-level estimates.

