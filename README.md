# Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia

This Repositories was used to include codes used in the project “Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia” (https://www.biorxiv.org/content/XXX).

Most of the code in this project is implemented using published software and existing R packages (including those from CRAN and Bioconductor).

## [Data sources](./Data_source/): 
Includes the RMP list and the sources of GWAS summary statistics utilized in this study.

## Gene mapping
[Gene_mapping.r](./Gene_mapping.r) annotated SNPs based on Gencode v26 (hg19 version) annotation and generated Manhattan plots and QQ plots for RMP genes. [zoom_local.r](./zoom_local.r) generated LocusZoom plots for genes exhibiting significant genomic associations.

## Fine mapping
[PolyFun.sh](./PolyFun.sh) performs SNP-level finemapping analysis based on [PolyFun](https://github.com/omerwe/polyfun). Visualization is available at [zoom_local.r](./zoom_local.r).

## Gene-based analysis
[MAGMA](https://cncr.nl/research/magma/) and [hMAGMA](https://github.com/thewonlab/H-MAGMA) analyses were performed using [MAGMA.sh](./MAGMA.sh) and [hMAGMA.sh](./hMAGMA.sh), respectively, with annotation files sourced from previously published projects. Visualization is performed by [Gene_based.r](./Gene_based.r).

## Gene-based analysis
The heritability estimation is implemented by [GCTB_sbayes.sh](GCTB_sbayes.sh), while [sbayes.R](sbayes.R) performs the aggregation and visualization of SNP-level estimates.

## SMR
[SMR](https://yanglab.westlake.edu.cn/software/smr/#Overview) analysis are implemented by [SMR.sh](./SMR.sh), with visualization performed by [SMR.r](./SMR.r).  Both cis-eQTL and cis-sQTL summary data are sourced from [BrainMeta v2](https://yanglab.westlake.edu.cn/software/smr/#sQTLsummarydata).

## FUSION
[FUSION-based TWAS](http://gusevlab.org/projects/fusion/) analysis are implemented by [FUSION.sh](./FUSION.sh), with visualization performed by [FUSION.r](./FUSION.r). The precomputed weights are provided by FUSION.

## TSMR & coloc
TwoSample-MR and gene colocalization analysis are based on the customized code [MR_and_coloc.R](MR_and_coloc.R), which primarily implements functionality by calling the [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR) and [coloc](https://github.com/chr1swallace/coloc) packages. [TSMR_and_coloc.R](TSMR_and_coloc.R) applies this customized code to the project.

## RNA-seq
[BrainSeq_SCZ.r](./BrainSeq_SCZ.r) and [BrainSeq_BP.r](./BrainSeq_BP.r) analyzed RNA-seq data from schizophrenia and bipolar disorder samples in the [BrainSeq](https://eqtl.brainseq.org/phase2/) database, respectively.

