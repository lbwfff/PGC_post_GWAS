# Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia

This Repositories was used to include codes used in the project “Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia” (https://www.biorxiv.org/content/XXX).

Most of the code in this project is implemented using published software and existing R packages (including those from CRAN and Bioconductor).

## [Data sources](./Data_source/): 
Includes the RMP list and the sources of GWAS summary statistics utilized in this study.

## STEP 0 Pre-processing of data 
The script [s0_data_prep.R](./s0_data_prep.R) is used for data preprocessing, we reformatted the NONCODE annotations to make them compatible as input for the Ribo-seq analysis tool and removed redundancies with the GENCODE annotations using the GffCompare tool.

## STEP1 Ribo-seq analysis
The [Ribo-seq folder](./Ribo-seq/) contains the command lines we used for the Ribo-seq analysis, including trim for reads, removal of rRNA, mapping, and prediction of activated ORFs. In [S1_RIbo_and_ORFprediction.R](./S1_RIbo_and_ORFprediction.R) Quality control was performed on the Ribo-seq.

In [s2_build_peptide_index.R](./s2_build_peptide_index.R) we constructed the indexes needed for proteomic analysis based on the results of Ribo-seq.

In [s3_ORF_Characterization.R](./s3_ORF_Characterization.R) we perform a preliminary characterization of the obtained ORFs


In [s11_for_3peptide.R](./s11_for_3peptide.R) we visualized three representative lncRNA-derived peptides.



