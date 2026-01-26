# Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia

This Repositories was used to include codes used in the project “Functional Genomic Analyses Reveal transfer RNA Modification Enzymes as Risk Genes for Bipolar I Disorder and Schizophrenia” (https://www.biorxiv.org/content/XXX).

Most of the code in this project is implemented by existing R packages (including Cran and bioconductor)

## STEP 0 Pre-processing of data 
The script [s0_data_prep.R](./s0_data_prep.R) is used for data preprocessing, we reformatted the NONCODE annotations to make them compatible as input for the Ribo-seq analysis tool and removed redundancies with the GENCODE annotations using the GffCompare tool.

## STEP1 Ribo-seq analysis
The [Ribo-seq folder](./Ribo-seq/) contains the command lines we used for the Ribo-seq analysis, including trim for reads, removal of rRNA, mapping, and prediction of activated ORFs. In [S1_RIbo_and_ORFprediction.R](./S1_RIbo_and_ORFprediction.R) Quality control was performed on the Ribo-seq.

In [s2_build_peptide_index.R](./s2_build_peptide_index.R) we constructed the indexes needed for proteomic analysis based on the results of Ribo-seq.

In [s3_ORF_Characterization.R](./s3_ORF_Characterization.R) we perform a preliminary characterization of the obtained ORFs

## STEP2 Computational Proteomics Analysis
We used the [Fragpipe](https://github.com/Nesvilab/FragPipe) platform to perform computational proteomics analysis, and the final parameters used can be seen in the fragpipe folder, which includes both sample annotation information and run logs.

[s4_MS_fdr_control.R](./s4_MS_fdr_control.R) is used to visualize the results of the FDR control of the [philosopher](https://github.com/Nesvilab/philosopher) and to count the number of mass spectrometry spectra obtained by different methods.

## STEP3 Estimation of protein expression
[s5_get_MS_expression.R](./s5_get_MS_expression.R) was used to make an estimate of protein expression levels, which was achieved through the use of [MSstats](https://bioconductor.org/packages/release/bioc/html/MSstats.html) and [MSstatsTMT](https://www.bioconductor.org/packages/release/bioc/html/MSstatsTMT.html)

## STEP4 peptide Characterization
In [s6_peptide_Characterization.R](./s6_peptide_Characterization.R) we trace the origin of the peptides, in addition we build a machine learning model for exploring why a large number of Ribo-seq peptides are not detected in protein mass spectrometry.

## STEP5 biostatistical analysis
Before that we used [s7_QC_for_protein_expression.R](./s7_QC_for_protein_expression.R) to process the expression matrices obtained from reanalyzing the [Jiang et al 2019 dataset](https://www.nature.com/articles/s41586-019-0987-8) and the [Gao et al 2019 dataset](https://www.sciencedirect.com/science/article/pii/S0092867419310037?via%3Dihub), and we checked for the presence of batch effects.

In [s8_biomark.R](./s8_biomark.R) we compared the predictive performance of models using canonical protein expression and lncRNA-derived peptides expression, and both were combined as biomarkers. These models included cancer tissue differentiation, survival prognosis prediction, and recurrence prognosis.

In [s9_biostatistic_for_peptide.R](./s9_biostatistic_for_peptide.R) we counted the biostatistical significance of separate lncRNA-derived peptides in the two datasets, and in [s10_TCGA.R](s10_TCGA.R) we discuss the lncRNAs from which these lncRNA-derived peptides originated in the TCGA-LIHC dataset for biostatistical significance.

In [s11_for_3peptide.R](./s11_for_3peptide.R) we visualized three representative lncRNA-derived peptides.



