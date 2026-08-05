#!/bin/bash
#SBATCH -p compute
#SBATCH -c 10
#SBATCH --mem 20000
#SBATCH -t 24:00:00

GWAS_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input"
out_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/hMAGMA_out"

for file in $GWAS_dir/*; do
    filename=$(basename "$file")

    /scratch/lb4489/project/GWAS/MAGMA/magma \
        --bfile /scratch/lb4489/project/GWAS/MAGMA/g1000_eur \
        --gene-annot /scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/Annotation_Files/Adultbrain.transcript.annot \
        --pval "$GWAS_dir/$filename" ncol='NOBS' \
        --gene-model snp-wise=mean \
        --out "$out_dir/${filename}_singlegene_hmagma"
    
    
    geneset=/scratch/lb4489/project/GWAS/GWAS_from_PGC/RNAMEgeneset_hMAGMA.txt

    /scratch/lb4489/project/GWAS/MAGMA/magma \
        --gene-results "$out_dir/${filename}_singlegene_hmagma.genes.raw" \
        --set-annot ${geneset} \
        --out "$out_dir/${filename}_geneset_RMRgene"
    
    echo "Processed: $filename"
	
done