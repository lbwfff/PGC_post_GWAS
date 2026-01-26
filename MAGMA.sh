#!/bin/bash
#SBATCH -p compute
#SBATCH -c 5
#SBATCH --mem 20000
#SBATCH -t 08:00:00

GWAS_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input"
out_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_out"

for file in $GWAS_dir/*; do

    filename=$(basename "$file")

    awk 'BEGIN {FS=OFS="\t"} NR>1 {print $1, $2, $3}' "$GWAS_dir/$filename" > "$out_dir/${filename}snp-loc"

    /scratch/lb4489/project/GWAS/MAGMA/magma  --annotate --snp-loc "$out_dir/${filename}snp-loc" --gene-loc /scratch/lb4489/project/GWAS/MAGMA/NCBI37.gene.loc --out "$out_dir/${filename}_magma"

    /scratch/lb4489/project/GWAS/MAGMA/magma \
        --bfile /scratch/lb4489/project/GWAS/MAGMA/g1000_eur \
        --gene-annot "$out_dir/${filename}_magma.genes.annot" \
        --pval "$GWAS_dir/$filename" ncol='NOBS' \
        --gene-model snp-wise=mean \
        --out "$out_dir/${filename}_singlegene"

    geneset=/scratch/lb4489/project/GWAS/GWAS_from_PGC/RNAMEgeneset.txt

    /scratch/lb4489/project/GWAS/MAGMA/magma \
        --gene-results "$out_dir/${filename}_singlegene.genes.raw" \
        --set-annot ${geneset} \
        --out "$out_dir/${filename}_geneset_RMRgene"

    echo "Processed: $filename"
	
done
