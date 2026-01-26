#!/bin/bash
#SBATCH -p compute
#SBATCH -c 20
#SBATCH -t 20:00:00

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate

conda activate r

GWAS_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS"


for file in $GWAS_dir/*; do
    
    if [ -f "$file" ]; then  # Check if it's a regular file (not a directory)
        filename=$(basename "$file")
        echo "$filename"
    fi
    
    for chr in {1..22}; do

    Rscript /scratch/lb4489/project/GWAS/fusion_twas-master/FUSION.assoc_test.R \
        --sumstats /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/"$filename" \
        --weights /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/weights/CMC.BRAIN.RNASEQ.pos \
        --weights_dir /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/weights/ \
        --ref_ld_chr /scratch/lb4489/project/GWAS/fusion_twas-master/LDREF/1000G.EUR. \
        --chr $chr \
        --out /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/TWAS_result_brain/"${filename}_TWAS_${chr}.dat"

        echo "Processed: ${filename} for chromosome ${chr}"
    done
	
done


for file in $GWAS_dir/*; do
    
    if [ -f "$file" ]; then  # Check if it's a regular file (not a directory)
        filename=$(basename "$file")
        echo "$filename"
    fi
    
    for chr in {1..22}; do

    Rscript /scratch/lb4489/project/GWAS/fusion_twas-master/FUSION.assoc_test.R \
        --sumstats /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/"$filename" \
        --weights /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/weights/GTEx.Brain_Cortex.pos \
        --weights_dir /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/weights/ \
        --ref_ld_chr /scratch/lb4489/project/GWAS/fusion_twas-master/LDREF/1000G.EUR. \
        --chr $chr \
        --out /scratch/lb4489/project/GWAS/GWAS_from_PGC/TWAS/TWAS_result_GTEx_cortex/"${filename}_TWAS_${chr}.dat"

        echo "Processed: ${filename} for chromosome ${chr}"
    done
	
done