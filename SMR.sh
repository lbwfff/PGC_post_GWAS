#!/bin/bash
#SBATCH -p compute
#SBATCH -c 10
#SBATCH --mem 50G
#SBATCH -t 48:00:00

GWAS_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/SMR_GWAS_input"


EQTL_dir="/scratch/lb4489/project/GWAS/GTEx"
OUTPUT_dir="/scratch/lb4489/project/GWAS/GWAS_from_PGC/SMR_result"

for gwas_file in $GWAS_dir/*; do 
    gwas_filename=$(basename "$gwas_file")
    
    for eqtl_file in $EQTL_dir/*.besd; do
        eqtl_filename=$(basename "$eqtl_file" .besd)
        
        echo "Processed: $gwas_filename with $eqtl_filename"

	    /scratch/lb4489/project/GWAS/smr-1.3.1-linux-x86_64/smr \
    	--bfile /scratch/lb4489/project/GWAS/csmr/data/1kg/EUR \
    	--maf 0.01 \
    	--gwas-summary "$gwas_file" \
    	--beqtl-summary "$EQTL_dir/$eqtl_filename" \
    	--out "$OUTPUT_dir/${gwas_filename}_$eqtl_filename" \
	    --thread-num 10 
    done
done

