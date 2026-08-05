#!/bin/bash
#SBATCH -p compute
#SBATCH -c 10
#SBATCH -t 4:00:00

source /share/apps/NYUAD5/miniconda/3-4.11.0/bin/activate
conda activate polyfun

mkdir -p output

python ./polyfun/extract_snpvar.py --sumstats ./PGC3_for_polyfun.txt --out output/snps_with_var.gz --allow-missing

python ./polyfun/extract_snpvar.py --sumstats ./PGC_BIPI_for_polyfun.txt --out output/BIPI_with_var.gz --allow-missing
python ./polyfun/extract_snpvar.py --sumstats ./PGC_SCZ_for_polyfun.txt --out output/SCZ_with_var.gz --allow-missing
python ./polyfun/extract_snpvar.py --sumstats ./PGC_AD_for_polyfun.txt --out output/AD_with_var.gz --allow-missing
python ./polyfun/extract_snpvar.py --sumstats ./PGC_ADHD_for_polyfun.txt --out output/ADHD_with_var.gz --allow-missing
python ./polyfun/extract_snpvar.py --sumstats ./PGC_ED_for_polyfun.txt --out output/ED_with_var.gz --allow-missing

cat output/snps_with_var.gz | zcat | head

python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr2_198000001_201000001 --sumstats ./output/SCZ_with_var.gz --n 130644 --chr 2 --start 200544698 --end 201072460 --method susie --max-num-causal 2 --out output/TYW5_SCZ_UKB.gz
python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr3_51000001_54000001 --sumstats ./output/SCZ_with_var.gz --n 130644 --chr 3 --start 53067447 --end 53633655 --method susie --max-num-causal 2 --out output/DCP1A_SCZ_UKB.gz
python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr10_18000001_21000001 --sumstats ./output/SCZ_with_var.gz --n 130644 --chr 10 --start 18584490 --end 19192552 --method susie --max-num-causal 2 --out output/NSUN6_SCZ_UKB.gz
python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr14_102000001_105000001 --sumstats ./output/SCZ_with_var.gz --n 130644 --chr 14 --start 103743521 --end 104253411 --method susie --max-num-causal 2 --out output/TRMT61A_SCZ_UKB.gz
python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr5_6000001_9000001 --sumstats ./output/BIPI_with_var.gz --n 475038 --chr 5 --start 6349352 --end 6885405 --method susie --max-num-causal 2 --out output/NSUN2_BIPI_UKB.gz
python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr17_33000001_36000001 --sumstats ./output/BIPI_with_var.gz --n 475038 --chr 17 --start 34706001 --end 35215408 --method susie --max-num-causal 2 --out output/MRM1_BIPI_UKB.gz
python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/chr19_9000001_12000001 --sumstats ./output/BIPI_with_var.gz --n 475038 --chr 19 --start 9000001 --end 11074114 --method susie --max-num-causal 2 --out output/QTRT1_BIPI_UKB.gz
    

    