#!/bin/bash
#SBATCH -p compute
#SBATCH -c 50
#SBATCH --mem 20G
#SBATCH -t 10:00:00

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_SCZ.ma \
  --impute-summary \
  --out PGC_SCZ_imputed

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --sbayes R --scale -1 \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_SCZ_imputed.imputed.ma \
  --chain-length 20000 \
  --burn-in 5000 \
  --out sbayesr_SCZ


/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_BIPI.ma \
  --impute-summary \
  --out PGC_BIPI_imputed

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --sbayes R --scale -1 \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_BIPI_imputed.imputed.ma \
  --chain-length 20000 \
  --burn-in 5000 \
  --out sbayesr_BIPI

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_AD.ma \
  --impute-summary \
  --out PGC_AD_imputed

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --sbayes R  \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_AD_imputed.imputed.ma \
  --chain-length 20000 \
  --burn-in 5000 \
  --out sbayesr_AD
  
/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_ADHD.ma \
  --impute-summary \
  --out PGC_ADHD_imputed

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --sbayes R  \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_ADHD_imputed.imputed.ma \
  --chain-length 20000 \
  --burn-in 5000 \
  --out sbayesr_ADHD
  
/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_ED.ma \
  --impute-summary \
  --out PGC_ED_imputed

/scratch/lb4489/project/GWAS/ldsc/gctb_2.5.4_Linux/gctb \
  --sbayes R --scale -1 \
  --ldm-eigen ukbEUR_HM3 \
  --gwas-summary PGC_ED_imputed.imputed.ma \
  --chain-length 20000 \
  --burn-in 5000 \
 --out sbayesr_ED
  
  
  