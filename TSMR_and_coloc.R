
#还是可以试一下还是进行finemap的coloc

source('MR_and_coloc.R')

all_outcome<-data.table::fread('SMR_plot/PGC_BIPI')
all_outcome<-data.frame(rsids=c(all_outcome$SNP),
                        af_alt=c(all_outcome$Freq))

result <- MR_and_coloc('./SMR_plot/NSUN2_BIPI.ENSG00000037474.15.txt','BIPI','./SMR_plot/',
                       2865,420000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_BIPI.ENSG00000037474.15.txt')
write.csv(readsmr[["GWAS"]],file = 'NSUN2_BIPI_SNP.csv')

source('SMRlocus.R')
readsmr<-ReadSMRData('./SMR_plot/NSUN2_SCZ.ENSG00000037474.15.txt')
SMREffectPlot_old(readsmr)
SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_BIPI.ENSG00000037474.15.txt')
SMREffectPlot_old(readsmr)

all_outcome<-data.table::fread('SMR_plot/PGC_SCZ')
all_outcome<-data.frame(rsids=c(all_outcome$SNP),
                        af_alt=c(all_outcome$Freq))

result <- MR_and_coloc('./SMR_plot/NSUN2_SCZ.ENSG00000037474.15.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6) #每个snp的样本数都不太一样怎么说


readsmr<-ReadSMRData('./SMR_plot/NSUN2_SCZ.ENSG00000037474.15.txt')
write.csv(readsmr[["GWAS"]],file = 'NSUN2_SCZ_SNP.csv')


result <- MR_and_coloc('./SMR_plot/NAT10_SCZ.ENSG00000104408.11.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/EIF3E_SCZ.ENSG00000135372.9.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)


all_outcome<-data.table::fread('SMR_plot/PGC_SCZ')
all_outcome<-data.frame(rsids=c(all_outcome$SNP),
                        af_alt=c(all_outcome$Freq))

result <- MR_and_coloc('./SMR_plot/TRMT61B_SCZ_sQTL.chr2_29084174_29084778_clu_89159_.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./myplotTRMT61A.ENSG00000166166.13.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

#TRMT61B倒是可以说一下，其他的就算了吧


source('MR_and_coloc.R')

all_outcome<-data.table::fread('SMR_plot/PGC_BIPI')
all_outcome<-data.frame(rsids=c(all_outcome$SNP),
                        af_alt=c(all_outcome$Freq))

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL_BIP.chr5_6632869_6633744_clu_41364_.txt','BIPI','./SMR_plot',
                       2865,420000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL2_BIP.chr5_6605488_6606933_clu_41360_.txt','BIPI','./SMR_plot',
                       2865,420000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL3_BIP.chr5_6620411_6623327_clu_41362_.txt','BIPI','./SMR_plot',
                       2865,420000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL3_BIP.chr5_6621857_6623327_clu_41362_.txt','BIPI','./SMR_plot',
                       2865,420000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)


all_outcome<-data.table::fread('SMR_plot/PGC_SCZ')
all_outcome<-data.frame(rsids=c(all_outcome$SNP),
                        af_alt=c(all_outcome$Freq))

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL_SCZ.chr5_6621857_6623327_clu_41362_.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL2_SCZ.chr5_6632869_6633744_clu_41364_.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/NSUN2_sQTL3_SCZ.chr5_6605488_6606933_clu_41360_.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/EIF3F_sQTL_SCZ.chr8_109245901_109247227_clu_69235_.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)

result <- MR_and_coloc('./SMR_plot/EIF3F_sQTL2_SCZ.chr8_109241424_109247227_clu_69235_.txt','SCZ','./SMR_plot',
                       2865,120000,37,all_outcome,'../PTSD/clump/1kg.v3/EUR',clump_r2=0.2,exposure_p=10^-6)


####

source ('SMRlocus.R')

p<-list()

SMRData = ReadSMRData("SMR_plot/NSUN2_sQTL_SCZ.chr5_6621857_6623327_clu_41362_.txt")

p[[1]]<-SMREffectPlot(data=SMRData, trait_name="SCZ") 

SMRData = ReadSMRData("SMR_plot/NSUN2_sQTL2_SCZ.chr5_6632869_6633744_clu_41364_.txt")

p[[2]]<-SMREffectPlot(data=SMRData, trait_name="SCZ")

SMRData = ReadSMRData("SMR_plot/NSUN2_sQTL3_SCZ.chr5_6605488_6606933_clu_41360_.txt")

p[[3]]<-SMREffectPlot(data=SMRData, trait_name="BIPI") 


SMRData = ReadSMRData("SMR_plot/NSUN2_sQTL_BIP.chr5_6632869_6633744_clu_41364_.txt")

p[[4]]<-SMREffectPlot(data=SMRData, trait_name="SCZ") 

SMRData = ReadSMRData("SMR_plot/NSUN2_sQTL2_BIP.chr5_6605488_6606933_clu_41360_.txt")

p[[5]]<-SMREffectPlot(data=SMRData, trait_name="SCZ") 

SMRData = ReadSMRData("SMR_plot/NSUN2_sQTL3_BIP.chr5_6621857_6623327_clu_41362_.txt")

p[[6]]<-SMREffectPlot(data=SMRData, trait_name="SCZ") 

pdf("SMR_effect_plot_NSUN2splicing.pdf",width = 12,height = 8)
patchwork::wrap_plots(p,nrow=2) 
dev.off()


SMRData = ReadSMRData("SMR_plot/TRMT61B_SCZ_sQTL.chr2_29084174_29084778_clu_89159_.txt")

p<-SMREffectPlot(data=SMRData, trait_name="SCZ") 

pdf("SMR_TRMT61B_SCZ.pdf",width = 5,height = 4)
p
dev.off()



