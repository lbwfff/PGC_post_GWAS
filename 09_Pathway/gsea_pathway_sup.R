
#对应主图2a的内容，为不同分析方法得到的RMP的数量统计

#我设想的2a是一个数量上的统计，MAGMA，HMAGMA，SMR，FUSION（DLPFC & blood），时间空间细胞类型的共定位
#全部放在一起

fig2a<-list()

#MAGMA

list<-list.files('MAGMA_singlegene/')
anno<-read.table('NCBI37.gene.loc',header = F)

filesave<-data.frame()

for (i in list){
  
  inf<-read.table(paste0('MAGMA_singlegene/',i),header = T) 
  
  disease<-unlist(strsplit(i,'_'))[2]
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-merge(inf,anno,by.x=c('GENE'),by.y=c('V1'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$V6 %in% hla$V6),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  gsea_results<-data.frame(
                    RMPsign=c(length(inf$V6[(inf$V6 %in% cache[["RMP"]]) & inf$fdr<0.05])),
                    NonRMPsign=c(length(inf$V6[!(inf$V6 %in% cache[["RMP"]]) & inf$fdr<0.05])))
  gsea_results$disorder<-disease
  
  
  filesave<-rbind(filesave,gsea_results)
}


filesave$methods<-c('MAGMA')
fig2a[['MAGMA']]<-filesave

#HMAGMA


list<-list.files('hMAGMA/hMAGMA_out/')
list<-list[grep('genes.out',list)]
load('hMAGMA/geneAnno_allgenes.rda')

filesave<-data.frame()

for (i in list){
  
  inf<-read.table(paste0('hMAGMA/hMAGMA_out/',i),header = T)
  disease<-unlist(strsplit(i,'_'))[2]
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  
  inf<-merge(inf,geneAnno1,by.x=c('GENE'),by.y=c('ensembl_gene_id'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$hgnc_symbol %in% hla$hgnc_symbol),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$hgnc_symbol[(inf$hgnc_symbol %in% cache[["RMP"]]) & inf$fdr<0.05])),
    NonRMPsign=c(length(inf$hgnc_symbol[!(inf$hgnc_symbol %in% cache[["RMP"]]) & inf$fdr<0.05])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave$methods<-c('HMAGMA')
fig2a[['HMAGMA']]<-filesave


#SMR

file<-list.files('eQTL_BrainMeta/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

filesave<-data.frame()

for (i in list){
  
  combsmr<-data.frame()
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  for (j in 1:22){
    inf<-read.table(paste0('eQTL_BrainMeta/PGC_',i,'_BrainMeta_cis_eQTL_chr',j,'.smr'),header = T)
    combsmr<-rbind(combsmr,inf)
  }
  
  combsmr<-combsmr[combsmr$p_HEIDI>0.05  & !is.na(combsmr$p_HEIDI),]
  combsmr$fdr<-p.adjust(combsmr$p_SMR,method = 'BH')
  
  inf<-combsmr
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$Gene),]
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$Gene[(inf$Gene %in% cache[["RMP"]]) & inf$fdr<0.05])),
    NonRMPsign=c(length(inf$Gene[!(inf$Gene %in% cache[["RMP"]]) & inf$fdr<0.05])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave$methods<-c('SMR_DPLFC')
fig2a[['SMR_DLPFC']]<-filesave

#


file<-list.files('TWAS_result/TWAS_result_GTEx_cortex/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

filesave<-data.frame()

for (i in list){
  combsmr<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('TWAS_result/TWAS_result_GTEx_cortex/PGC_',i,'_for_fusion_TWAS_',j,'.dat'),header = T)
    combsmr<-rbind(combsmr,inf)
  }
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  combsmr<-combsmr[!is.na(combsmr$TWAS.P),]
  combsmr$fdr<-p.adjust(combsmr$TWAS.P,method = 'BH')
  
  inf<-combsmr[,c(3,7,11,14:21)]
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$ID),]
  rownames(inf)<-inf$order
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$ID[(inf$ID %in% cache[["RMP"]]) & inf$fdr<0.05])),
    NonRMPsign=c(length(inf$ID[!(inf$ID %in% cache[["RMP"]]) & inf$fdr<0.05])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave$methods<-c('FUSION_DPLFC')
fig2a[['FUSION_DLPFC']]<-filesave

#

file<-list.files('TWAS_result/gtex_blood/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

filesave<-data.frame()

for (i in list){
  combsmr<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('TWAS_result/gtex_blood/PGC_',i,'_for_fusion_TWAS_',j,'.dat'),header = T)
    combsmr<-rbind(combsmr,inf)
  }
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  combsmr<-combsmr[!is.na(combsmr$TWAS.P),]
  combsmr$fdr<-p.adjust(combsmr$TWAS.P,method = 'BH')
  
  inf<-combsmr[,c(3,7,11,14:21)]
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$ID),]
  rownames(inf)<-inf$order
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$ID[(inf$ID %in% cache[["RMP"]]) & inf$fdr<0.05])),
    NonRMPsign=c(length(inf$ID[!(inf$ID %in% cache[["RMP"]]) & inf$fdr<0.05])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)

}

filesave$methods<-c('FUSION_blood')
fig2a[['FUSION_blood']]<-filesave

#然后就是三个property分析

file<-read.csv('property/dev_property_cor.csv')

filesave<-data.frame()

for (i in unique(file$disorder)){

  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-file[file$disorder==i,]
  inf$fdr<-p.adjust(inf$Pvalue,method = 'BH')
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$GENE),]
  rownames(inf)<-inf$order
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$GENE[(inf$GENE %in% cache[["RMP"]]) & inf$fdr<0.05 & abs(inf$Rvalue) > 0.6])),
    NonRMPsign=c(length(inf$GENE[!(inf$GENE %in% cache[["RMP"]]) & inf$fdr<0.05 & abs(inf$Rvalue) > 0.6])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave$methods<-c('BrainSpan_coloc')

fig2a[['BrainSpan_coloc']]<-filesave


file<-read.csv('property/region_property_cor.csv')

filesave<-data.frame()

for (i in unique(file$disorder)){
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-file[file$disorder==i,]
  inf$fdr<-p.adjust(inf$Pvalue,method = 'BH')
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$GENE),]
  rownames(inf)<-inf$order
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$GENE[(inf$GENE %in% cache[["RMP"]]) & inf$fdr<0.05 & abs(inf$Rvalue) > 0.6])),
    NonRMPsign=c(length(inf$GENE[!(inf$GENE %in% cache[["RMP"]]) & inf$fdr<0.05 & abs(inf$Rvalue) > 0.6])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave$methods<-c('BrainRegion_coloc')

fig2a[['BrainRegion_coloc']]<-filesave

file<-read.csv('property/celltype_property_cor.csv')

filesave<-data.frame()

for (i in unique(file$disorder)){
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-file[file$disorder==i,]
  inf$fdr<-p.adjust(inf$Pvalue,method = 'BH')
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$GENE),]
  rownames(inf)<-inf$order
  
  gsea_results<-data.frame(
    RMPsign=c(length(inf$GENE[(inf$GENE %in% cache[["RMP"]]) & inf$fdr<0.05 & abs(inf$Rvalue) > 0.6])),
    NonRMPsign=c(length(inf$GENE[!(inf$GENE %in% cache[["RMP"]]) & inf$fdr<0.05 & abs(inf$Rvalue) > 0.6])))
  gsea_results$disorder<-disease
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave$methods<-c('BrainCell_coloc')

fig2a[['BrainCell_coloc']]<-filesave

####

#先不考虑coloc的情况下看一下可视化

library(ggplot2)
library(tidyr)
library(dplyr)

p<-list()

for (i in 1:5) {

df_wide <- fig2a[[i]]
df_wide <- df_wide %>%
  mutate(label_text = paste0(NonRMPsign, " / ", RMPsign))
df_wide$label_text[df_wide$label_text=='0 / 0']<-c('')

df_long <- df_wide %>%
  pivot_longer(
    cols = c(RMPsign, NonRMPsign),        
    names_to = "Type",   
    values_to = "Value" )


df_long$disorder<-factor(df_long$disorder,levels = rev(gsub('MDD','PDD',rownames(heatmap_mat))))
df_long$Type<-factor(df_long$Type,levels = rev(c('NonRMPsign','RMPsign')))

if(i==1) {

p[[i]]<-
ggplot(df_long, aes(x = Value, y = disorder, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("NonRMPsign" = "#4A90E2", "RMPsign" = "#FF6B6B")) +
  geom_text(data = filter(df_wide, label_text != ""),
            aes(x = max(NonRMPsign + RMPsign) / 2, 
              y = disorder, 
              label = label_text),
            color = "black",
            hjust = 0.2, color = "black",  size = 3,   
            inherit.aes = FALSE) +
  theme_bw() +
  labs( x = "Number") +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())+
  theme(aspect.ratio=3) } else{
  
    p[[i]]<-
      ggplot(df_long, aes(x = Value, y = disorder, fill = Type)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_y_discrete(drop = FALSE)+
      scale_fill_manual(values = c("NonRMPsign" = "#4A90E2", "RMPsign" = "#FF6B6B")) +
      geom_text(data = filter(df_wide, label_text != ""),
                aes(x = max(NonRMPsign + RMPsign) / 2, 
                    y = disorder, 
                    label = label_text),
                color = "black",
                hjust = 0.2, color = "black",  size = 3,   
                inherit.aes = FALSE) +
      theme_bw() +
      labs( x = "Number") +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom",
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank())+
      theme(aspect.ratio=3)
    
  }

}


pdf("./property/GSEA_num.pdf", width = 12, height = 8)
patchwork::wrap_plots(p, nrow = 1,guides='collect')
dev.off()

#效果还行

