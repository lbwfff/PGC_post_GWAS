
#包括GSEA的分析，对应主图2b-2c
#还包括一部分为splsda的input部分的整理

library(fgsea)
####geneset

modlist<-read.csv('geneset2.csv')
modlist<-modlist[modlist$Gene.Symbol !='',]
unique(modlist$Gene.Symbol)
cache<-list()
cache[['RMP']]<-c(modlist$Gene.Symbol)

motype<-paste0(unique(modlist$Modification),collapse = ',')
motype<-unlist(strsplit(motype,'[,]'))
motype<-unique(motype)

for (i in motype){
  
  inf<-modlist[grep(i,modlist$Modification),]
  cache[[i]]<-c(inf$Gene.Symbol)
  
}

motype<-paste0(unique(modlist$RNA),collapse = ',')
motype<-unlist(strsplit(motype,'[,]'))
motype<-unique(motype)

for (i in motype){
  
  inf<-modlist[grep(i,modlist$RNA),]
  cache[[i]]<-c(inf$Gene.Symbol)
  
}

motype<-paste0(unique(modlist$Subcellular),collapse = ',')
motype<-unlist(strsplit(motype,'[,]'))
motype<-unique(motype)

for (i in motype){
  
  inf<-modlist[grep(i,modlist$Subcellular),]
  cache[[i]]<-c(inf$Gene.Symbol)
  
}

cache <- Filter(function(x) length(x) > 4, cache)


###################################################
#把不同类型分析的结果按照disorder-based的方式统计出来
#把FDR的计算统一下再重新绘图

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
  
  gene_ranks <- (- log10(inf$P))
  names(gene_ranks) <- inf$V6
  
  gsea_results <- fgsea(
    pathways = cache, 
    stats = gene_ranks,
    minSize = 5,
    maxSize = 500,
    scoreType = "pos" 
  )
  
  gsea_results<-as.data.frame(gsea_results)
  gsea_results$disorder<-c(disease)
  
  filesave<-rbind(filesave,gsea_results)
}

filesave<-filesave[,-8]
filesave<-as.data.frame(filesave)
filesave$padj<-p.adjust(filesave$pval,method = 'BH')

write.csv(filesave,file = 'property/MAGMA_GSEA.csv')


#hMAGMA

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
  
  inf<-inf[!duplicated(inf$hgnc_symbol),]
  gene_ranks <- (- log10(inf$P))
  names(gene_ranks) <- inf$hgnc_symbol
  
  gsea_results <- fgsea(
    pathways = cache, 
    stats = gene_ranks,
    minSize = 5,
    maxSize = 500,
    scoreType = "pos" 
  )
  
  gsea_results<-as.data.frame(gsea_results)
  gsea_results$disorder<-c(disease)
  
  filesave<-rbind(filesave,gsea_results)
  
}

filesave<-filesave[,-8]
filesave<-as.data.frame(filesave)
filesave$padj<-p.adjust(filesave$pval,method = 'BH')

write.csv(filesave,file = 'property/HMAGMA_GSEA.csv')

#SMR_DLPFC

geneset<-read.csv('geneset2.csv')

file<-list.files('eQTL_BrainMeta/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

filesave<-data.frame()

for (i in list){
  
  combsmr<-data.frame()
  
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
  gene_ranks <- (- log10(inf$p_SMR)*sign(inf$b_SMR))
  names(gene_ranks) <- inf$Gene
  
  gsea_results <- fgsea(
    pathways = cache, 
    stats = gene_ranks,
    minSize = 5,
    maxSize = 500,
    scoreType = "std" 
  )
  
  gsea_results<-as.data.frame(gsea_results)
  gsea_results$disorder<-c(i)
  
  filesave<-rbind(filesave,gsea_results)
  
  
}

filesave<-filesave[,-8]
filesave<-as.data.frame(filesave)
filesave$padj<-p.adjust(filesave$pval,method = 'BH')

write.csv(filesave,file = 'property/SMR_DLPFC_GSEA.csv')

#FUSION_DLPFC

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
  
  combsmr<-combsmr[!is.na(combsmr$TWAS.P),]
  combsmr$fdr<-p.adjust(combsmr$TWAS.P,method = 'BH')
  
  inf<-combsmr[,c(3,7,11,14:21)]
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1

  inf<-inf[!duplicated(inf$ID),]
  rownames(inf)<-inf$order
  
  gene_ranks <- inf$TWAS.Z
  names(gene_ranks) <- inf$ID
  
  gsea_results <- fgsea(
    pathways = cache, 
    stats = gene_ranks,
    minSize = 5,
    maxSize = 500,
    scoreType = "std" 
  )
  
  gsea_results<-as.data.frame(gsea_results)
  gsea_results$disorder<-c(i)
  
  filesave<-rbind(filesave,gsea_results)
}

filesave<-filesave[,-8]
filesave<-as.data.frame(filesave)
filesave$padj<-p.adjust(filesave$pval,method = 'BH')

write.csv(filesave,file = 'property/FUSION_DLPFC_GSEA.csv')


#FUSION_blood

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
  
  combsmr<-combsmr[!is.na(combsmr$TWAS.P),]
  combsmr$fdr<-p.adjust(combsmr$TWAS.P,method = 'BH')
  
  inf<-combsmr[,c(3,7,11,14:21)]
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  inf<-inf[!duplicated(inf$ID),]
  rownames(inf)<-inf$order
  
  gene_ranks <- inf$TWAS.Z
  names(gene_ranks) <- inf$ID
  
  gsea_results <- fgsea(
    pathways = cache, 
    stats = gene_ranks,
    minSize = 5,
    maxSize = 500,
    scoreType = "std" 
  )
  
  gsea_results<-as.data.frame(gsea_results)
  gsea_results$disorder<-c(i)
  
  filesave<-rbind(filesave,gsea_results)
}

filesave<-filesave[,-8]
filesave<-as.data.frame(filesave)
filesave$padj<-p.adjust(filesave$pval,method = 'BH')

write.csv(filesave,file = 'property/FUSION_blood_GSEA.csv')


#########
#把magma原生的geneset放进去

list<-list.files('for_list2/list2_magma/')
list<-list[grep('geneset.gsa.out',list)]

rmrgenset<-data.frame()

for (i in list){
  inf<-read.table(paste0('for_list2/list2_magma/',i),header = T) 
  
  disease<-unlist(strsplit(i,'_'))[2]

  inf$or_uci95<-exp(inf$BETA + 1.96 * inf$SE)
  inf$or_lci95<-exp(inf$BETA - 1.96 * inf$SE)
  inf$or<-exp(inf$BETA)
  
  inf<-inf[inf$NGENES>4,]
  inf$disorder<-disease
  rmrgenset<-rbind(rmrgenset,inf)
}

write.csv(rmrgenset,file = 'property/MAGMA_COMP.csv')

#hMAGMA


list<-list.files('for_list2/list2_magma/')
list<-list[grep('geneset_hmagma.gsa.out',list)]
rmrgenset<-data.frame()

for (i in list){
  inf<-read.table(paste0('for_list2/list2_magma/',i),header = T) 
  
  disease<-unlist(strsplit(i,'_'))[2]
  #森林
  inf$or_uci95<-exp(inf$BETA + 1.96 * inf$SE)
  inf$or_lci95<-exp(inf$BETA - 1.96 * inf$SE)
  inf$or<-exp(inf$BETA)
  
  inf<-inf[inf$NGENES>4,]

  inf$disorder<-disease

  rmrgenset<-rbind(rmrgenset,inf)

}

write.csv(rmrgenset,file = 'property/hMAGMA_COMP.csv')


###

#可以先写一个demo，后续再把MAGMA的geneset什么的加进来
mgs<-read.csv('property/MAGMA_GSEA.csv')
mgs$ana<-c('MAGMA')
mgs2<-read.csv('property/HMAGMA_GSEA.csv')
mgs2$ana<-c('hMAGMA')
mgs3<-read.csv('property/SMR_DLPFC_GSEA.csv')
mgs3$ana<-c('SMR_DLPFC')
mgs4<-read.csv('property/FUSION_DLPFC_GSEA.csv')
mgs4$ana<-c('FUSION_DLPFC')
mgs5<-read.csv('property/FUSION_blood_GSEA.csv')
mgs5$ana<-c('FUSION_blood')

mgs6<-read.csv('property/MAGMA_COMP.csv')
mgs6$ana<-c('MAGMA_COMP')
mgs7<-read.csv('property/hMAGMA_COMP.csv')
mgs7$ana<-c('HMAGMA_COMP')

combi<-rbind(mgs,mgs2,mgs3,mgs4,mgs5)
combi$disorder[combi$disorder=='PDD']<-c('MDD')
plot1<-combi[combi$pathway=='RMP',]
plot1$padj<-p.adjust(plot1$pval,method = 'BH')

combi2<-rbind(mgs6,mgs7)
combi2$disorder[combi2$disorder=='PDD']<-c('MDD')
plot1s<-combi2[combi2$VARIABLE=='RMP',]
plot1s$fdr<-p.adjust(plot1s$P,method = 'BH')
##

library(dplyr)
library(tidyr)
library(tibble)

##left

wide_df <- plot1s %>%
  select(ana, disorder, BETA) %>%
  pivot_wider(names_from = ana, values_from = BETA)


heatmap_matrix <- wide_df %>%
  column_to_rownames("disorder") %>% 
  as.matrix()       

wide_df2 <- plot1s %>%
  select(ana, disorder, fdr) %>%
  pivot_wider(names_from = ana, values_from = fdr)


fdr <- wide_df2 %>%
  column_to_rownames("disorder") %>% 
  as.matrix()    
fdr[is.na(fdr)]<-1

library(circlize)
library(ComplexHeatmap)

col_fun1 = colorRamp2(c(-0.5,0,0.5), c('#22B5AF','white','#F57F17'))

pdf('./property/GSEA_left.pdf',width = 8,height = 8)

Heatmap(
  heatmap_matrix,col = col_fun1,
  name = "Beta Value",na_col = "gray90",
  cluster_columns = FALSE,  
  cluster_rows = F,    
  show_column_names = T,
  width = ncol(heatmap_matrix)*unit(6, "mm"),
  height = nrow(heatmap_matrix)*unit(4, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    current_fdr <- fdr[i, j]
    if (!is.na(current_fdr)) {
      if  (current_fdr < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "plain",alpha=0.8))
      } }})

dev.off()

##right

wide_df <- plot1 %>%
  select(ana, disorder, NES) %>%
  pivot_wider(names_from = ana, values_from = NES)


heatmap_matrix <- wide_df %>%
  column_to_rownames("disorder") %>% 
  as.matrix()       

wide_df2 <- plot1 %>%
  select(ana, disorder, padj) %>%
  pivot_wider(names_from = ana, values_from = padj)


fdr <- wide_df2 %>%
  column_to_rownames("disorder") %>% 
  as.matrix()    
fdr[is.na(fdr)]<-1

col_fun1 = colorRamp2(c(-2,0,2), c('#9d9dc7','white','#e3aba7'))

pdf('./property/GSEA_right.pdf',width = 8,height = 8)

Heatmap(
  heatmap_matrix,col = col_fun1,
  name = "NES score",na_col = "gray90",
  cluster_columns = FALSE,  
  cluster_rows = F,    
  show_column_names = T,
  width = ncol(heatmap_matrix)*unit(6, "mm"),
  height = nrow(heatmap_matrix)*unit(4, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    current_fdr <- fdr[i, j]
    if (!is.na(current_fdr)) {
      if  (current_fdr < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 12, fontface = "plain",alpha=0.8))
      } }})

dev.off()

#第一张热图的FDR是没问题的

###

plot2<-combi[combi$pathway %in% c('mRNA','tRNA','rRNA','snRNA','miRNA'),]
plot2s<-combi2[combi2$VARIABLE %in% c('mRNA','tRNA','rRNA','snRNA','miRNA'),]

plot2<-data.frame(disorder=c(plot2$disorder,plot2s$disorder),
                  pvalue=c(plot2$pval,plot2s$P),
                  pathway=c(plot2$pathway,plot2s$VARIABLE),
                  ana=c(plot2$ana,plot2s$ana))
plot2$padj<-p.adjust(plot2$pvalue,method = 'BH')

plot2$disorder<-factor(plot2$disorder,levels = rev(rownames(heatmap_mat)))

library(ggplot2)
library(MetBrewer)

global_shape_mapping <- c(
  'MAGMA_COMP'=10,
  'HMAGMA_COMP'=13,
  "MAGMA" = 16, 
  "hMAGMA" = 17, 
  "SMR_DLPFC" = 7,  
  "FUSION_DLPFC" = 15, 
  "FUSION_blood" = 3   )

p<-list()

p[[1]]<-
  ggplot(plot2, aes(y = disorder, x = -log10(padj), color  = pathway)) +
  geom_jitter(aes(shape = ana), height = 0.2, size = 2) +
  labs(y = NULL, x = "-log10(FDR)") +
  theme_bw() +
  scale_color_manual(values=c(met.brewer("Egypt", 5)))+
  scale_shape_manual(values = global_shape_mapping, guide = "none") +
  guides(color = guide_legend(title = NULL, ncol = 3)) +
  xlim(0, 1.2) +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())+
  theme(aspect.ratio=2)


plot3<-combi[combi$pathway %in% c('cytosol','mitochondrion','nucleus'),]
plot3s<-combi2[combi2$VARIABLE %in% c('cytosol','mitochondrion','nucleus'),]

plot3<-data.frame(disorder=c(plot3$disorder,plot3s$disorder),
                  pvalue=c(plot3$pval,plot3s$P),
                  pathway=c(plot3$pathway,plot3s$VARIABLE),
                  ana=c(plot3$ana,plot3s$ana))
plot3$padj<-p.adjust(plot3$pvalue,method = 'BH')
plot3$disorder<-factor(plot3$disorder,levels = rev(rownames(heatmap_mat)))

p[[2]]<-
  ggplot(plot3, aes(y = disorder, x = -log10(padj), color  = pathway)) +
  geom_jitter(aes(shape = ana), height = 0.2, size = 2) +
  labs(y = NULL, x = "-log10(FDR)") +
  theme_bw() +
  scale_color_manual(values=c(met.brewer("Egypt", 3)))+
  scale_shape_manual(values = global_shape_mapping, guide = "none") +
  guides(color = guide_legend(title = NULL, ncol = 3)) +
  xlim(0, 1.2) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())+
  theme(aspect.ratio=2)


plot4<-combi[!(combi$pathway %in% c('RMP','cytosol','mitochondrion','nucleus',
                                    'mRNA','tRNA','rRNA','snRNA','miRNA')),]
plot4s<-combi2[!(combi2$VARIABLE %in% c('RMP','cytosol','mitochondrion','nucleus',
                                        'mRNA','tRNA','rRNA','snRNA','miRNA')),]

plot4<-data.frame(disorder=c(plot4$disorder,plot4s$disorder),
                  pvalue=c(plot4$pval,plot4s$P),
                  pathway=c(plot4$pathway,plot4s$VARIABLE),
                  ana=c(plot4$ana,plot4s$ana))
plot4<-plot4[!(plot4$pathway %in% c('mcm5U','ncm5U','cm5U')),]
plot4$pathway[plot4$pathway=='mcm5s2U']<-c('U34')
plot4$padj<-p.adjust(plot4$pvalue,method = 'BH')
plot4$disorder<-factor(plot4$disorder,levels = rev(rownames(heatmap_mat)))


p[[3]]<-
  ggplot(plot4, aes(y = disorder, x = -log10(padj), color  = pathway)) +
  geom_jitter(aes(shape = ana), height = 0.2, size = 2) +
  labs(y = NULL, x = "-log10(FDR)") +
  theme_bw() +
  scale_color_manual(values=c(met.brewer("Egypt", 11)))+
  scale_shape_manual(values = global_shape_mapping, guide = "none") +
  guides(color = guide_legend(title = NULL, ncol = 3)) +
  xlim(0, 1.2) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank())+
  theme(aspect.ratio=2)

pdf("./property/GSEA_sub.pdf", width = 10, height = 7)
patchwork::wrap_plots(p, nrow = 1)
dev.off()


###

source("combine_gsea_tables.R")


#########################
#


scoretable<-read.csv('for_list2/scoretable_insup.csv') 

scoretable<-scoretable[scoretable$total.score>1,]
colnames(scoretable)[4]<-c('PIP')
scoretable<-scoretable[,-17]

df_wide <- scoretable %>%
  pivot_wider(
    names_from = Disorder,              
    values_from = GWAS:RNA.seq.isoform.expression,
    names_glue = "{Disorder}_{.value}"
  )

df_wide[is.na(df_wide)]<-0


min_max_scale <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(0) # 防止分母为0
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

df_wide <- df_wide %>%
  select_if(~ !is.numeric(.) || sum(. != 0, na.rm = TRUE) > 0) %>%
  mutate(across(-RMP, min_max_scale))

pca_matrix <- df_wide %>% 
  tibble::column_to_rownames("RMP") %>% 
  as.matrix()

pca_result <- prcomp(pca_matrix, center = TRUE, scale. = FALSE)

library(ggplot2)

pca_df <- as.data.frame(pca_result$x) %>%
  tibble::rownames_to_column("RMP")

pc1_var <- round(summary(pca_result)$importance[2, 1] * 100, 1) 
pc2_var <- round(summary(pca_result)$importance[2, 2] * 100, 1)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "#4A90E2", size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = RMP), size = 3, max.overlaps = 15) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(aspect.ratio=1)

scoretable<-read.csv('for_list2/scoretable_insup.csv') 

#大概是这种感觉，但是用这张表做总感觉怪怪的，要不把原始的统计整理出来？
#或者我想看一下两个维度是不是分别和BIPI和SCZ的风险相关？


#或者和一些已知的因子做对照？这样的话或许可以做一些非经验的权重？
#但是这样做的话会多出很多乱七八糟的分析我不太想加
#但是不加入twosampleMR什么的乱七八糟的分析又好像是可行的

#AD的话APOE，BIN1，SCZ：DRD2，GRIN2A，BIPI：ANK3，CACNA1C，ADHD：SLC6A3，DRD4
#对照的话：TYR，FABP4，ALB之类的在脑中不怎么影响的因子


#我想的是这里的PCA可以把figrue2的几个分析加上celltype，brianregion，一起放进去然后分析？

scoretable<-read.csv('for_list2/scoretable_insup.csv') 
pcagene<-unique(scoretable$RMP)
pcagene<-c(pcagene,'APOE','APP','PSEN1','PSEN2','BIN1','TREM2',
           'DRD2','GRIN2A','HTR2A','SYN2',
           'CACNA1C','ANK3','TRANK1','SYNE1','TENM4',
           'SLC6A3','DRD4','FOXP2','ST3GAL3',
           'TYR','FABP4','ALB','TNNT3','OR5M3')

#MAGMA

list<-list.files('MAGMA_singlegene/')
anno<-read.table('NCBI37.gene.loc',header = F)
anno$V6[anno$V6=='FTSJ2']<-'MRM2'

filesave<-data.frame()

for (i in list){
  
  inf<-read.table(paste0('MAGMA_singlegene/',i),header = T) 
  
  disease<-unlist(strsplit(i,'_'))[2]
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-merge(inf,anno,by.x=c('GENE'),by.y=c('V1'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$V6 %in% hla$V6),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[match(pcagene,inf$V6),]
  inf<-data.frame(Gene=c(pcagene),
                  Score=c(-log10(inf$fdr)), #inf$ZSTAT
                  Disorder=c(disease),
                  Method=c('MAGMA'))
  
  filesave<-rbind(filesave,inf)
  
}

#HMAGMA

list<-list.files('hMAGMA/hMAGMA_out/')
list<-list[grep('genes.out',list)]
load('hMAGMA/geneAnno_allgenes.rda')
geneAnno1$hgnc_symbol[geneAnno1$hgnc_symbol=='FTSJ2']<-c('MRM2')

for (i in list){
  
  inf<-read.table(paste0('hMAGMA/hMAGMA_out/',i),header = T)
  disease<-unlist(strsplit(i,'_'))[2]
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  
  inf<-merge(inf,geneAnno1,by.x=c('GENE'),by.y=c('ensembl_gene_id'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$hgnc_symbol %in% hla$hgnc_symbol),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[match(pcagene,inf$hgnc_symbol),]
  inf<-data.frame(Gene=c(pcagene),
                  Score=c(-log10(inf$fdr)), #inf$ZSTAT
                  Disorder=c(disease),
                  Method=c('HMAGMA'))
  
  filesave<-rbind(filesave,inf)
  
}

#SMR

file<-list.files('eQTL_BrainMeta/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

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
  
  inf<-combsmr[match(pcagene,combsmr$Gene),]
  
  
  inf<-data.frame(Gene=c(pcagene),
                  Score=c(log10(inf$fdr)*sign(inf$b_SMR)*-1),
                  Disorder=c(disease),
                  Method=c('SMR'))

  filesave<-rbind(filesave,inf)
  
}


#

file<-list.files('TWAS_result/TWAS_result_GTEx_cortex/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)


for (i in list){
  combsmr<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('TWAS_result/TWAS_result_GTEx_cortex/PGC_',i,'_for_fusion_TWAS_',j,'.dat'),header = T)
    combsmr<-rbind(combsmr,inf)
  }
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  combsmr<-combsmr[!is.na(combsmr$TWAS.P),]
  inf<-combsmr[match(pcagene,combsmr$ID),]
  inf<-data.frame(Gene=c(pcagene),
                  Score=c(inf$TWAS.Z),
                  Disorder=c(disease),
                  Method=c('FUSION'))
  
  filesave<-rbind(filesave,inf)
  
}


#

file<-list.files('TWAS_result/gtex_blood/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

for (i in list){
  combsmr<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('TWAS_result/gtex_blood/PGC_',i,'_for_fusion_TWAS_',j,'.dat'),header = T)
    combsmr<-rbind(combsmr,inf)
  }
  
  disease<-i
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  combsmr<-combsmr[!is.na(combsmr$TWAS.P),]
  inf<-combsmr[match(pcagene,combsmr$ID),]
  inf<-data.frame(Gene=c(pcagene),
                  Score=c(inf$TWAS.Z),
                  Disorder=c(disease),
                  Method=c('FUSION_blood'))
  
  filesave<-rbind(filesave,inf)
  
}

#然后就是celltype和brain region了

#SMR subregion

list<-list.files('brainregion/')
list<-list[grep('smr',list)]

library(stringr)
library(tidyr)

task1 <- unique(str_split_fixed(list, "_", 3)[, 2])
task2 <- unique(str_extract(list, "(?<=Brain_).*(?=\\.smr)"))

for (i in task1){
  
  multireg<-data.frame()
  
  for (j in task2){
    
    subsmr<-read.table(paste0('brainregion/PGC_',i,'_Brain_',j,'.smr'),header = T)
    subsmr<-subsmr[!is.na(subsmr$p_HEIDI),]
    subsmr<-subsmr[subsmr$p_HEIDI>0.05,]
    subsmr<-subsmr[match(pcagene,subsmr$Gene),]
    subsmr$Gene<-pcagene
    
    subsmr$region<-j
    multireg<-rbind(multireg,subsmr)
  }
  
  disease<-i
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-data.frame(Gene=c(multireg$Gene),
                  Score=c(c(log10(multireg$p_SMR)*sign(multireg$b_SMR)*-1)),
                  Disorder=c(disease),
                  Method=c(paste0('SMR_',multireg$region)))
  
  filesave<-rbind(filesave,inf)
}

#SMR cell type

list<-list.files('celltype/')
list<-list[grep('smr',list)]

task1 <- unique(str_split_fixed(list, "_", 3)[, 2])
task2 <- str_split(list, "_")
task2 <- sapply(task2, function(x) {
  combined <- paste(x[3:length(x)], collapse = "_")
  str_remove(combined, "\\.smr$")
})
task2 <- unique(task2)

for (i in task1){
  
  multireg<-data.frame()
  
  for (j in task2){
    
    subsmr<-read.table(paste0('celltype/PGC_',i,'_',j,'.smr'),header = T)
    subsmr<-subsmr[!is.na(subsmr$p_HEIDI),]
    subsmr<-subsmr[subsmr$p_HEIDI>0.05,]
    subsmr<-subsmr[match(pcagene,subsmr$Gene),]
    subsmr$Gene<-pcagene
    
    subsmr$region<-j
    multireg<-rbind(multireg,subsmr)
  }
  
  disease<-i
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-data.frame(Gene=c(multireg$Gene),
                  Score=c(c(log10(multireg$p_SMR)*sign(multireg$b_SMR)*-1)),
                  Disorder=c(disease),
                  Method=c(paste0('SMR_',multireg$region)))
  
  filesave<-rbind(filesave,inf)
  
}


supinf<-read.csv('property/snpforfigure3.csv') #基于原始GWAS信号的统计
supinf<-supinf[,-1]

filesave<-rbind(filesave,supinf)

table(filesave$Method)


################################################

#折腾了一天了这个结果，但依然不怎么理想我真是

filesave_wide <- filesave[filesave$Disorder %in% c('SCZ'),] %>%
  pivot_wider(
    id_cols = Gene,
    names_from = c(Disorder, Method),
    names_sep = "_",
    values_from = Score)

# filesave_wide<-filesave_wide[,-5]
print(filesave_wide[1:2, ])

filesave_wide <- filesave_wide %>%
  select_if(~ !is.numeric(.) || (sum(!is.na(.)) > 13 && sum(. != 0, na.rm = TRUE) > 0))

library(pcaMethods)
pca_matrix<-as.data.frame(filesave_wide)
rownames(pca_matrix)<-pca_matrix$Gene
pca_matrix<-pca_matrix[,-1]

sczh2<-read.csv('sczh2_result.csv')
sczh2<-sczh2[match(rownames(pca_matrix),sczh2$GENE),]
sczh2$ratio<-sczh2$H2/sczh2$length
pca_matrix$sch2<-scale(sczh2$pve)
pca_matrix$sch2ratio<-scale(sczh2$ratio)
# ppca_result <- pcaMethods::pca(pca_matrix, method = "ppca", nPcs = 5)
# pca_df <- as.data.frame(scores(ppca_result))
# 
# library(ggplot2)
# 
# pca_df <- as.data.frame(pca_df) %>%
#   tibble::rownames_to_column("RMP")

#AD的话APOE，BIN1，SCZ：DRD2，GRIN2A，BIPI：ANK3，CACNA1C，ADHD：SLC6A3，DRD4

# 'APOE','APP','PSEN1','PSEN2','BIN1','TREM2',
# 'DRD2','GRIN2A',
# 'CACNA1C','ANK3','TRANK1',
# 'SLC6A3','DRD4','FOXP2','ST3GAL3',
# 对照的话: 'TYR','FABP4','ALB','TNNT3','OR5M3'

#for SCZ

gene_annotation <- data.frame(
  Gene = c(pcagene),
  Group = c(rep("Unknown", 42), rep("Other", 6),rep("Positive", 2),
            rep("Other", 2),rep("Positive",1),
            rep("Other", 8),rep("Negative", 5)))

# pca_matrix[is.na(pca_matrix)]<-0
# pca_matrix<-abs(pca_matrix)
ppca_result <- pcaMethods::pca(pca_matrix, method = "ppca", nPcs = 2) #,scale='pareto',center = T

completed_data <- as.data.frame(ppca_result@completeObs) %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(gene_annotation, by = "Gene")

train_data <- completed_data %>% filter(Group %in% c("Positive", "Negative"))
test_data <- completed_data %>% filter(Group == "Unknown")

library(mixOmics)
library(dplyr)

X_train <- train_data %>% dplyr::select(-Gene, -Group) %>% as.matrix()
Y_train <- train_data$Group 

plsda_model <- splsda(X_train, Y_train, ncomp = 3)
predict(plsda_model, newdata = X_train)$predict[, "Positive", ncomp = 3]

X_test <- test_data %>% dplyr::select(-Gene, -Group) %>% as.matrix()
plsda_pred <- predict(plsda_model, newdata = X_test)

unknown_scores <- plsda_pred$predict[, "Positive", ncomp = 3]

plsda_ranking <- data.frame(
  Gene = test_data$Gene,
  Pos_Score = unknown_scores
) %>%
  arrange(desc(Pos_Score)) 

print(plsda_ranking)

#


#AD的话APOE，BIN1，SCZ：DRD2，GRIN2A，BIPI：ANK3，CACNA1C，ADHD：SLC6A3，DRD4

# 'APOE','APP','PSEN1','PSEN2','BIN1','TREM2',
# 'DRD2','GRIN2A',‘DRD3’
# 'CACNA1C','ANK3','TRANK1',
# 'SLC6A3','DRD4','FOXP2','ST3GAL3',
# 对照的话: 'TYR','FABP4','ALB','TNNT3','OR5M3'

filesave_wide <- filesave[filesave$Disorder %in% c('BIPI'),] %>%
  pivot_wider(
    id_cols = Gene,
    names_from = c(Disorder, Method),
    names_sep = "_",
    values_from = Score)

# filesave_wide<-filesave_wide[,-5]
print(filesave_wide[1:2, ])

filesave_wide <- filesave_wide %>%
  select_if(~ !is.numeric(.) || (sum(!is.na(.)) > 8 && sum(. != 0, na.rm = TRUE) > 0))

pca_matrix<-as.data.frame(filesave_wide)
rownames(pca_matrix)<-pca_matrix$Gene
pca_matrix<-pca_matrix[,-1]

sczh2<-read.csv('BIPIh2_result.csv')
sczh2<-sczh2[match(rownames(pca_matrix),sczh2$GENE),]
sczh2$ratio<-(sczh2$H2/sczh2$length)
pca_matrix$sch2<-sczh2$pve
pca_matrix$sch2ratio<-sczh2$ratio

#for BIPI

gene_annotation <- data.frame(
  Gene = c(pcagene),
  Group = c(rep("Unknown", 42), rep("Other", 10),
            rep("Positive", 5), rep("Other", 4),rep("Negative", 5))
)

# pca_matrix[is.na(pca_matrix)]<-0
# pca_matrix<-abs(pca_matrix)
ppca_result <- pcaMethods::pca(pca_matrix, method = "ppca", nPcs = 2) #,scale = "pareto", center = T

completed_data <- as.data.frame(ppca_result@completeObs) %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(gene_annotation, by = "Gene")

train_data <- completed_data %>% filter(Group %in% c("Positive", "Negative"))
test_data <- completed_data %>% filter(Group == "Unknown")

X_train <- train_data %>% dplyr::select(-Gene, -Group) %>% as.matrix()
Y_train <- train_data$Group 

plsda_model <- splsda(X_train, Y_train, ncomp = 3)
predict(plsda_model, newdata = X_train)$predict[, "Positive", ncomp = 3]

X_test <- test_data %>% dplyr::select(-Gene, -Group) %>% as.matrix()
plsda_pred <- predict(plsda_model, newdata = X_test)

unknown_scores <- plsda_pred$predict[, "Positive", ncomp = 3]


plsda_ranking2 <- data.frame(
  Gene = test_data$Gene,
  Pos_Score = unknown_scores
) %>%
  arrange(desc(Pos_Score)) 

print(plsda_ranking2)

plot<-merge(plsda_ranking,plsda_ranking2,by = c('Gene'))
colnames(plot)<-c('Gene','SCZscore','BIPIscore')

####


ggplot(plot[plot$Gene %in% scoretable$RMP[scoretable$total.score>1] ,], aes(x = SCZscore, y = BIPIscore)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = Gene), size = 3) +
  theme_bw() +
  labs(title = "")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)+
  theme(aspect.ratio=1)




