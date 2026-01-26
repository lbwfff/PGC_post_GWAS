
geneset<-read.csv('geneset2.csv')

file<-list.files('eQTL_BrainMeta/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

p<-list()

for (i in list){
  cache<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('eQTL_BrainMeta/PGC_',i,'_BrainMeta_cis_eQTL_chr',j,'.smr'),header = T)
    cache<-rbind(cache,inf)
  }
  
  cache<-cache[cache$p_HEIDI>0.05  & !is.na(cache$p_HEIDI),]
  cache$fdr<-p.adjust(cache$p_SMR,method = 'BH')
  
  inf<-cache
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  
  p[[which(list==i)]]<-
    ggplot(inf[inf$order<5000,], aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$Gene %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
               color = '#0f7bbb', size = 3) +  
    geom_text_repel(data = inf[inf$Gene %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
                    aes(label = Gene),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE, 
                    nudge_x = 250,  nudge_y = 1, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('SMR analyses for ',i),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "top",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    theme(aspect.ratio=1)
  
}

pdf("for_list2/SMR_eQTL.pdf",width = 16,height = 16)
wrap_plots(p,nrow=3) 
dev.off()



#######


file<-list.files('sQTL/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

p<-list()

for (i in list){
  cache<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('sQTL/PGC_',i,'_BrainMeta_cis_sQTL_chr',j,'.smr'),header = T)
    cache<-rbind(cache,inf)
  }
  
  cache<-cache[cache$p_HEIDI>0.05  & !is.na(cache$p_HEIDI),]
  cache$fdr<-p.adjust(cache$p_SMR,method = 'BH')
  
  inf<-cache
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  
  p[[which(list==i)]]<-
    ggplot(inf[inf$order<2000,], aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$Gene %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
               color = '#0f7bbb', size = 3) +  
    geom_text_repel(data = inf[inf$Gene %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
                    aes(label = Gene),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE,  max.overlaps = Inf,
                    nudge_x = 200,  nudge_y =0.5, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('SMR analyses for ',i),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "top",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    theme(aspect.ratio=1)
  
}

library(patchwork)

pdf("for_list2/SMR_sQTL.pdf",width = 16,height = 16)
wrap_plots(p,nrow=3) 
dev.off()


#####NSUN2的另外三个 splicing
# chr5:6605521:6606933, 0.1168480
# chr5:6618137:6620219, 0.1500020



#chr5:6620411:6623327 -0.11

#########
#NSUN2的两个eQTL和9个sQTL事件都画一下effect size plot

source('SMRlocus.R')

p<-list()

readsmr<-ReadSMRData('./SMR_plot/NSUN2_BIPI.ENSG00000037474.15.txt')
p[[1]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_SCZ.ENSG00000037474.15.txt')
p[[2]]<-SMREffectPlot_v2(readsmr)

pdf("SMR_effect_plot_NSUN2_eQTL.pdf",width = 8,height = 4)
patchwork::wrap_plots(p,nrow=1) 
dev.off()

p<-list()

readsmr<-ReadSMRData('./SMR_plot/NSUN2_sQTL_BIP.chr5_6632869_6633744_clu_41364_.txt')
p[[1]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_sQTL3_BIP.chr5_6621857_6623327_clu_41362_.txt')
p[[2]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_sQTL2_BIP.chr5_6605488_6606933_clu_41360_.txt')
p[[3]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_sQTL2_SCZ.chr5_6632869_6633744_clu_41364_.txt')
p[[4]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_sQTL_SCZ.chr5_6621857_6623327_clu_41362_.txt')
p[[5]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./SMR_plot/NSUN2_sQTL3_SCZ.chr5_6605488_6606933_clu_41360_.txt')
p[[6]]<-SMREffectPlot_v2(readsmr)

pdf("SMR_effect_plot_NSUN2_sQTL.pdf",width = 12,height = 8)
patchwork::wrap_plots(p,nrow=2) 
dev.off()

p<-list()

readsmr<-ReadSMRData('./NSUN2_sQTL_BIPI_add2.chr5_6620411_6623327_clu_41362_.txt')
p[[1]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./NSUN2_sQTL_BIPI_add3.chr5_6618137_6620219_clu_41361_.txt')
p[[2]]<-SMREffectPlot_v2(readsmr)

readsmr<-ReadSMRData('./NSUN2_sQTL_BIPI_add1.chr5_6605521_6606933_clu_41360_.txt')
p[[3]]<-SMREffectPlot_v2(readsmr)

pdf("SMR_effect_plot_NSUN2_sQTL_sup.pdf",width = 12,height = 4)
patchwork::wrap_plots(p,nrow=1) 
dev.off()
