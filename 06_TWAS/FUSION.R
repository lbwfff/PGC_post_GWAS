geneset<-read.csv('geneset2.csv')

file<-list.files('TWAS_result/TWAS_result_GTEx_cortex/')
list<-sapply(file, function(x) unlist(strsplit(x, '_'))[2])
list<-unique(list)

p<-list()

for (i in list){
  cache<-data.frame()
  
  for (j in 1:22){
    inf<-read.table(paste0('TWAS_result/TWAS_result_GTEx_cortex/PGC_',i,'_for_fusion_TWAS_',j,'.dat'),header = T)
    cache<-rbind(cache,inf)
  }
  
  cache<-cache[!is.na(cache$TWAS.P),]
  cache$fdr<-p.adjust(cache$TWAS.P,method = 'BH')
  
  inf<-cache[,c(3,7,11,14:21)]
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  inf<-inf[!duplicated(inf$ID),] #为什么会有重复的基因
  rownames(inf)<-inf$order
  
  p[[which(list==i)]]<-
    ggplot(inf, aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$ID %in% geneset$Gene.Symbol & inf$group == 'sign', ],
               color = '#0f7bbb', size = 3) +
    geom_text_repel(data = inf[inf$ID %in% geneset$Gene.Symbol & inf$group == 'sign', ],
                    aes(label = ID),color = '#0f7bbb',size = 3,segment.color = "black",
                    show.legend = FALSE,
                    nudge_x = 80,  nudge_y = 0.2,
                    box.padding = 0.5, point.padding = 0.5) +
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('TWAS analyses for ',i),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "top",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    theme(aspect.ratio=1)
  
}

library('patchwork')

pdf("for_list2/FUSION.pdf",width = 16,height = 16)
wrap_plots(p,nrow=4) 
dev.off()



