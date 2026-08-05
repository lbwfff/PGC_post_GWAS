library(tidyverse)
library(ggtext)
library(normentR)

################MAGMA#########################

modlist<-read.csv('geneset2.csv')
modlist<-modlist[modlist$Gene.Symbol !='',]

gpos<-anno[anno$V6 %in% modlist$Gene.Symbol,]

library(table1)
table1(~ Effector | Modification, data=modlist,
       overall=c(left="Total"))

table(modlist$Modification,modlist$RNA)
unique(modlist$Modification)
##

library(cowplot)

geneset<-read.csv('geneset2.csv')
geneset<-geneset[geneset$Gene.Symbol !='',]

list<-list.files('MAGMA_singlegene/')
anno<-read.table('NCBI37.gene.loc',header = F)

genelist<-list(RNAME=as.character(geneset$Gene.Symbol))

p<-list()
p2<-list()

for (i in list){
  
  inf<-read.table(paste0('MAGMA_singlegene/',i),header = T) #这个数据形式.....是不是可以做GSEA？
  
  disease<-unlist(strsplit(i,'_'))[2]
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-merge(inf,anno,by.x=c('GENE'),by.y=c('V1'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$V6 %in% hla$V6),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  # gsea<-inf$ZSTAT
  # names(gsea)<-inf$V6
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  
  p[[which(list==i)]]<-
    ggplot(inf[inf$order<5000,], aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$V6 %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
               color = '#0f7bbb', size = 3) +  
    geom_text_repel(data = inf[inf$V6 %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
                    aes(label = V6),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE, 
                    nudge_x = 350,  nudge_y = 0.85, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('Gene-based analyses for ',disease),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "None",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    theme(aspect.ratio=1)
  
  #QQ plot
  inf<-inf[order(inf$P),]
  ci <- 0.95
  nsnp <- nrow(inf)
  
  plotdata1 <- tibble(
    observed = -log10(sort(inf$P)),
    expected = -log10(ppoints(nsnp)),
    clower = -log10(qbeta(
      p = (1 - ci) / 2,
      shape1 = seq(nsnp),
      shape2 = rev(seq(nsnp)))),
    cupper = -log10(qbeta(
      p = (1 + ci) / 2,
      shape1 = seq(nsnp),
      shape2 = rev(seq(nsnp)))))
  plotdata1$gene<-inf$V6
  
  p2[[which(list==i)]]<-
    ggplot(plotdata1,aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cupper, ymin = clower),
                fill = "grey30", alpha = 0.5) +
    geom_step(data = plotdata1, 
              linewidth = 0.5, direction = "vh") +
    geom_segment(data = plotdata1 %>% filter(expected == max(expected)),
                 aes(x = 0, xend = expected, y = 0, yend = expected),
                 linewidth = 1.25, alpha = 0.5,
                 color = "grey30", lineend = "round") +
    labs(
      x = str_glue("Expected -log10(P)"),
      y = str_glue("Observed -log10(P)"),
      title = paste0('Gene-based analyses for ',disease)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10) ,
          panel.grid = element_blank()) +
    geom_point(data = plotdata1[plotdata1$gene %in% inf$V6[inf$group=='sign'] &
                                  plotdata1$gene %in% geneset$Gene.Symbol,], 
               color = '#0f7bbb', size = 1.5) + 
    geom_text_repel(data = plotdata1[plotdata1$gene %in% inf$V6[inf$group=='sign'] &
                                       plotdata1$gene %in% geneset$Gene.Symbol,], 
                    aes(label = gene),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE,  nudge_y = 1, nudge_x = -1, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme(aspect.ratio=1)
  
}

library(patchwork)

pdf("./for_list2/MAGMA_genelevel_withouthla.pdf",width = 16,height = 16)
wrap_plots(p,nrow=4) 
dev.off()

pdf("./for_list2/MAGMA_genelevel_qq_withouthla.pdf",width = 16,height = 16)
wrap_plots(p2,nrow=4) 
dev.off()

#for SCZ

for (i in list[15]){
  
  inf<-read.table(paste0('MAGMA_singlegene/',i),header = T) 
  disease<-unlist(strsplit(i,'_'))[2]
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  inf<-merge(inf,anno,by.x=c('GENE'),by.y=c('V1'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$V6 %in% hla$V6),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  
  pdf("for_list2/SCZ_MAGMA_genelevel_2.pdf",width = 6,height = 6) 
  ggplot(inf, aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$V6 %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
               color = '#0f7bbb', size = 3) +  
    geom_text_repel(data = inf[inf$V6 %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
                    aes(label = V6),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE, max.overlaps = Inf,
                    nudge_x = -350,  nudge_y = 0.8, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('Gene-based analyses for ',disease),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "top",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    coord_cartesian(ylim=c(0,10),xlim=c(0,3000))+
    theme(aspect.ratio=1)
  dev.off()
  
}


##################################
######geneset

list<-list.files('for_list2/list2_magma/')
list<-list[grep('geneset.gsa.out',list)]

p<-list()
p2<-list()
p3<-list()
rmrgenset<-data.frame()

for (i in list){
  inf<-read.table(paste0('for_list2/list2_magma/',i),header = T) 
  
  disease<-unlist(strsplit(i,'_'))[2]
  #森林
  inf$or_uci95<-exp(inf$BETA + 1.96 * inf$SE)
  inf$or_lci95<-exp(inf$BETA - 1.96 * inf$SE)
  inf$or<-exp(inf$BETA)
  
  inf<-inf[inf$NGENES>4,]
  rmrp<-inf[inf$VARIABLE=='RMP',]
  rmrp$disease<-disease
  rmrgenset<-rbind(rmrgenset,rmrp)
  
  plot1<-inf [inf$VARIABLE %in% c('I','m1A','m3C','m6A','cm5U',
                                  '2′-O-Me','Cap','Ψ','m5C','m5U'),]
  plot1$VARIABLE<-sub('Ψ','Psi',plot1$VARIABLE)
  plot1$VARIABLE<-sub('cm5U','U34 modifications',plot1$VARIABLE)
  plot2<-inf [inf$VARIABLE %in% c('mRNA','tRNA','miRNA','rRNA','snRNA'),]
  plot3<-inf [inf$VARIABLE %in% c('cytosol','nucleus','mitochondrion'),]
  
  plot1$FDR<-p.adjust(plot1$P,method = 'BH')
  plot1$sign<-ifelse(plot1$FDR<0.05,'*','')
  plot1$VARIABLE<-paste(plot1$VARIABLE,plot1$sign)
  plot1$VARIABLE<-factor(plot1$VARIABLE,levels = c(plot1[order(plot1$or,decreasing = F),]$VARIABLE))
  
  p[[which(list==i)]]<-
    ggplot(plot1, aes(y=reorder(VARIABLE,or), x=or,  colour=VARIABLE)) + 
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", nrow(plot1),type ='continuous'))+
    geom_vline(xintercept=1, linetype=3) +
    theme_minimal_hgrid(10, rel_small = 1) +
    labs(color = "",y = "", x = "Odds Ratio",
         title= paste0("Gene-set level analysis for ",
                       disease ,", 95% CI") )+
    theme(legend.position = "none")+
    theme(aspect.ratio=1)
  
  plot2$FDR<-p.adjust(plot2$P,method = 'BH')
  plot2$sign<-ifelse(plot2$FDR<0.05,'*','')
  plot2$VARIABLE<-paste(plot2$VARIABLE,plot2$sign)
  plot2$VARIABLE<-factor(plot2$VARIABLE,levels = c(plot2[order(plot2$or,decreasing = F),]$VARIABLE))
  
  p2[[which(list==i)]]<-
    ggplot(plot2, aes(y=reorder(VARIABLE,or), x=or,  colour=VARIABLE)) + 
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", nrow(plot2),type ='continuous'))+
    geom_vline(xintercept=1, linetype=3) +
    theme_minimal_hgrid(10, rel_small = 1) +
    labs(color = "",y = "", x = "Odds Ratio",
         title= paste0("Gene-set level analysis for ",
                       disease ,", 95% CI") )+
    theme(legend.position = "none")+
    theme(aspect.ratio=1)
  
  plot3$FDR<-p.adjust(plot3$P,method = 'BH')
  plot3$sign<-ifelse(plot3$FDR<0.05,'*','')
  plot3$VARIABLE<-paste(plot3$VARIABLE,plot3$sign)
  plot3$VARIABLE<-factor(plot3$VARIABLE,levels = c(plot3[order(plot3$or,decreasing = F),]$VARIABLE))
  
  p3[[which(list==i)]]<-
    ggplot(plot3, aes(y=reorder(VARIABLE,or), x=or,  colour=VARIABLE)) + 
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", nrow(plot3),type ='continuous'))+
    geom_vline(xintercept=1, linetype=3) +
    theme_minimal_hgrid(10, rel_small = 1) +
    labs(color = "",y = "", x = "Odds Ratio",
         title= paste0("Gene-set level analysis for ",
                       disease ,", 95% CI") )+
    theme(legend.position = "none")+
    theme(aspect.ratio=1)
  
}


pdf("for_list2/MAGMA_genesetlevel_1.pdf",width = 16,height = 16)
wrap_plots(p,nrow=4) 
dev.off()

pdf("for_list2/MAGMA_genesetlevel_2.pdf",width = 12,height = 12)
wrap_plots(p2,nrow=4) 
dev.off()

pdf("for_list2/MAGMA_genesetlevel_3.pdf",width = 12,height = 12)
wrap_plots(p3,nrow=4) 
dev.off()

rmrgenset$fdr<-p.adjust(rmrgenset$P,method = 'BH')
rmrgenset$sign<-ifelse(rmrgenset$fdr<0.05,'*','')
rmrgenset$disease<-paste(rmrgenset$disease,rmrgenset$sign)

pdf("for_list2/MAGMA_RMRset.pdf",width = 6,height = 6)
ggplot(rmrgenset, aes(y=reorder(disease,or), x=or,  colour=disease)) + 
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 16,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds Ratio",
       title= paste0("Gene-set level analysis for RMR, 95% CI") )+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)
dev.off()

################hMAGMA#########################

library(cowplot)

geneset<-read.csv('geneset2.csv')
geneset<-geneset[geneset$Gene.Symbol !='',]

list<-list.files('hMAGMA/hMAGMA_out/')
list<-list[grep('genes.out',list)]
load('hMAGMA/geneAnno_allgenes.rda')

genelist<-list(RNAME=as.character(geneset$Gene.Symbol))

p<-list()
p2<-list()
# p3<-list()
# p4<-list()

for (i in list){
  
  inf<-read.table(paste0('hMAGMA/hMAGMA_out/',i),header = T) #这个数据形式.....是不是可以做GSEA？
  disease<-unlist(strsplit(i,'_'))[2]
  
  if(disease=='MDD') {disease<-'PDD'} else{disease<-disease}
  
  
  inf<-merge(inf,geneAnno1,by.x=c('GENE'),by.y=c('ensembl_gene_id'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$hgnc_symbol %in% hla$hgnc_symbol),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  
  p[[which(list==i)]]<-
    ggplot(inf[inf$order<5000,,], aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$hgnc_symbol %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
               color = '#0f7bbb', size = 3) +  
    geom_text_repel(data = inf[inf$hgnc_symbol %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
                    aes(label = hgnc_symbol),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE, 
                    nudge_x = 350,  nudge_y = 0.85, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('Gene-based analyses for ',disease),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "None",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    theme(aspect.ratio=1)
  
  #QQ plot
  inf<-inf[order(inf$P),]
  ci <- 0.95
  nsnp <- nrow(inf)
  
  plotdata1 <- tibble(
    observed = -log10(sort(inf$P)),
    expected = -log10(ppoints(nsnp)),
    clower = -log10(qbeta(
      p = (1 - ci) / 2,
      shape1 = seq(nsnp),
      shape2 = rev(seq(nsnp)))),
    cupper = -log10(qbeta(
      p = (1 + ci) / 2,
      shape1 = seq(nsnp),
      shape2 = rev(seq(nsnp)))))
  plotdata1$gene<-inf$hgnc_symbol
  
  p2[[which(list==i)]]<-
    ggplot(plotdata1,aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cupper, ymin = clower),
                fill = "grey30", alpha = 0.5) +
    geom_step(data = plotdata1, 
              linewidth = 0.5, direction = "vh") +
    geom_segment(data = plotdata1 %>% filter(expected == max(expected)),
                 aes(x = 0, xend = expected, y = 0, yend = expected),
                 linewidth = 1.25, alpha = 0.5,
                 color = "grey30", lineend = "round") +
    labs(
      x = str_glue("Expected -log10(P)"),
      y = str_glue("Observed -log10(P)"),
      title = paste0('Gene-based analyses for ',disease)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10) ,
          panel.grid = element_blank()) +
    geom_point(data = plotdata1[plotdata1$gene %in% inf$hgnc_symbol[inf$group=='sign'] &
                                  plotdata1$gene %in% geneset$Gene.Symbol,], 
               color = '#0f7bbb', size = 1.5) + 
    geom_text_repel(data = plotdata1[plotdata1$gene %in% inf$hgnc_symbol[inf$group=='sign'] &
                                       plotdata1$gene %in% geneset$Gene.Symbol,], 
                    aes(label = gene),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE,  nudge_y = 1, nudge_x = -1, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme(aspect.ratio=1)
  
  # inf<-inf[inf$V6 %in% geneset$Gene.Symbol,]
  
}

library(patchwork)

pdf("for_list2/hMAGMA_genelevel_withouthla.pdf",width = 16,height = 16)
wrap_plots(p,ncol=4) 
dev.off()

pdf("for_list2/hMAGMA_genelevel_qq_withouthla.pdf",width = 16,height = 16)
wrap_plots(p2,ncol=4) 
dev.off()

#for SCZ

for (i in list[15]){
  
  inf<-read.table(paste0('hMAGMA/hMAGMA_out/',i),header = T) 
  disease<-unlist(strsplit(i,'_'))[2]
  inf<-merge(inf,geneAnno1,by.x=c('GENE'),by.y=c('ensembl_gene_id'))
  
  hla<-inf[inf$CHR==6 & inf$START > 25000000 & inf$STOP<33000000,]
  inf<-inf[!(inf$hgnc_symbol %in% hla$hgnc_symbol),]
  
  inf$fdr<-p.adjust(inf$P,method = c('BH'))
  inf<-inf[order(inf$fdr,decreasing = T),]
  inf$order<-nrow(inf):1
  
  
  library('ggrepel')
  inf$group<-ifelse(inf$fdr<0.05,'sign','non-sign')
  
  pdf("SCZ_hMAGMA_genelevel_2.pdf",width = 6,height = 6) 
  ggplot(inf, aes(x = order, y = -log10(fdr), color = group)) +
    geom_point(alpha = 0.7, size = 3) +  
    scale_color_manual(values = c("#A0A0A0", "#96bfce")) +
    geom_point(data = inf[inf$hgnc_symbol %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
               color = '#0f7bbb', size = 3) +  
    geom_text_repel(data = inf[inf$hgnc_symbol %in% geneset$Gene.Symbol & inf$group == 'sign', ], 
                    aes(label = hgnc_symbol),color = '#0f7bbb',size = 3,segment.color = "black", 
                    show.legend = FALSE, max.overlaps = Inf,
                    nudge_x = -2600,  nudge_y = 1.2, 
                    box.padding = 0.5, point.padding = 0.5) + 
    theme_bw() +
    theme_classic(base_size = 15) + 
    labs(title = paste0('Gene-based analyses for ',disease),
         y = '-log10 (FDR)', x = NULL,color = 'Group') + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.position = "top",  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10)  ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5)+
    coord_cartesian(ylim=c(0,10),xlim=c(250,4250))+
    theme(aspect.ratio=1)
  dev.off()
  
  inf<-inf[order(inf$ZSTAT,decreasing = T),]
  gsea<-inf$ZSTAT
  names(gsea)<-inf$V6
  
  #QQ plot
  inf<-inf[order(inf$P),]
  ci <- 0.95
  nsnp <- nrow(inf)
  
  plotdata1 <- tibble(
    observed = -log10(sort(inf$P)),
    expected = -log10(ppoints(nsnp)),
    clower = -log10(qbeta(
      p = (1 - ci) / 2,
      shape1 = seq(nsnp),
      shape2 = rev(seq(nsnp)))),
    cupper = -log10(qbeta(
      p = (1 + ci) / 2,
      shape1 = seq(nsnp),
      shape2 = rev(seq(nsnp)))))
  plotdata1$gene<-inf$hgnc_symbol
  
  pdf("SCZ_hMAGMA_genelevel_qq_2.pdf",width = 6,height = 6) 
  ggplot(plotdata1,aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymax = cupper, ymin = clower),
                fill = "grey30", alpha = 0.5) +
    geom_step(data = plotdata1, 
              linewidth = 0.5, direction = "vh") +
    geom_segment(data = plotdata1 %>% filter(expected == max(expected)),
                 aes(x = 0, xend = expected, y = 0, yend = expected),
                 linewidth = 1.25, alpha = 0.5,
                 color = "grey30", lineend = "round") +
    labs(
      x = str_glue("Expected -log10(P)"),
      y = str_glue("Observed -log10(P)"),
      title = paste0('Gene-based analyses for ',disease)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
          legend.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 10) ,
          panel.grid = element_blank()) +
    geom_point(data = plotdata1[plotdata1$gene %in% inf$hgnc_symbol[inf$group=='sign'] &
                                  plotdata1$gene %in% geneset$Gene.Symbol,], 
               color = '#8282aa', size = 1.5) + 
    geom_text_repel(data = plotdata1[plotdata1$gene %in% inf$hgnc_symbol[inf$group=='sign'] &
                                       plotdata1$gene %in% geneset$Gene.Symbol,], 
                    aes(label = gene),color = '#8282aa',size = 3,segment.color = "black", 
                    show.legend = FALSE,  nudge_y = 1, nudge_x = -1, 
                    box.padding = 0.5, point.padding = 0.5) + 
    coord_cartesian(ylim=c(0,6),xlim=c(0,3))+
    theme(aspect.ratio=1)
  dev.off()
}


######geneset

list<-list.files('for_list2/list2_magma/')
list<-list[grep('geneset_hmagma.gsa.out',list)]

p<-list()
p2<-list()
p3<-list()
rmrgenset<-data.frame()

for (i in list){
  inf<-read.table(paste0('for_list2/list2_magma/',i),header = T) 
  
  disease<-unlist(strsplit(i,'_'))[2]
  #森林
  inf$or_uci95<-exp(inf$BETA + 1.96 * inf$SE)
  inf$or_lci95<-exp(inf$BETA - 1.96 * inf$SE)
  inf$or<-exp(inf$BETA)
  
  inf<-inf[inf$NGENES>4,]
  rmrp<-inf[inf$VARIABLE=='RMP',]
  rmrp$disease<-disease
  rmrgenset<-rbind(rmrgenset,rmrp)
  
  plot1<-inf [inf$VARIABLE %in% c('I','m1A','m3C','m6A','cm5U',
                                  '2′-O-Me','Cap','Ψ','m5C','m5U'),]
  plot1$VARIABLE<-sub('Ψ','Psi',plot1$VARIABLE)
  plot1$VARIABLE<-sub('cm5U','U34 modifications',plot1$VARIABLE)
  plot2<-inf [inf$VARIABLE %in% c('mRNA','tRNA','miRNA','rRNA','snRNA'),]
  plot3<-inf [inf$VARIABLE %in% c('cytosol','nucleus','mitochondrion'),]
  
  plot1$FDR<-p.adjust(plot1$P,method = 'BH')
  plot1$sign<-ifelse(plot1$FDR<0.05,'*','')
  plot1$VARIABLE<-paste(plot1$VARIABLE,plot1$sign)
  plot1$VARIABLE<-factor(plot1$VARIABLE,levels = c(plot1[order(plot1$or,decreasing = F),]$VARIABLE))
  
  p[[which(list==i)]]<-
    ggplot(plot1, aes(y=reorder(VARIABLE,or), x=or,  colour=VARIABLE)) + 
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", nrow(plot1),type ='continuous'))+
    geom_vline(xintercept=1, linetype=3) +
    theme_minimal_hgrid(10, rel_small = 1) +
    labs(color = "",y = "", x = "Odds Ratio",
         title= paste0("Gene-set level analysis for ",
                       disease ,", 95% CI") )+
    theme(legend.position = "none")+
    theme(aspect.ratio=1)
  
  plot2$FDR<-p.adjust(plot2$P,method = 'BH')
  plot2$sign<-ifelse(plot2$FDR<0.05,'*','')
  plot2$VARIABLE<-paste(plot2$VARIABLE,plot2$sign)
  plot2$VARIABLE<-factor(plot2$VARIABLE,levels = c(plot2[order(plot2$or,decreasing = F),]$VARIABLE))
  
  p2[[which(list==i)]]<-
    ggplot(plot2, aes(y=reorder(VARIABLE,or), x=or,  colour=VARIABLE)) + 
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", nrow(plot2),type ='continuous'))+
    geom_vline(xintercept=1, linetype=3) +
    theme_minimal_hgrid(10, rel_small = 1) +
    labs(color = "",y = "", x = "Odds Ratio",
         title= paste0("Gene-set level analysis for ",
                       disease ,", 95% CI") )+
    theme(legend.position = "none")+
    theme(aspect.ratio=1)
  
  plot3$FDR<-p.adjust(plot3$P,method = 'BH')
  plot3$sign<-ifelse(plot3$FDR<0.05,'*','')
  plot3$VARIABLE<-paste(plot3$VARIABLE,plot3$sign)
  plot3$VARIABLE<-factor(plot3$VARIABLE,levels = c(plot3[order(plot3$or,decreasing = F),]$VARIABLE))
  
  p3[[which(list==i)]]<-
    ggplot(plot3, aes(y=reorder(VARIABLE,or), x=or,  colour=VARIABLE)) + 
    geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
    geom_point(size=2)+
    scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", nrow(plot3),type ='continuous'))+
    geom_vline(xintercept=1, linetype=3) +
    theme_minimal_hgrid(10, rel_small = 1) +
    labs(color = "",y = "", x = "Odds Ratio",
         title= paste0("Gene-set level analysis for ",
                       disease ,", 95% CI") )+
    theme(legend.position = "none")+
    theme(aspect.ratio=1)
  
}


pdf("for_list2/hMAGMA_genesetlevel_1.pdf",width = 16,height = 16)
wrap_plots(p,nrow=4) 
dev.off()

pdf("for_list2/hMAGMA_genesetlevel_2.pdf",width = 12,height = 12)
wrap_plots(p2,nrow=4) 
dev.off()

pdf("for_list2/hMAGMA_genesetlevel_3.pdf",width = 12,height = 12)
wrap_plots(p3,nrow=4) 
dev.off()

rmrgenset$fdr<-p.adjust(rmrgenset$P,method = 'BH')
rmrgenset$sign<-ifelse(rmrgenset$fdr<0.05,'*','')
rmrgenset$disease<-paste(rmrgenset$disease,rmrgenset$sign)

pdf("for_list2/hMAGMA_RMRset.pdf",width = 6,height = 6)
ggplot(rmrgenset, aes(y=reorder(disease,or), x=or,  colour=disease)) + 
  geom_errorbarh(aes(xmin=or_lci95, xmax=or_uci95), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 16,type ='continuous'))+
  geom_vline(xintercept=1, linetype=3) +
  theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Odds Ratio",
       title= paste0("Gene-set level analysis for RMR, 95% CI") )+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)
dev.off()

#######


