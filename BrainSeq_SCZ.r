
library(DESeq2)
library(patchwork)

load('brainseq/methprop_pd.Rdata')
load('brainseq/rse_tx_unfiltered.Rdata')
load('brainseq/rse_gene_unfiltered.Rdata')

pheno<-as.data.frame(pd)
pheno<-data.frame(id=c(rownames(pheno)),
                  gender=pheno$Sex,
                  Region=pheno$Region,
                  dx=c(pheno$Dx),
                  ageStage=pheno$ageGroup,
                  MMR=pheno$mean_mitoRate,
                  TSGR=pheno$mean_totalAssignedGene,
                  RIN=c(pheno$mean_RIN),
                  ETHN=c(pheno$Race),
                  age=c(pheno$ageStage))

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)

mod = model.matrix(~Dx + Age + Sex + mitoRate + Region + rRNA_rate + totalAssignedGene + RIN + snpPC1 + snpPC2 +snpPC3 + snpPC4 + snpPC5 + Race,
                   data = colData(rse_gene))


#########
#gene expression

test<-assays(rse_gene)$rpkm

rpkm_to_tpm <- function(rpkm) {
  sumRPKM <- sum(rpkm, na.rm = TRUE)
  tpm <- rpkm / sumRPKM * 1e6
  return(tpm)
}

tpm_matrix <- apply(test, 2, rpkm_to_tpm)
tpm_matrix[1:5,1:5]

# ENSG00000037474|ENSG00000171103

test<-tpm_matrix[grep('ENSG00000037474|ENSG00000122687|ENSG00000144034|ENSG00000135372|ENSG00000104907',rownames(tpm_matrix)),]
test<-log2(test+1)


library(ggpubr)

theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.8, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-ifelse(long_matrix$Var1=='ENSG00000037474.14','NSUN2',
                             ifelse(long_matrix$Var1=='ENSG00000122687.17','MRM2',
                                            ifelse(long_matrix$Var1=='ENSG00000144034.14','TPRKB',
                                                   ifelse(long_matrix$Var1=='ENSG00000135372.8','NAT10','TRMT1'))))

    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=value),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=1)

  }
}

pdf("five_gene_scz.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2, guides="collect")
dev.off()

#森林和meta分析结果

library(ggpubr)

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$name<-ifelse(allresult$name=='ENSG00000037474.14','NSUN2',
                         ifelse(allresult$name=='ENSG00000122687.17','MRM2',
                                ifelse(allresult$name=='ENSG00000144034.14','TPRKB',
                                       ifelse(allresult$name=='ENSG00000135372.8','NAT10','TRMT1'))))


p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',],
         aes(y=reorder(name,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(-3.3,3.3))+
  theme(aspect.ratio=1)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',],
         aes(y=reorder(name,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=1)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',],
         aes(y=reorder(name,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(-3.8,3.8))+
  theme(aspect.ratio=1)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',],
         aes(y=reorder(name,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(-2.2,2.2))+
  theme(aspect.ratio=1)

library(patchwork)

pdf("Brain_seq_SCZ_gene_forest.pdf",width = 8,height = 8)
wrap_plots(p,nrow=2)
dev.off()


library(metafor)

metar<-data.frame()

for (i in unique(allresult$name)){

  forme<-allresult[allresult$name==i,]
  res <- rma(yi = forme$beta,
             sei = (forme$beta_upper - forme$beta_lower) / (2 * 1.96),
             method = "REML")

  cache<- data.frame( name=i,
                      beta=c(res$b[1]),
                      beta_lower=c(res$ci.lb),
                      beta_upper=c(res$ci.ub),
                      pvalue=c(res$pval))

  metar<-rbind(metar,cache)
}


pdf("Brain_seq_SCZ_gene_forest_meta.pdf",width = 4,height = 4)

ggplot(metar,
       aes(y=reorder(name,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("BIP LM (95% CI)") )+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)

dev.off()
#TRMT11的方向都不对了

#对于NSUN2的话最好还是在isoform水平的分析

######################################
#isoform

test<-assays(rse_tx)$tpm

test<-test[grep('ENST00000345094|ENST00000452837|ENST00000515662|ENST00000416603|ENST00000417036|ENST00000441127|ENST00000419437|ENST00000461636|ENST00000464045|ENST00000484006',rownames(test)),]

library(ggpubr)

theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.8, 0.8),#图例在绘图区域的位置
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}

p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-factor(long_matrix$Var1,levels =unique(long_matrix$Var1))


    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "lopg2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=0.5)

  }
}

pdf("Brainseq_THUMPD3.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect")
dev.off()

####森林

test<-assays(rse_tx)$tpm
test<-test[grep('ENST00000345094|ENST00000452837|ENST00000515662|ENST00000416603|ENST00000417036|ENST00000441127|ENST00000419437|ENST00000461636|ENST00000464045|ENST00000484006',rownames(test)),]

pheno<-as.data.frame(pd)
pheno<-data.frame(id=c(rownames(pheno)),
                  gender=pheno$Sex,
                  Region=pheno$Region,
                  dx=c(pheno$Dx),
                  age=c(pheno$ageStage))

p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000345094.7','ENST00000452837.6','ENST00000515662.6'),'Protein coding',
                          ifelse(allresult$name %in% c('ENST00000416603.1','ENST00000417036.1','ENST00000441127.5','ENST00000419437.5'),'CDS incomplete','NMD or Non-coding'))

allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding','CDS incomplete','NMD or Non-coding'))


p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-2.5,2.5))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4.8,4.8))+
  theme(aspect.ratio=0.75)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-3.2,3.2))+
  theme(aspect.ratio=0.75)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-1.8,1.8))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_THUMPD3.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2,guides='collect')
dev.off()

#目前的问题在于这样的分析结果还是很抽象，需要把抽象的指标变成具体的指标

r1<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='Protein coding'],]
r2<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='CDS incomplete'],]
r3<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='NMD or Non-coding'],]

funlos<-data.frame(coding=c(apply(r1,2,sum)),
                   incomplete=c(apply(r2,2,sum)),
                   non=c(apply(r3,2,sum)))
sum<-apply(funlos,1,sum)
funlos <- funlos / rowSums(funlos)
funlos$sum<-sum

allresult<-data.frame()

library(dplyr)
library(tidyr)

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, data = selextexp)

    lmresult<-data.frame(name=c('THUMPD3'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_',j))

    allresult<-rbind(allresult,lmresult)

  }
}

#这个模型和我的设想比较接近

allresult

res <- rma(yi = allresult$beta,
           sei = (allresult$beta_upper - allresult$beta_lower) / (2 * 1.96),
           method = "REML")

cache<- data.frame( name='THUMPD3',
                    beta=c(res$b[1]),
                    beta_lower=c(res$ci.lb),
                    beta_upper=c(res$ci.ub),
                    pvalue=c(res$pval),
                    group=c('meta'))

allresult<-rbind(allresult,cache) #为什么更加显著了，不应该平均一些吗？
allresult$group<-factor(allresult$group,levels = c("meta",'M_HIPPO','M_DLPFC',
                                                   "F_HIPPO","F_DLPFC"))


pdf("Brain_seq_THUMPD3_isoform_function.pdf",width = 4,height = 4)

ggplot(allresult, aes(y=group, x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-2.7,2))+
  theme(aspect.ratio=0.75)

dev.off()

#####
#NSUN2

test<-assays(rse_tx)$tpm
test[1:5,1:5]
test<-test[grep('ENST00000264670|ENST00000506139|ENST00000504374|ENST00000514127|ENST00000505264|ENST00000505892|ENST00000513888|ENST00000507888|ENST00000502932',rownames(test)),]

#box
p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-factor(long_matrix$Var1,levels =unique(long_matrix$Var1))


    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=0.5)

  }
}

pdf("Brainseq_NSUN2_box.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect")
dev.off()

#forest
p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000264670.10','ENST00000506139.5'),'Protein coding',
                          ifelse(allresult$name %in% c('ENST00000504374.5','ENST00000514127.1'),'CDS incomplete','NMD or Non-coding'))

allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding','CDS incomplete','NMD or Non-coding'))


p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-1.8,1.8))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4.8,4.8))+
  theme(aspect.ratio=0.75)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-3.2,3.2))+
  theme(aspect.ratio=0.75)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-1.6,1.6))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_NSUN2_forest.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2,guides='collect')
dev.off()

#function

r1<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='Protein coding'],]
r2<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='CDS incomplete'],]
r3<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='NMD or Non-coding'],]

funlos<-data.frame(coding=c(apply(r1,2,sum)),
                   incomplete=c(apply(r2,2,sum)),
                   non=c(apply(r3,2,sum)))
sum<-apply(funlos,1,sum)
funlos <- funlos / rowSums(funlos)
funlos$sum<-sum

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, data = selextexp)

    lmresult<-data.frame(name=c('NSUN2'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_',j))


    allresult<-rbind(allresult,lmresult)

  }
}

#这个模型和我的设想比较接近

allresult

res <- rma(yi = allresult$beta,
           sei = (allresult$beta_upper - allresult$beta_lower) / (2 * 1.96),
           method = "REML")

cache<- data.frame( name='NSUN2',
                    beta=c(res$b[1]),
                    beta_lower=c(res$ci.lb),
                    beta_upper=c(res$ci.ub),
                    pvalue=c(res$pval),
                    group=c('meta'))

allresult<-rbind(allresult,cache)
allresult$group<-factor(allresult$group,levels = c("meta",'M_HIPPO','M_DLPFC',
                                                   "F_HIPPO","F_DLPFC"))


pdf("Brain_seq_NSUN2_isoform_function_SCZ.pdf",width = 4,height = 4)

ggplot(allresult, aes(y=group, x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4,1))+
  theme(aspect.ratio=0.75)

dev.off()


##

###
#TRMT61B

gc()
test<-assays(rse_tx)$tpm
test<-test[grep('ENST00000306108|ENST00000419999|ENST00000439947|ENST00000484060|ENST00000490390',rownames(test)),]

p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-factor(long_matrix$Var1,levels =unique(long_matrix$Var1))


    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=0.5)

  }
}

pdf("Brainseq_TRMT61B.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect")
dev.off()

p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000306108.9'),'Protein coding',
                          ifelse(allresult$name %in% c('ENST00000419999.5'),'CDS incomplete','NMD or Non-coding'))

allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding','CDS incomplete','NMD or Non-coding'))

p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-1.5,1.5))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4.2,4.2))+
  theme(aspect.ratio=0.75)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-3.4,3.4))+
  theme(aspect.ratio=0.75)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-2.0,2.0))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_TRMT61B_forest.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2,guides='collect')
dev.off()


r1<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='Protein coding'],]
r2<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='CDS incomplete'],]
r3<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='NMD or Non-coding'],]

funlos<-data.frame(coding=c(r1),
                   incomplete=c(r2),
                   non=c(apply(r3,2,sum)))
sum<-apply(funlos,1,sum)
funlos <- funlos / rowSums(funlos)
funlos$sum<-sum

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, data = selextexp)

    lmresult<-data.frame(name=c('NSUN2'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_',j))


    allresult<-rbind(allresult,lmresult)

  }
}

#这个模型和我的设想比较接近

allresult

res <- rma(yi = allresult$beta,
           sei = (allresult$beta_upper - allresult$beta_lower) / (2 * 1.96),
           method = "REML")

cache<- data.frame( name='THUMPD3',
                    beta=c(res$b[1]),
                    beta_lower=c(res$ci.lb),
                    beta_upper=c(res$ci.ub),
                    pvalue=c(res$pval),
                    group=c('meta'))

allresult<-rbind(allresult,cache) #为什么更加显著了，不应该平均一些吗？
allresult$group<-factor(allresult$group,levels = c("meta",'M_HIPPO','M_DLPFC',
                                                   "F_HIPPO","F_DLPFC"))


pdf("Brain_seq_TRMT61B_isoform_function.pdf",width = 4,height = 4)

ggplot(allresult, aes(y=group, x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-2.5,1.2))+
  theme(aspect.ratio=0.75)

dev.off()

#现在的模型正常多了

######
#把剩下几个sQTl的基因也看了得了

###
gc()
test<-assays(rse_tx)$tpm
test<-test[grep('ENST00000199320|ENST00000506390|ENST00000514911|ENST00000509182|ENST00000514605',rownames(test)),]


p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-factor(long_matrix$Var1,levels =unique(long_matrix$Var1))


    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=0.5)

  }
}

pdf("Brainseq_DIMT1_box.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect")
dev.off()

p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000199320.8','ENST00000506390.5'),
                          'Protein coding','NMD or Non-coding')

allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding','CDS incomplete','NMD or Non-coding'))

p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4,4))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-5,5))+
  theme(aspect.ratio=0.75)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4,4))+
  theme(aspect.ratio=0.75)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4,4))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_DIMT1_forest.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2,guides='collect')
dev.off()

r1<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='Protein coding'],]
r2<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='CDS incomplete'],]
r3<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='NMD or Non-coding'],]

funlos<-data.frame(coding=c(apply(r1,2,sum)),
                   incomplete=c(0),
                   non=c(apply(r3,2,sum)))
sum<-apply(funlos,1,sum)
funlos <- funlos / rowSums(funlos)
funlos$sum<-sum

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, data = selextexp)

    lmresult<-data.frame(name=c('NSUN2'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_',j))


    allresult<-rbind(allresult,lmresult)

  }
}

#这个模型和我的设想比较接近

allresult

res <- rma(yi = allresult$beta,
           sei = (allresult$beta_upper - allresult$beta_lower) / (2 * 1.96),
           method = "REML")

cache<- data.frame( name='THUMPD3',
                    beta=c(res$b[1]),
                    beta_lower=c(res$ci.lb),
                    beta_upper=c(res$ci.ub),
                    pvalue=c(res$pval),
                    group=c('meta'))

allresult<-rbind(allresult,cache) #为什么更加显著了，不应该平均一些吗？
allresult$group<-factor(allresult$group,levels = c("meta",'M_HIPPO','M_DLPFC',
                                                   "F_HIPPO","F_DLPFC"))


pdf("Brain_seq_DIMT1_isoform_function.pdf",width = 4,height = 4)

ggplot(allresult, aes(y=group, x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-3,1.5))+
  theme(aspect.ratio=0.75)

dev.off()

####
#DTWD1

gc()
test<-assays(rse_tx)$tpm
test<-test[grep('ENST00000251250|ENST00000415425|ENST00000403028|ENST00000558653|ENST00000559164|ENST00000559405|ENST00000560632|ENST00000329873|ENST00000557988|ENST00000557968|ENST00000559223|ENST00000560735|ENST00000561188',rownames(test)),]

p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-factor(long_matrix$Var1,levels =unique(long_matrix$Var1))


    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=0.5)

  }
}

pdf("Brainseq_DTWD1_box.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect")
dev.off()


p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000251250.7','ENST00000403028.7',
                                                'ENST00000558653.5'),'Protein coding',
                          ifelse(allresult$name %in% c('ENST00000559164.5','ENST00000559405.5',
                                                       'ENST00000560632.5'),'CDS incomplete',
                                 'NMD or Non-coding'))


allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding',
                                                       'CDS incomplete','NMD or Non-coding'))

p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4,4))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=0.75)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=0.75)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_DTWD1_forest.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2,guides='collect')
dev.off()


r1<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='Protein coding'],]
r2<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='CDS incomplete'],]
r3<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='NMD or Non-coding'],]

funlos<-data.frame(coding=c(apply(r1,2,sum)),
                   incomplete=c(apply(r2,2,sum)),
                   non=c(apply(r3,2,sum)))
sum<-apply(funlos,1,sum)
funlos <- funlos / rowSums(funlos)
funlos$sum<-sum

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, data = selextexp)

    lmresult<-data.frame(name=c('NSUN2'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_',j))


    allresult<-rbind(allresult,lmresult)

  }
}

allresult

res <- rma(yi = allresult$beta,
           sei = (allresult$beta_upper - allresult$beta_lower) / (2 * 1.96),
           method = "REML")

cache<- data.frame( name='THUMPD3',
                    beta=c(res$b[1]),
                    beta_lower=c(res$ci.lb),
                    beta_upper=c(res$ci.ub),
                    pvalue=c(res$pval),
                    group=c('meta'))

allresult<-rbind(allresult,cache) #为什么更加显著了，不应该平均一些吗？
allresult$group<-factor(allresult$group,levels = c("meta",'M_HIPPO','M_DLPFC',
                                                   "F_HIPPO","F_DLPFC"))


pdf("Brain_seq_DTWD1_isoform_function.pdf",width = 4,height = 4)

ggplot(allresult, aes(y=group, x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  theme(aspect.ratio=0.75)

dev.off()

############################


gc()
test<-assays(rse_tx)$tpm
test<-test[grep('ENST00000280665|ENST00000397173|ENST00000540622|ENST00000543381|ENST00000541700|ENST00000536665|ENST00000537725|ENST00000535873',rownames(test)),]

p<-list()

for (i in c('M','F')){
  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-factor(long_matrix$Var1,levels =unique(long_matrix$Var1))


    library(ggpubr)

    p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' ',j))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
      theme(aspect.ratio=0.5)

  }
}

pdf("Brainseq_DCP1B_box.pdf",width = 18,height = 12)
wrap_plots(p,nrow=2, guides="collect")
dev.off()


p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_',j))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }
}

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000280665.10'),'Protein coding',
                          ifelse(allresult$name %in% c('ENST00000540622.1'),'CDS incomplete',
                                 'NMD or Non-coding'))


allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding',
                                                       'CDS incomplete','NMD or Non-coding'))

p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-4,4))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_HIPPO',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F HIPPO LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=0.75)

p[[4]]<-
  ggplot(allresult[allresult$group=='F_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt", 3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=0.75)

p[[3]]<-
  ggplot(allresult[allresult$group=='M_DLPFC',], aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-6,6))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_DCP1B_forest.pdf",width = 10,height = 8)
wrap_plots(p,nrow=2,guides='collect')
dev.off()


r1<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='Protein coding'],]
r2<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='CDS incomplete'],]
r3<-test[rownames(test) %in% allresult$name[allresult$Biotypr=='NMD or Non-coding'],]

funlos<-data.frame(coding=c(r1),
                   incomplete=c(r2),
                   non=c(apply(r3,2,sum)))
sum<-apply(funlos,1,sum)
funlos <- funlos / rowSums(funlos)
funlos$sum<-sum

allresult<-data.frame()

for (i in c('M','F')){

  for (j in c('HIPPO','DLPFC')){

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$age=='Adult',]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(mod[match(rownames(selextexp),rownames(mod)),])

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$RaceAS+corv$RaceCAUC+corv$RaceHISP, data = selextexp)

    lmresult<-data.frame(name=c('NSUN2'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_',j))


    allresult<-rbind(allresult,lmresult)

  }
}

allresult

res <- rma(yi = allresult$beta,
           sei = (allresult$beta_upper - allresult$beta_lower) / (2 * 1.96),
           method = "REML")

cache<- data.frame( name='THUMPD3',
                    beta=c(res$b[1]),
                    beta_lower=c(res$ci.lb),
                    beta_upper=c(res$ci.ub),
                    pvalue=c(res$pval),
                    group=c('meta'))

allresult<-rbind(allresult,cache) #为什么更加显著了，不应该平均一些吗？
allresult$group<-factor(allresult$group,levels = c("meta",'M_HIPPO','M_DLPFC',
                                                   "F_HIPPO","F_DLPFC"))


pdf("Brain_seq_DCP1B_isoform_function.pdf",width = 4,height = 4)

ggplot(allresult, aes(y=group, x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M DLPFC LM (95% CI)") )+
  theme(legend.position = "right")+
  theme(aspect.ratio=0.75)

dev.off()
