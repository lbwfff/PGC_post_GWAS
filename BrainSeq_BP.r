#然后就是在BIP上把对于SCZ做的分析如法炮制

###我还是觉得对于BIP的统计可以捡起来，哪怕是阴性的结果也无所谓

library(DESeq2)

load('caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda')
phase3ge<-assays(rse_gene)$counts

phase3meta = data.frame(sample=c(rownames(colData(rse_gene))),
                        RIN=c(sapply(colData(rse_gene)$RIN,"[",1)),
                        totalAssignedGene=c(sapply(colData(rse_gene)$totalAssignedGene, mean)),
                        mitoRate=c(sapply(colData(rse_gene)$mitoRate,mean)),
                        overallMapRate=c(sapply(colData(rse_gene)$overallMapRate, mean)),
                        rRNA_rate=c(sapply(colData(rse_gene)$rRNA_rate,mean)),
                        ERCCsumLogErr=c(sapply(colData(rse_gene)$ERCCsumLogErr,mean)),
                        Region=c(colData(rse_gene)$Region),
                        Race=c(colData(rse_gene)$Race),
                        Age=c(colData(rse_gene)$Age),
                        Sex=c(colData(rse_gene)$Sex),
                        Dx=c(colData(rse_gene)$Dx))

lib_size <- colSums(phase3ge)
tpm_matrix  <- sweep(phase3ge, 2, lib_size, "/") * 1e6
tpm_matrix[1:5,1:5]

test<-tpm_matrix[grep('ENSG00000037474|ENSG00000122687',rownames(tpm_matrix)),]
test<-log2(test+1)
pheno<-as.data.frame(phase3meta)
pheno<-data.frame(id=c(pheno$sample),
                  gender=pheno$Sex,
                  Region=pheno$Region,
                  dx=c(pheno$Dx),
                  age=c(pheno$Age))

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

    selectsample<-pheno[pheno$gender==i & pheno$dx!='Schizo' & pheno$age>18,]
    selextexp<-test[,colnames(test) %in% selectsample$id]
    long_matrix <- reshape2::melt(selextexp)
    match<-selectsample[match(long_matrix$Var2,selectsample$id),]
    long_matrix$group<-match$dx

    long_matrix$Var1<-ifelse(long_matrix$Var1=='ENSG00000037474.14','NSUN2','MRM2')
    long_matrix$group<-factor(long_matrix$group,level=c('Control','Bipolar'))

    library(ggpubr)

    p[[paste0(i)]] <- ggplot(long_matrix,aes(x=Var1,y=value),palette = "jco", add = "jitter")+
      geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
      labs(x = NULL, y = "log2(TPM)",title = paste0(i,' Caudate'))+
      scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
      theme_zg()+ stat_compare_means(aes(group = group),
                                     label = "p.signif",vjust=0.5,
                                     method='wilcox.test')+
      theme(aspect.ratio=1)


}

pdf("two_gene_BIP.pdf",width = 8,height = 6)
wrap_plots(p,nrow=1, guides="collect")
dev.off()

#####################

#然后就是森林图，然后就是isoform水平的分析

library(ggpubr)

p<-list()

allresult<-data.frame()

for (i in c('M','F')){

    selectsample<-pheno[pheno$gender==i & pheno$Region=='Caudate' & pheno$age>18,]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(phase3meta[match(rownames(selextexp),phase3meta$sample),])
    corv$Race<-ifelse(corv$Race=='AA',0,1)


    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$Race, family = binomial)
      summary(model1)$coefficients

      cache<-data.frame(name=c(colnames(selextexp)[k]),
                        beta=c(coef(model1)[2]),
                        beta_lower=confint(model1,level = 0.95)[2,1],
                        beta_upper=confint(model1,level = 0.95)[2,2],
                        pvalue=c(summary(model1)$coefficients[2,4]),
                        group=paste0(i,'_Caudate'))

      lmresult<-rbind(lmresult,cache)
    }

    allresult<-rbind(allresult,lmresult)

  }

allresult$name<-ifelse(allresult$name=='ENSG00000122687.17','MRM2','NSUN2')
allresult$fdr<-p.adjust(allresult$pvalue,method = 'BH')
allresult$sign<-ifelse(allresult$fdr<0.05,'*','')
allresult$usname<-paste0(allresult$name,allresult$sign)

p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_Caudate',],
         aes(y=reorder(usname,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M Caudate LM (95% CI)") )+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_Caudate',],
         aes(y=reorder(usname,beta), x=beta)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("VanGogh1", 2,type ='continuous'))+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F Caudate LM (95% CI)") )+
  theme(legend.position = "none")+
  theme(aspect.ratio=1)

library(patchwork)

pdf("Brain_seq3_BIP.pdf",width = 12,height = 20)
wrap_plots(p,nrow=1)
dev.off()

#做一个meta分析把，对于SCZ也可以这样做

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


pdf("Brain_seq_BIP_gene_forest_meta.pdf",width = 4,height = 4)
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

#########################
#isoform 水平的分析
#只用看NSUN2了

load('caudate_brainseq_phase3_hg38_rseTx_merged_n464.rda')
phase3ge<-assays(rse_tx)$tpm


test<-phase3ge[grep('ENST00000264670|ENST00000506139|ENST00000504374|ENST00000514127|ENST00000505264|ENST00000505892|ENST00000513888|ENST00000507888|ENST00000502932',rownames(phase3ge)),]
test<-log2(test+1)
pheno<-as.data.frame(phase3meta)
pheno<-data.frame(id=c(pheno$sample),
                  gender=pheno$Sex,
                  Region=pheno$Region,
                  dx=c(pheno$Dx),
                  age=c(pheno$Age))

library(ggpubr)

p<-list()

for (i in c('M','F')){
  j<-'Caudate'
  selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$dx!='Schizo' & pheno$age>18,]
  selextexp<-test[,colnames(test) %in% selectsample$id]
  long_matrix <- reshape2::melt(selextexp)
  match<-selectsample[match(long_matrix$Var2,selectsample$id),]
  long_matrix$group<-match$dx

  long_matrix$Var1<-factor(long_matrix$Var1,levels =c('ENST00000264670.10','ENST00000506139.5','ENST00000504374.5','ENST00000514127.1','ENST00000505264.1','ENST00000505892.5','ENST00000513888.5','ENST00000507888.1','ENST00000502932.1'))
  long_matrix$group<-factor(long_matrix$group,level=c('Control','Bipolar'))


  library(ggpubr)

  p[[paste0(i,'_',j)]] <- ggplot(long_matrix,aes(x=Var1,y=log2(value+1)),palette = "jco", add = "jitter")+
    geom_boxplot(aes(fill=group),position=position_dodge(0.5),width=0.6)+
    labs(x = NULL, y = "lopg2(TPM)",title = paste0(i,' ',j))+
    scale_fill_manual(values = c("#66c2a5",'#f46d43')) +
    theme_zg()+ stat_compare_means(aes(group = group),label = "p.signif",vjust=0.5,method='wilcox.test')+
    theme(aspect.ratio=0.5)


}

pdf("NSUN2_isoform_BIP.pdf",width = 14,height = 8)
wrap_plots(p,nrow=1, guides="collect")
dev.off()

#森林

p<-list()

allresult<-data.frame()

for (i in c('M','F')){

  j<-'Caudate'

    selectsample<-pheno[pheno$gender==i & pheno$Region==j & pheno$dx!='Schizo' & pheno$age>18,]
    selextexp<-t(test[,colnames(test) %in% selectsample$id])

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    corv<-as.data.frame(phase3meta[match(rownames(selextexp),phase3meta$sample),])
    corv$Race<-ifelse(corv$Race=='AA',0,1)

    lmresult<-data.frame()

    for (k in 1:ncol(selextexp)) {
      model1 <-glm(ifelse(selectsample$dx=='Control',0,1) ~ selextexp[,k]+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$Race, family = binomial)
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

allresult$Biotypr<-ifelse(allresult$name %in% c('ENST00000264670.10','ENST00000506139.5'),'Protein coding',
                          ifelse(allresult$name %in% c('ENST00000504374.5','ENST00000514127.1'),'CDS incomplete','NMD or Non-coding'))

allresult$Biotypr<-factor(allresult$Biotypr,levels = c('Protein coding','CDS incomplete','NMD or Non-coding'))


p<-list()

p[[1]]<-
  ggplot(allresult[allresult$group=='M_Caudate',],
         aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("M Caudate LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-3.5,3.5))+
  theme(aspect.ratio=0.75)

p[[2]]<-
  ggplot(allresult[allresult$group=='F_Caudate',],
         aes(y=reorder(name,beta), x=beta,  colour=Biotypr)) +
  geom_errorbarh(aes(xmin=beta_lower, xmax=beta_upper), height=.3) +
  geom_point(size=2)+
  scale_color_manual(values=MetBrewer::met.brewer("Egypt",3,type ='continuous')[c(2,3,1)])+
  geom_vline(xintercept=0, linetype=3) +
  cowplot::theme_minimal_hgrid(10, rel_small = 1) +
  labs(color = "",y = "", x = "Beta value",
       title= paste0("F Caudate LM (95% CI)") )+
  theme(legend.position = "right")+
  coord_cartesian(xlim=c(-8,8))+
  theme(aspect.ratio=0.75)

library(patchwork)

pdf("Brain_seq_NSUN2isoform_BIP.pdf",width = 10,height = 4)
wrap_plots(p,nrow=1,guides='collect')
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

library(dplyr)
library(tidyr)

for (i in c('M','F')){

    selectsample<-pheno[pheno$gender==i &pheno$dx!='Schizo' & pheno$age>18,]
    selextexp<-funlos[rownames(funlos) %in% selectsample$id,]

    selectsample<-selectsample[match(rownames(selextexp),selectsample$id),]
    selextexp$group<-selectsample$dx
    corv<-as.data.frame(phase3meta[match(rownames(selextexp),phase3meta$sample),])
    corv$Race<-ifelse(corv$Race=='AA',0,1)

    selextexp$severity_score <-
      (selextexp$coding + 0.5*selextexp$incomplete)

    model1<-lm(ifelse(selectsample$dx=='Control',0,1) ~ severity_score + sum+corv$mitoRate+corv$Age+corv$rRNA_rate+corv$totalAssignedGene+corv$RIN+corv$Race, data = selextexp)

    lmresult<-data.frame(name=c('THUMPD3'),
                         beta=c(coef(model1)[2]),
                         beta_lower=confint(model1,level = 0.95)[2,1],
                         beta_upper=confint(model1,level = 0.95)[2,2],
                         pvalue=c(summary(model1)$coefficients[2,4]),
                         group=paste0(i,'_','Caudate'))

    allresult<-rbind(allresult,lmresult)


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
allresult$group<-factor(allresult$group,levels = c("meta",'M_Caudate',
                                                   "F_Caudate"))


pdf("Brain_seq_NSUN2_isoform_function.pdf",width = 4,height = 4)

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


