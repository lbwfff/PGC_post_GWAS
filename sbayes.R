
empirical_p <- function(D, x) {
  D <- D[!is.na(D)]
  n <- length(D)
  p_high <- (sum(D > x) + 1) / (n + 1)
  p_low  <- (sum(D < x) + 1) / (n + 1)
  p_two  <- min(1, 2 * min(p_high, p_low))
  list(p_high = p_high, p_low = p_low, p_two = p_two)
}

#遗传力的分析
#brainmeta给的基因有点太复杂化了，是不是缩减一些会好一些？

library(data.table)
library(GenomicRanges)

sczh2<-fread('./sbayes/sbayesr_output.snpRes')

#hg19的坐标

load('./sbayes/geneAnno_allgenes.rda')
gene<-read.table("./sbayes/Gencode26_gene.bed")
promoter<-read.table("./sbayes/Gencode26_promoter.bed")

generanges <- GRanges(c(gene[,1],paste0('chr',promoter[,1])),
                      IRanges(c(gene[,2],promoter[,2]),c(gene[,3],promoter[,3])),
                      gene=c(gene[,4],promoter[,4]))

grl <- split(generanges, generanges$gene)
gr_span <- range(grl) 
keep <- elementNROWS(gr_span) == 1L
gr_span1 <- gr_span[keep]
gr_clean <- unlist(gr_span1, use.names = FALSE)
gr_clean$gene <- names(gr_span1)
generanges<-gr_clean

snps <- GRanges(paste0('chr',sczh2$Chrom), IRanges(sczh2$Position, sczh2$Position), rsid=sczh2$Name)

olap <- findOverlaps(snps,generanges)
gene_to_snps <- split(queryHits(olap), subjectHits(olap))

geneh2<-data.frame()

for (i in 1:length(gene_to_snps)){
  index<-as.numeric(names(gene_to_snps)[i])
  id<-generanges$gene[index]
  gene<-NA
  
  snplist<-sczh2[gene_to_snps[[i]],]
  length<-nrow(snplist)
  mergh2<-sum(snplist$VarExplained)
  
  cache<-data.frame(ID=c(id),
             GENE=c(gene),
             H2=c(mergh2),
             length=c(length))
  geneh2<-rbind(geneh2,cache)
}


library(rtracklayer)

# gtf <- data.frame(import.gff('../biodata/gencode.v47.annotation.gtf',
#                              format = 'gtf')) 
# gtf<-gtf[gtf$type=='gene',]
# gtf<-data.frame(id=c(gtf$gene_id),
#                 gene=c(gtf$gene_name))
save(gtf,file = 'i2g.RData')
match<-gtf[match(geneh2$ID,substr(gtf$id,1,15)),]
geneh2$GENE<-match$gene

sczh2<-geneh2

library(dplyr)
library(ggplot2)
library(ggrepel)

sczh2$pve<-sczh2$H2/sum(sczh2$H2)

modlist<-read.csv('geneset2.csv')
modlist<-modlist[modlist$Gene.Symbol !='',]

dens <- density(log2(sczh2$pve), na.rm = TRUE)

lab_df <- sczh2 %>%
  filter(GENE %in% modlist$Gene.Symbol,
         pve > quantile(pve, 0.95, na.rm = TRUE)) %>%
  mutate(
    log2pve = log2(pve),
    dens_y = approx(dens$x, dens$y, xout = log2pve, rule = 2)$y
  )

risklist<-read.csv('INT-18_SCZ_Risk_Gene_List.csv')
sczh2risk<-sczh2[sczh2$GENE %in% risklist$sczgenenames,]
sczh2risk<-sczh2risk[!is.na(sczh2risk$GENE),]


pdf('SCZ_gene_h2.pdf',width = 6,height = 5)
ggplot(sczh2, aes(x = log2(pve))) +
  geom_density(aes(y = after_stat(density)),
               adjust = 1/2, size = 1) +
  geom_density(data = sczh2risk, aes(x = log2(pve)), adjust = 1/2,
               color = "#96bfce", size = 0.5)+
  theme_classic(base_size = 12) +
  labs(title = "SCZ",
       x = "-log2 (Relative heritability contribution)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        aspect.ratio = 0.75) +
  geom_vline(xintercept = log2(quantile(sczh2$pve, 0.95, na.rm = TRUE)),
             linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(
    data = lab_df,aes(x = log2pve, y = dens_y, label = GENE),
    color = "#0f7bbb",size = 3,
    segment.color = "black",
    show.legend = FALSE,nudge_y = 0.02,nudge_x = 2, 
    box.padding = 0.5,point.padding = 0.5)
dev.off()

write.csv(sczh2,file='sczh2_result.csv')

##########################

biph2<-fread('./sbayes/sbayesr_BIPI.snpRes')


snps <- GRanges(paste0('chr',biph2$Chrom), IRanges(biph2$Position, biph2$Position), rsid=sczh2$Name)

olap <- findOverlaps(snps,generanges)
gene_to_snps <- split(queryHits(olap), subjectHits(olap))

geneh2<-data.frame()

for (i in 1:length(gene_to_snps)){
  index<-as.numeric(names(gene_to_snps)[i])
  id<-generanges$gene[index]
  gene<-NA
  
  snplist<-biph2[gene_to_snps[[i]],]
  mergh2<-sum(snplist$VarExplained)
  length<-nrow(snplist)
  
  cache<-data.frame(ID=c(id),
                    GENE=c(gene),
                    H2=c(mergh2),
                    length=length)
  geneh2<-rbind(geneh2,cache)
}

match<-gtf[match(geneh2$ID,substr(gtf$id,1,15)),]
geneh2$GENE<-match$gene

geneh2$pve<-geneh2$H2/sum(geneh2$H2)

library(ggplot2)

dens <- density(log2(geneh2$pve), na.rm = TRUE)

lab_df <- geneh2 %>%
  filter(GENE %in% modlist$Gene.Symbol,
         pve > quantile(pve, 0.95, na.rm = TRUE)) %>%
  mutate(
    log2pve = log2(pve),
    dens_y = approx(dens$x, dens$y, xout = log2pve, rule = 2)$y
  )

risklist<-read.csv('bipolar_disorder_gandal_et_al.csv',header = F)
risklist<-as.character(risklist[1,3:1120])
geneh2risk<-geneh2[geneh2$GENE %in% risklist,]
geneh2risk<-geneh2risk[!is.na(geneh2risk$GENE),]

pdf('BIP_gene_h2.pdf',width = 6,height = 5)
ggplot(geneh2, aes(x = log2(pve))) +
  geom_density(aes(y = after_stat(density)),
               adjust = 1/2, size = 1) +
  geom_density(data = geneh2risk, aes(x = log2(pve)), adjust = 1/2,
               color = "#96bfce", size = 0.5)+
  theme_classic(base_size = 12) +
  labs(title = "SCZ",
       x = "-log2 (Relative heritability contribution)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        aspect.ratio = 0.75) +
  geom_vline(xintercept = log2(quantile(geneh2$pve, 0.95, na.rm = TRUE)),
             linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(
    data = lab_df,aes(x = log2pve, y = dens_y, label = GENE),
    color = "#0f7bbb",size = 3,
    segment.color = "black",
    show.legend = FALSE,nudge_y = 0.06,nudge_x = 2, 
    box.padding = 0.5,point.padding = 0.5)
dev.off()

write.csv(geneh2,file='BIPIh2_result.csv')

idx <- which(geneh2$GENE == "NSUN2")
rank(-(geneh2$pve))[idx]

##############################
#然后我想做一下随机抽样的结果

sczsnph2<-fread('./sbayes/sbayesr_output.snpRes')

premuh2<-data.frame()

for (i in 1:10000){
  
  sample<-sample(1:nrow(sczsnph2),37,replace=F)
  premu<-sczsnph2[sample,]
  
  length<-nrow(premu)
  mergh2<-sum(premu$VarExplained)
  
  cache<-data.frame(ID=c(paste0('premutation_',i)),
                    GENE=c(gene),
                    H2=c(mergh2),
                    length=c(length))
  
  premuh2<-rbind(premuh2,cache)
}

ptext<-paste0('Empirical P high = ',
              round(as.numeric(empirical_p(premuh2$H2, 
                                           sczh2$H2[sczh2$GENE=='NSUN2' & !is.na(sczh2$GENE)])[1]),3))


pdf('SCZ_premuta_h2.pdf',width = 6,height = 5)
ggplot(premuh2, aes(log2(H2), after_stat(density))) +
  geom_density(adjust = 1/2, size=1)+ 
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(size=12))+
  theme_bw() +
  theme_classic(base_size = 12) + 
  labs(title = paste0('SCZ'),
       x = '-log2 (Relative heritability contribution)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
        legend.position = "None",  
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10)  ) +
  geom_vline(xintercept = log2(sczh2$H2[sczh2$GENE=='NSUN2']), 
             linetype = "dashed", color = "black", size = 0.5)+
  theme(aspect.ratio=0.75)+
  annotate('text',x=-14,y=1.4,label=ptext,colour = "black",size=4)
dev.off()



biph2<-fread('./sbayes/sbayesr_BIPI.snpRes')

premuh2<-data.frame()

for (i in 1:10000){
  
  sample<-sample(1:nrow(biph2),37,replace=F)
  premu<-biph2[sample,]
  
  length<-nrow(premu)
  mergh2<-sum(premu$VarExplained)
  
  cache<-data.frame(ID=c(paste0('premutation_',i)),
                    GENE=c(gene),
                    H2=c(mergh2),
                    length=c(length))
  
  premuh2<-rbind(premuh2,cache)
}

ptext<-paste0('Empirical P high = ',
              round(as.numeric(empirical_p(premuh2$H2, 
                                           geneh2$H2[geneh2$GENE=='NSUN2' & !is.na(geneh2$GENE)])[1]),3))


pdf('BIPI_premuta_h2.pdf',width = 6,height = 5)
ggplot(premuh2, aes(log2(H2), after_stat(density))) +
  geom_density(adjust = 1/2, size=1)+ 
  theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(size=12))+
  theme_bw() +
  theme_classic(base_size = 12) + 
  labs(title = paste0('BIPI'),
       x = '-log2 (heritability)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  
        legend.position = "None",  
        legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10)  ) +
  geom_vline(xintercept = log2(sczh2$H2[sczh2$GENE=='NSUN2']), 
             linetype = "dashed", color = "black", size = 0.5)+
  theme(aspect.ratio=0.75)+
  annotate('text',x=-16,y=1.4,label=ptext,colour = "black",size=4)
dev.off()


##################
#AD 和ADHD

library(data.table)

ADh2<-fread('./sbayes/sbayesr_AD.snpRes')

snps <- GRanges(paste0('chr',ADh2$Chrom), IRanges(ADh2$Position, ADh2$Position), rsid=ADh2$Name)

olap <- findOverlaps(snps,generanges)
gene_to_snps <- split(queryHits(olap), subjectHits(olap))

geneh2<-data.frame()

for (i in 1:length(gene_to_snps)){
  index<-as.numeric(names(gene_to_snps)[i])
  id<-generanges$gene[index]
  gene<-NA
  
  snplist<-ADh2[gene_to_snps[[i]],]
  mergh2<-sum(snplist$VarExplained)
  length<-nrow(snplist)
  
  cache<-data.frame(ID=c(id),
                    GENE=c(gene),
                    H2=c(mergh2),
                    length=length)
  geneh2<-rbind(geneh2,cache)
}

match<-gtf[match(geneh2$ID,substr(gtf$id,1,15)),]
geneh2$GENE<-match$gene

geneh2$pve<-geneh2$H2/sum(geneh2$H2)

library(ggplot2)

dens <- density(log2(geneh2$pve), na.rm = TRUE)

lab_df <- geneh2 %>%
  filter(GENE %in% modlist$Gene.Symbol,
         pve > quantile(pve, 0.95, na.rm = TRUE)) %>%
  mutate(
    log2pve = log2(pve),
    dens_y = approx(dens$x, dens$y, xout = log2pve, rule = 2)$y
  )



pdf('AD_gene_h2.pdf',width = 6,height = 5)
ggplot(geneh2, aes(x = log2(pve))) +
  geom_density(aes(y = after_stat(density)),
               adjust = 1/2, size = 1) +
  theme_classic(base_size = 12) +
  labs(title = "AD",
       x = "-log2 (Relative heritability contribution)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        aspect.ratio = 0.75) +
  geom_vline(xintercept = log2(quantile(geneh2$pve, 0.95, na.rm = TRUE)),
             linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(
    data = lab_df,aes(x = log2pve, y = dens_y, label = GENE),
    color = "#0f7bbb",size = 3,
    segment.color = "black",
    show.legend = FALSE,nudge_y = 0.06,nudge_x = 2, 
    box.padding = 0.5,point.padding = 0.5)
dev.off()


##

ADh2<-fread('./sbayes/sbayesr_ADHD.snpRes')

snps <- GRanges(paste0('chr',ADh2$Chrom), IRanges(ADh2$Position, ADh2$Position), rsid=ADh2$Name)

olap <- findOverlaps(snps,generanges)
gene_to_snps <- split(queryHits(olap), subjectHits(olap))

geneh2<-data.frame()

for (i in 1:length(gene_to_snps)){
  index<-as.numeric(names(gene_to_snps)[i])
  id<-generanges$gene[index]
  gene<-NA
  
  snplist<-ADh2[gene_to_snps[[i]],]
  mergh2<-sum(snplist$VarExplained)
  length<-nrow(snplist)
  
  cache<-data.frame(ID=c(id),
                    GENE=c(gene),
                    H2=c(mergh2),
                    length=length)
  geneh2<-rbind(geneh2,cache)
}

match<-gtf[match(geneh2$ID,substr(gtf$id,1,15)),]
geneh2$GENE<-match$gene

geneh2$pve<-geneh2$H2/sum(geneh2$H2)

library(ggplot2)

dens <- density(log2(geneh2$pve), na.rm = TRUE)

lab_df <- geneh2 %>%
  filter(GENE %in% modlist$Gene.Symbol,
         pve > quantile(pve, 0.95, na.rm = TRUE)) %>%
  mutate(
    log2pve = log2(pve),
    dens_y = approx(dens$x, dens$y, xout = log2pve, rule = 2)$y
  )



pdf('ADHD_gene_h2.pdf',width = 6,height = 5)
ggplot(geneh2, aes(x = log2(pve))) +
  geom_density(aes(y = after_stat(density)),
               adjust = 1/2, size = 1) +
  theme_classic(base_size = 12) +
  labs(title = "ADHD",
       x = "-log2 (Relative heritability contribution)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "none",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        aspect.ratio = 0.75) +
  geom_vline(xintercept = log2(quantile(geneh2$pve, 0.95, na.rm = TRUE)),
             linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(
    data = lab_df,aes(x = log2pve, y = dens_y, label = GENE),
    color = "#0f7bbb",size = 3,
    segment.color = "black",
    show.legend = FALSE,nudge_y = 0.06,nudge_x = 2, 
    box.padding = 0.5,point.padding = 0.5)
dev.off()


