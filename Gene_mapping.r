library(dplyr)
library(geni.plots)

list<-list.files('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/')

genelist<-read.csv('/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/geneset2.csv')

load('/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/geneAnno_allgenes.rda')
geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol %in% genelist$Gene.Symbol,]
geneAnno1<-geneAnno1[geneAnno1$hgnc_symbol!='',]

exon <- read.table("/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/Gencode26_exon.bed")
exon$V1 <- sub("^", "chr", exon$V1)

promoter <- read.table("/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/Gencode26_promoter.bed")
promoter$V1 <- sub("^", "chr", promoter$V1)

exonranges <- GRanges(exon[,1],IRanges(exon[,2],exon[,3]),gene=exon[,4])
promoterranges <- GRanges(promoter[,1], IRanges(promoter[,2], promoter[,3]), gene=promoter[,4])

library(GenomicRanges)

grl <- split(exonranges, exonranges$gene)
gene_span <- unlist(range(grl, ignore.strand = TRUE))
mcols(gene_span)$gene <- names(gene_span)
names(gene_span) <- NULL
gene_span <-gene_span[gene_span$gene %in% geneAnno1$ensembl_gene_id]

for (i in list){
  
  result<-fread(paste0('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/',i))
  
  snps <- result[,c(2,1,3)]
  colnames(snps) <- c("chr","SNP","Position")
  snps$chr <- sub("^", "chr", snps$chr)
  
  snps <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$SNP)
  
  olap <- findOverlaps(snps,gene_span)
  snpgene <- snps[queryHits(olap)]
  
  result<-result[result$SNP %in% snpgene$rsid,]
  fil<-result[result$CHR==6,]
  fil<-fil[fil$BP > 28477797 & fil$BP< 33448354,]
  result<-result[!(result$SNP %in% fil$SNP),]
  
  manhattan<-data.frame(chr=c(result$CHR),
                        pos=c(result$BP),
                        pvalue=c(result$P),
                        label=c(""))
  manhattanp<-manhattan[manhattan$pvalue<0.8,]
  p<-list()
  
  p[[1]]<-
  fig_manhattan(
    data = manhattanp,
    thresh = c(1e-06, 5e-08),
  )  +theme(aspect.ratio=0.5)
  
  p[[2]]<-
  fig_qq(
    pvalues = manhattan$pvalue
  ) +theme(aspect.ratio=1)
  
  path<-paste0('/scratch/lb4489/project/GWAS/',i,'ma&qq.pdf')
  
  pdf(path,width = 15,height = 6)
  print(wrap_plots(p,col=2) )
  dev.off()
  
  gc()
}

allrmr<-as.data.frame(gene_span)

####
#for SCZ
i<-list[15]

result<-fread(paste0('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/',i))

snps <- result[,c(2,1,3)]
colnames(snps) <- c("chr","SNP","Position")
snps$chr <- sub("^", "chr", snps$chr)

snps <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$SNP)

olap <- findOverlaps(snps,gene_span)
snpgene <- snps[queryHits(olap)]

result<-result[match( snpgene$rsid,result$SNP),]
result$gene<-gene_span$gene[subjectHits(olap)]
match<-geneAnno1[match(result$gene,geneAnno1$ensembl_gene_id),]
result$gene<-match$hgnc_symbol

fil<-result[result$CHR==6,]
fil<-fil[fil$BP > 28477797 & fil$BP< 33448354,]

result<-result[!(result$SNP %in% fil$SNP),]

manhattan<-data.frame(sid=c(result$SNP),
                      chr=c(result$CHR),
                      pos=c(result$BP),
                      pvalue=c(result$P))
manhattan$highlight<-c(0)
manhattan$highlight_shape<-c(0)
manhattan$label<-c("")

for (i in unique(result$gene)){
  gsn<-result[result$gene==i]
  p<-gsn$P[which.min(gsn$P)]
  sid<-gsn$SNP[which.min(gsn$P)]
  
  if(p < 1e-06) {
    
    manhattan$highlight[which(manhattan$sid==sid)] <-1
    manhattan$highlight_shape[which(manhattan$sid==sid)] <-5
    manhattan$label[which(manhattan$sid==sid)] <-i
    
    } else { gsn < NULL }
  
}

# manhattan$label[which.min(manhattan$pvalue)]<-'rs4445544:NSUN6'
manhattan<-manhattan[,-1]
manhattanp<-manhattan[manhattan$pvalue<0.8,]

p<-list()

p[[1]]<-
  fig_manhattan(
    data = manhattanp,
    thresh = c(1e-06, 5e-08),
  )  +theme(aspect.ratio=0.5)

p[[2]]<-
  fig_qq(
    pvalues = manhattan$pvalue
  ) +theme(aspect.ratio=1)


i<-list[15]
path<-paste0('/scratch/lb4489/project/GWAS/',i,'ma&qq.pdf')

pdf(path,width = 15,height = 6)
print(wrap_plots(p,col=2) )
dev.off()


#BIPI

i<-list[6]

result<-fread(paste0('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/',i))

snps <- result[,c(2,1,3)]
colnames(snps) <- c("chr","SNP","Position")
snps$chr <- sub("^", "chr", snps$chr)

snps <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$SNP)

olap <- findOverlaps(snps,gene_span)
snpgene <- snps[queryHits(olap)]

result<-result[match( snpgene$rsid,result$SNP),]
result$gene<-gene_span$gene[subjectHits(olap)]
match<-geneAnno1[match(result$gene,geneAnno1$ensembl_gene_id),]
result$gene<-match$hgnc_symbol

fil<-result[result$CHR==6,]
fil<-fil[fil$BP > 28477797 & fil$BP< 33448354,]

result<-result[!(result$SNP %in% fil$SNP),]

manhattan<-data.frame(sid=c(result$SNP),
                      chr=c(result$CHR),
                      pos=c(result$BP),
                      pvalue=c(result$P))
manhattan$highlight<-c(0)
manhattan$highlight_shape<-c(0)
manhattan$label<-c("")

for (i in unique(result$gene)){
  gsn<-result[result$gene==i]
  p<-gsn$P[which.min(gsn$P)]
  sid<-gsn$SNP[which.min(gsn$P)]
  
  if(p < 1e-06) {
    
    manhattan$highlight[which(manhattan$sid==sid)] <-1
    manhattan$highlight_shape[which(manhattan$sid==sid)] <-5
    manhattan$label[which(manhattan$sid==sid)] <-i
    
  } else { gsn < NULL }
  
}

# manhattan$label[which.min(manhattan$pvalue)]<-'rs4445544:NSUN6'
manhattan<-manhattan[,-1]
manhattanp<-manhattan[manhattan$pvalue<0.8,]
manhattanp$chr<-as.character(manhattanp$chr)

p<-list()

p[[1]]<-
  fig_manhattan(
    data = manhattanp,
    thresh = c(1e-06, 5e-08)
  )  +theme(aspect.ratio=0.5)

p[[2]]<-
  fig_qq(
    pvalues = manhattan$pvalue
  ) +theme(aspect.ratio=1)


i<-list[6]
path<-paste0('/scratch/lb4489/project/GWAS/',i,'ma&qq.pdf')

pdf(path,width = 15,height = 6)
print(wrap_plots(p,col=2) )
dev.off()


#####ALCDEP

i<-list[3]

result<-fread(paste0('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/',i))

snps <- result[,c(2,1,3)]
colnames(snps) <- c("chr","SNP","Position")
snps$chr <- sub("^", "chr", snps$chr)

snps <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$SNP)

olap <- findOverlaps(snps,gene_span)
snpgene <- snps[queryHits(olap)]

result<-result[result$SNP %in% snpgene$rsid,]
result<-result[grep('rs',result$SNP),]
fil<-result[result$CHR==6,]
fil<-fil[fil$BP > 28477797 & fil$BP< 33448354,]
result<-result[!(result$SNP %in% fil$SNP),]

manhattan<-data.frame(chr=c(result$CHR),
                      pos=c(result$BP),
                      pvalue=c(result$P))
manhattan$label<-c("")
manhattanp<-manhattan[manhattan$pvalue<0.8,]

p<-list()

p[[1]]<-
  fig_manhattan(
    data = manhattan,
    thresh = c(1e-06, 5e-08)
  )  +theme(aspect.ratio=0.5)

p[[2]]<-
  fig_qq(
    pvalues = manhattan$pvalue
  ) +theme(aspect.ratio=1)



path<-paste0('/scratch/lb4489/project/GWAS/',i,'ma&qq.pdf')

pdf(path,width = 15,height = 6)
print(wrap_plots(p,col=2) )
dev.off()




