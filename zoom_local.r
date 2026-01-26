library(ggplot2)

library(patchwork)

# pair<-read.csv('/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/network_Mag&hmag.csv')
# pair<-pair[!duplicated(paste0(pair$V6,'_',pair$disease)),]

pair<-data.frame(V6=c('TYW5','DCP1A','NSUN6','TRMT61A',
                      'NSUN2','MRM1','QTRT1'),
                 disease=c('SCZ','SCZ','SCZ','SCZ',
                           'BIPI','BIPI','BIPI'))

load('/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/geneAnno_allgenes.rda')
exon <- read.table("/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/Gencode26_exon.bed")
exon$V1 <- sub("^", "chr", exon$V1)
promoter <- read.table("/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol/Gencode26_promoter.bed")
promoter$V1 <- sub("^", "chr", promoter$V1)

exonranges <- GRanges(exon[,1],IRanges(exon[,2],exon[,3]),gene=exon[,4])
promoterranges <- GRanges(promoter[,1], IRanges(promoter[,2], promoter[,3]), gene=promoter[,4])

polyfun<-c()

for (i in 1:nrow(pair)) {
  
  pgcd<-pair$disease[i]
  name<-pair$V6[i]
  gene<-geneAnno1$ensembl_gene_id[geneAnno1$hgnc_symbol==name]
  
  setwd("/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol")
  
  allsnps <- data.table::fread(paste0('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/PGC_',pgcd,'_for_magma'))
  
  nsunexon<-exonranges [exonranges$gene==gene]
  nsunpromote<-promoterranges [promoterranges$gene==gene]
  genebody <- range(c(nsunexon,nsunpromote))

  start<-genebody@ranges@start-250000
  end<-genebody@ranges@start+genebody@ranges@width+250000

  plot<-allsnps[paste0('chr',allsnps$CHR) == as.character(genebody@seqnames),]
  plot<-plot[plot$BP> start & plot$BP< end,]
  
  assoc<-plot[,c('SNP','CHR','BP','P')]
  colnames(assoc)<-c('marker','chr','pos','pvalue')
  
  corr<-ieugwasr::ld_matrix(assoc$marker, with_alleles = F, 
                            plink_bin = plinkbinr::get_plink_exe(),
                            bfile = '/scratch/lb4489/project/GWAS/MAGMA/g1000_eur')
  
  assoc<-assoc[match(rownames(corr),assoc$marker),]
  assoc$pos<-as.integer(assoc$pos)
  
  library(geni.plots)

  p<-
  fig_region(
      data = assoc,
      corr = corr,
      title = c(paste0(pgcd,'-',name)),
      build = 37)
  
  p[[1]]<-
  p[[1]]+annotate(
    geom = "rect",
    xmin = genebody@ranges@start,
    xmax = genebody@ranges@start+genebody@ranges@width,
    ymin = -Inf, 
    ymax = Inf,  
    fill = "grey",
    alpha = 0.3    
  )
  
  pdf(paste0(pgcd,'_',name,'.pdf'),width = 8,height = 6)
  print(p)
  dev.off()
  
#fine mapping region
  
  get_ld_block_max_overlap <- function(chr, start, end, block_size = 3e6) {
    block_start <- floor(start / block_size) * block_size + 1
    block_end <- floor(end / block_size) * block_size + 1
    blocks <- seq(block_start, block_end, by = block_size)
    
    block_ranges <- data.frame(
      start = blocks,
      end = blocks + block_size
    )
    
    overlap <- pmin(end, block_ranges$end) - pmax(start, block_ranges$start) + 1
    overlap[overlap < 0] <- 0
    
    max_idx <- which.max(overlap)
    block_start_max <- block_ranges$start[max_idx]
    block_end_max <- block_ranges$end[max_idx]
    
    ld_file <- sprintf("%s_%d_%d", chr, block_start_max, block_end_max)
    return(ld_file)
  }
  
  fileld<-get_ld_block_max_overlap(paste0('chr',plot$CHR[1]), start, end)
  
  start=max(start,unlist(strsplit(fileld,'_'))[2])
  stop=min(end,unlist(strsplit(fileld,'_'))[3])
  
  file=c(paste0(fileld,".gz"))
  file2=c(paste0(fileld,".npz"))
  
  cache<-paste0('python ./polyfun/finemapper.py --allow-missing --ld ./LD_blocks/',fileld,' --sumstats ./output/',pgcd,'_with_var.gz --n ',max(plot$NOBS),' --chr ',plot$CHR[1],' --start ',start ," --end ", end,' --method susie --max-num-causal 2 --out output/',name,'_',pgcd,'_UKB.gz')
  
  polyfun<-c(polyfun,cache)
  
}

polyfun<-unique(polyfun)
write.table(polyfun,file = '/scratch/lb4489/project/GWAS/PTSD_3/polyfun_PGC_RMRs.sh',
            quote = F,sep = '\t',row.names = F,col.names = F)


############
#把fine mapping放进来

pip<-data.frame()

for (i in 1:nrow(pair)) {
  
  pgcd<-pair$disease[i]
  name<-pair$V6[i]
  gene<-geneAnno1$ensembl_gene_id[geneAnno1$hgnc_symbol==name]
  
  setwd("/scratch/lb4489/project/GWAS/hMAGMA/HMAGMA_Protocol")
  
  allsnps <- data.table::fread(paste0('/scratch/lb4489/project/GWAS/GWAS_from_PGC/MAGMA_GWAS_input/PGC_',pgcd,'_for_magma'))
  
  nsunexon<-exonranges [exonranges$gene==gene]
  nsunpromote<-promoterranges [promoterranges$gene==gene]
  genebody <- range(c(nsunexon,nsunpromote))
  
  start<-genebody@ranges@start-250000
  end<-genebody@ranges@start+genebody@ranges@width+250000
  
  plot<-allsnps[paste0('chr',allsnps$CHR) == as.character(genebody@seqnames),]
  plot<-plot[plot$BP> start & plot$BP< end,]
  
  assoc<-plot[,c('SNP','CHR','BP','P')]
  colnames(assoc)<-c('marker','chr','pos','pvalue')
  
  finemapresult<-paste0('/scratch/lb4489/project/GWAS/PTSD_3/output/',name,'_',pgcd,'_UKB.gz')
  finemapresult<-fread(finemapresult)
  
  assoc2<-data.frame(marker=c(finemapresult$SNP),
                    chr=c(finemapresult$CHR),
                    pos=c(finemapresult$BP),
                    prob=c(as.numeric(finemapresult$PIP)))
  
  assoc<-assoc[assoc$marker %in% assoc2$marker,]
  assoc2<-assoc2[assoc2$marker %in% assoc$marker,]
  
  corr<-ieugwasr::ld_matrix(assoc$marker, with_alleles = F, 
                            plink_bin = plinkbinr::get_plink_exe(),
                            bfile = '/scratch/lb4489/project/GWAS/MAGMA/g1000_eur')
  
  assoc<-assoc[match(rownames(corr),assoc$marker),]
  assoc$pos<-as.integer(assoc$pos)
  assoc2<-assoc2[match(rownames(corr),assoc2$marker),]
  assoc2$pos<-as.integer(assoc2$pos)
  
  library(geni.plots)
  
  p<-list()
  
  p[[1]]<-
    fig_region(
      data = assoc,
      corr = corr,
      title = c(paste0(pgcd,'-',name)),
      build = 37,title_center = TRUE)
  
  p[[1]][[1]]<-
    p[[1]][[1]]+annotate(
      geom = "rect",
      xmin = genebody@ranges@start,
      xmax = genebody@ranges@start+genebody@ranges@width,
      ymin = -Inf, 
      ymax = Inf,  
      fill = "grey",
      alpha = 0.3    
    )
  
  p[[2]]<-
    fig_region(
      data = assoc2,
      corr = corr,
      build = 37,
      prob=T,
      title_center = TRUE)
  
  p[[2]][[1]]<-
    p[[2]][[1]]+annotate(
      geom = "rect",
      xmin = genebody@ranges@start,
      xmax = genebody@ranges@start+genebody@ranges@width,
      ymin = -Inf, 
      ymax = Inf,  
      fill = "grey",
      alpha = 0.3    
    )
  
  remove_legend <- function(p) {
    p + theme(legend.position = "none")
  }
  
  final_plot <- wrap_plots(
    remove_legend(p[[1]][[1]]),
    remove_legend(p[[2]][[1]]),
    p[[2]][[2]],
    ncol = 1,
    heights = c(1, 1, 0.5,0.75) 
  )
  
  # final_plot
  
  pdf(paste0(pgcd,'_',name,'_finemapping.pdf'),width = 8,height = 6)
  print(final_plot)
  dev.off()
  
  max<-finemapresult[finemapresult$BP > genebody@ranges@start & finemapresult$BP< (genebody@ranges@start+genebody@ranges@width) , ]
  maxpip<-data.frame(gene=c(name),
                     pgcd=c(pgcd),
                     value=c(max(max$PIP)))
  pip<-rbind(pip,maxpip)
  
  print(paste0(pgcd,'_',name,'_finemapping'))
  
  gc()
}

write.csv(pip,file = '/scratch/lb4489/project/GWAS/GIFT/maxpip.csv')


