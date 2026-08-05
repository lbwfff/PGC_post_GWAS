
#对colocdb结果的分析，此处为最新版本

library(ggplot2)
library(patchwork)
library(MetBrewer)
library(bbplot)

list<-list.files('property/coloc_tier12/')

list<-list[grep('COLOCdb',list)]

smr<-list[grep('SMR',list)]
coloc<-list[-grep('SMR',list)]

cols <- c(
  "Schizophrenia" = met.brewer("Tsimshian", 6)[1],
  "Bipolar disorder" = met.brewer("Tsimshian", 6)[2],
  "Other"          = "#A0A0A0")

smr_colocdb<-data.frame()
coloc_colocdb<-data.frame()

for (i in smr){
  
  gene<-unlist(strsplit(i, '_'))[3]
  label<-unlist(strsplit(gene, '[.]'))[1]
  
  smrresult<-read.csv(paste0('property/coloc_tier12/',i))
  colocresult<-read.csv(paste0('property/coloc_tier12/COLOCdb_',gene))
  
  smrresult$molecule<-factor(smrresult$molecule)
  smrresult$Trait<-ifelse(grepl('schizop',smrresult$trait_description,ignore.case = TRUE),'Schizophrenia',
                          ifelse(grepl('Bipolar',smrresult$trait_description,ignore.case = TRUE),'Bipolar disorder','Other'))
  smrresult$gene<-label
  smr_colocdb<-rbind(smr_colocdb,smrresult)

  
  colocresult$molecule<-factor(colocresult$molecule)
  colocresult$Trait<-ifelse(grepl('schizop',colocresult$trait_description,ignore.case = TRUE),'Schizophrenia',
                            ifelse(grepl('Bipolar',colocresult$trait_description,ignore.case = TRUE),'Bipolar disorder','Other'))
  
  if(nrow(colocresult)==0) { 
    print(paste0('PASS gene : ',label))
    
  } else{
    
    colocresult$gene<-label
    coloc_colocdb<-rbind(coloc_colocdb,colocresult)
    
    
  }
  
}

# smr_colocdb<-smr_colocdb[smr_colocdb$gene %in% table$RMP,]

smr_colocdb<-smr_colocdb[smr_colocdb$p_heidi>0.05,]

#############################

library(dplyr)

BIPstats_result <- smr_colocdb %>%
  group_by(gene) %>%
  summarise(
    Total_Significant_Traits = n(),
    total_z2 = sum(qnorm(p_smr / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    target_z2 = sum(qnorm(p_smr[Trait =='Bipolar disorder'] / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    Signal_Proportion = ifelse(total_z2 > 0, target_z2 / total_z2, 0),
    mean_target_z2 = mean(qnorm(p_smr[Trait =='Bipolar disorder'] / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    mean_other_z2 = mean(qnorm(p_smr[Trait =='Other'] / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    Enrichment_Fold = mean_target_z2 / mean_other_z2
  ) %>%
  mutate(across(everything(), ~replace(., is.nan(.), 0)))

print(stats_result)

library(ggplot2)
library(ggrepel)
library(viridis)

p<-list()

list<-c('NSUN2','QTRT1','MRM1','NSUN6','TRMT61A')
list2<-c("CACNA1C", "ANK3", "TRANK1","SYNE1","TENM4") 

BIPstats_result<-BIPstats_result[BIPstats_result$gene %in% c(list,list2),]
BIPstats_result$group<-ifelse(BIPstats_result$gene %in% list,'Tier1','Other')

p[[1]]<-
ggplot(BIPstats_result, aes(x = log10(Total_Significant_Traits), y = Signal_Proportion,shape=group)) +
  geom_point(data =BIPstats_result[BIPstats_result$gene %in% list,],
             aes(color = Enrichment_Fold,size = Enrichment_Fold),alpha = 0.6) + #size = Enrichment_Fold, 
  geom_point(data =BIPstats_result[!(BIPstats_result$gene %in% list),],
             aes(color = Enrichment_Fold,size = Enrichment_Fold, ) ,alpha = 0.6) +
  scale_color_viridis_c(option = "plasma", name = "Enrichment Fold",limits = c(0.5,1.25),oob = scales::squish) +
  # scale_size_continuous(range = c(2, 12), name = "Enrichment Fold") +
  geom_text_repel(data = head(BIPstats_result %>% arrange(desc(Signal_Proportion)), 8),
                  aes(label = gene), size = 3) +
  labs(x = "Number of Significant Traits(log10)",
       y = "BIP Signal Proportion") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")+
  theme(aspect.ratio=1)


SCZstats_result <- smr_colocdb %>%
  group_by(gene) %>%
  summarise(
    Total_Significant_Traits = n(),
    total_z2 = sum(qnorm(p_smr / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    target_z2 = sum(qnorm(p_smr[Trait =='Schizophrenia'] / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    Signal_Proportion = ifelse(total_z2 > 0, target_z2 / total_z2, 0),
    mean_target_z2 = mean(qnorm(p_smr[Trait =='Schizophrenia'] / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    mean_other_z2 = mean(qnorm(p_smr[Trait =='Other'] / 2, lower.tail = FALSE)^2, na.rm = TRUE),
    Enrichment_Fold = mean_target_z2 / mean_other_z2
  ) %>%
  mutate(across(everything(), ~replace(., is.nan(.), 0)))

list<-c('NSUN2','TYW5','NSUN6','TRMT61A','MRM1')
list2<-c("DRD2", "GRIN2A", "CACNA1C")

SCZstats_result<-SCZstats_result[SCZstats_result$gene %in% c(list,list2),]
SCZstats_result$group<-ifelse(SCZstats_result$gene %in% list,'Tier1','Other')


p[[2]]<-
ggplot(SCZstats_result, aes(x = log10(Total_Significant_Traits), y = Signal_Proportion,shape=group)) +
  geom_point(data =SCZstats_result[SCZstats_result$gene %in% list,],
             aes(color = Enrichment_Fold,size = Enrichment_Fold),alpha = 0.6) + #size = Enrichment_Fold, 
  geom_point(data =SCZstats_result[!(SCZstats_result$gene %in% list),],
             aes(color = Enrichment_Fold,size = Enrichment_Fold, ) ,alpha = 0.6) +
  scale_color_viridis_c(option = "plasma", name = "Enrichment Fold",limits = c(0.5,1.25),oob = scales::squish) +
  # scale_size_continuous(range = c(2, 12), name = "Enrichment Fold") +
  geom_text_repel(data = head(SCZstats_result %>% arrange(desc(Signal_Proportion)), 10),
                  aes(label = gene), size = 3) +
  labs(x = "Number of Significant Traits(log10)",
       y = "SCZ Signal Proportion") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")+
  coord_cartesian(ylim=c(0,0.48))+
  theme(aspect.ratio=1)


pdf("./property/SMR_meta_n.pdf",width = 10,height = 5)
patchwork::wrap_plots(p,nrow=1) 
dev.off()

######
#把coloc也整理出来吧，作为补充图

library(dplyr)

BIPstats_result <- coloc_colocdb %>%
  group_by(gene) %>%
  summarise(
    Total_Coloc_Traits = n(),
    total_pp_sum = sum(pp_h4_abf, na.rm = TRUE),
    target_pp_sum = sum(pp_h4_abf[Trait =='Bipolar disorder'], na.rm = TRUE),
    Coloc_Proportion = ifelse(total_pp_sum > 0, target_pp_sum / total_pp_sum, 0),
    mean_target_pp = mean(pp_h4_abf[Trait =='Bipolar disorder'], na.rm = TRUE),
    mean_other_pp = mean(pp_h4_abf[Trait == 'Other'], na.rm = TRUE),
    Coloc_Enrichment = ifelse(mean_other_pp > 0, mean_target_pp / mean_other_pp, Inf)
  ) %>%
  mutate(across(everything(), ~replace(., is.nan(.), 0)))

p<-list()

list<-c('NSUN2','QTRT1','MRM1','NSUN6','TRMT61A')
list2<-c("CACNA1C", "ANK3", "TRANK1","SYNE1","TENM4") 

BIPstats_result<-BIPstats_result[BIPstats_result$gene %in% c(list,list2),]
BIPstats_result$group<-ifelse(BIPstats_result$gene %in% list,'Tier1','Other')

p[[1]]<-

ggplot(BIPstats_result, aes(x = Total_Coloc_Traits, y = Coloc_Proportion)) +
  geom_point(data =BIPstats_result[BIPstats_result$gene %in% list,],
             aes(color = Coloc_Enrichment,size = Coloc_Enrichment),alpha = 0.6) + #size = Enrichment_Fold, 
  geom_point(data =BIPstats_result[!(BIPstats_result$gene %in% list),],
             aes(color = Coloc_Enrichment,size = Coloc_Enrichment, ) ,alpha = 0.6) +
  scale_color_viridis_c(option = "plasma", name = "Enrichment Fold",limits = c(0.5,1.25),oob = scales::squish) +
  geom_text_repel(data = head(BIPstats_result %>% arrange(desc(Coloc_Proportion)), 20),
                  aes(label = gene), size = 3) +
  labs(x = "Number of Significant Traits",
       y = "BIP Signal Proportion") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")+
  theme(aspect.ratio=1)


SCZstats_result <- coloc_colocdb %>%
  group_by(gene) %>%
  summarise(
    Total_Coloc_Traits = n(),
    total_pp_sum = sum(pp_h4_abf, na.rm = TRUE),
    target_pp_sum = sum(pp_h4_abf[Trait =='Schizophrenia'], na.rm = TRUE),
    Coloc_Proportion = ifelse(total_pp_sum > 0, target_pp_sum / total_pp_sum, 0),
    mean_target_pp = mean(pp_h4_abf[Trait =='Schizophrenia'], na.rm = TRUE),
    mean_other_pp = mean(pp_h4_abf[Trait == 'Other'], na.rm = TRUE),
    Coloc_Enrichment = ifelse(mean_other_pp > 0, mean_target_pp / mean_other_pp, Inf)
  ) %>%
  mutate(across(everything(), ~replace(., is.nan(.), 0)))

list<-c('NSUN2','TYW5','NSUN6','TRMT61A','MRM1')
list2<-c("DRD2", "GRIN2A", "CACNA1C")

SCZstats_result<-SCZstats_result[SCZstats_result$gene %in% c(list,list2),]
SCZstats_result$group<-ifelse(SCZstats_result$gene %in% list,'Tier1','Other')

p[[2]]<-
  
  ggplot(SCZstats_result, aes(x = Total_Coloc_Traits, y = Coloc_Proportion)) +
  geom_point(data =SCZstats_result[SCZstats_result$gene %in% list,],
             aes(color = Coloc_Enrichment,size = Coloc_Enrichment),alpha = 0.6) + #size = Enrichment_Fold, 
  geom_point(data =SCZstats_result[!(SCZstats_result$gene %in% list),],
             aes(color = Coloc_Enrichment,size = Coloc_Enrichment, ) ,alpha = 0.6) +
  scale_color_viridis_c(option = "plasma", name = "Enrichment Fold",limits = c(0.5,1.25),oob = scales::squish) +
  geom_text_repel(data = head(SCZstats_result %>% arrange(desc(Coloc_Proportion)), 20),
                  aes(label = gene), size = 3) +
  labs(x = "Number of Significant Traits",
       y = "BIP Signal Proportion") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")+
  theme(aspect.ratio=1)

pdf("./property/coloc_meta.pdf", width = 6, height = 5)
patchwork::wrap_plots(p,nrow=1) 
dev.off()


#########################
##画一张单独的SMR，再把coloc补上就行了


library(ggplot2)
library(patchwork)
library(MetBrewer)
library(bbplot)

list<-list.files('colocdb/')

list<-list[grep('COLOCdb',list)]

smr<-list[grep('SMR',list)]
coloc<-list[-grep('SMR',list)]

p1<-list()
p2<-list()

cols <- c(
  "Schizophrenia" = "#7570b3",
  "Bipolar disorder" = "#d95f02",
  "Other"          = "#A0A0A0")


smr_colocdb<-data.frame()
coloc_colocdb<-data.frame()

for (i in smr){
  
  gene<-unlist(strsplit(i, '_'))[3]
  label<-unlist(strsplit(gene, '[.]'))[1]
  
  smrresult<-read.csv(paste0('colocdb/',i))
  smrresult<-smrresult[smrresult$p_heidi>0.05,]
  colocresult<-read.csv(paste0('colocdb/COLOCdb_',gene))
  
  smrresult$qtl_type<-factor(smrresult$qtl_type)
  smrresult$Trait<-ifelse(grepl('schizop',smrresult$trait_description,ignore.case = TRUE),'Schizophrenia',
                          ifelse(grepl('Bipolar',smrresult$trait_description,ignore.case = TRUE),'Bipolar disorder',
                                 'Other'))
  smrresult$gene<-label
  smrresult<-smrresult[!is.na(smrresult$qtl_type),]
  smr_colocdb<-rbind(smr_colocdb,smrresult)
  
  p1[[which(smr==i)]]<-
    ggplot(smrresult, aes(x = -log10(p_smr), y = qtl_type, color  = Trait)) +
    geom_jitter(data = subset(smrresult, Trait == "Other"),
                aes(x = -log10(p_smr), y = qtl_type, color = Trait),
                height = 0.2, size = 2, alpha = 0.6) +
    geom_jitter(data = subset(smrresult, Trait != "Other"),
                aes(x = -log10(p_smr), y = qtl_type, color = Trait),
                height = 0.2, size = 2) +
    theme_bw() +
    scale_color_manual(values = cols)+
    theme(axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank())+
    # bbc_style()+
    labs(y = "xQTL", x= "-log10(P SMR)",title=paste0(label,' SMR')) +
    theme(aspect.ratio=1)
  
  colocresult$molecule<-factor(colocresult$molecule)
  colocresult$Trait<-ifelse(grepl('schizop',colocresult$trait_description,ignore.case = TRUE),'Schizophrenia',
                            ifelse(grepl('Bipolar',colocresult$trait_description,ignore.case = TRUE),'Bipolar disorder','Other'))
  
  if(nrow(colocresult)==0) { 
    print(paste0('PASS gene : ',label))
    
  } else{
    
    colocresult$gene<-label
    coloc_colocdb<-rbind(coloc_colocdb,colocresult)
    
    p2[[which(smr==i)]]<-
      ggplot(colocresult, aes(y = molecule, x = pp_h4_abf, color  = Trait)) +
      geom_jitter(data = subset(colocresult, Trait == "Other"),
                  aes(x = pp_h4_abf, y = molecule, color = Trait),
                  height = 0.2, size = 2, alpha = 0.6) +
      geom_jitter(data = subset(colocresult, Trait != "Other"),
                  aes(x = pp_h4_abf, y = molecule, color = Trait),
                  height = 0.2, size = 2) +
      labs(y = "xQTL", x = "pp h4 abf",title=paste0(label,' Coloc')) +
      theme_bw() +
      scale_color_manual(values=cols)+
      theme(axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y  = element_blank())+
      theme(aspect.ratio=1)
    
  }
  
}


pdf("./property/colocdb_SMR.pdf",width = 12,height = 8)
wrap_plots(p1,nrow=2) 
dev.off()

clean_plot_list <- lapply(p2, function(p) {
  tryCatch({
    ggplot_build(p) 
    return(p) 
  }, error = function(e) {
    return(ggplot() + 
             theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "Data Missing", color = "red"))
  })
})

pdf("./property/colocdb_coloc.pdf", width = 12, height = 8)
wrap_plots(clean_plot_list, nrow = 2)
dev.off()

#加一个eQTL only的，带方向的散点图，和转录组数据放在一起

sczeqtl<-smr_colocdb[smr_colocdb$qtl_type=='eQTL' & smr_colocdb$Trait=='Schizophrenia',]
bipeqtl<-smr_colocdb[smr_colocdb$qtl_type=='eQTL' & smr_colocdb$Trait=='Bipolar disorder',]

library(ggplot2)

slopes_df <- sczeqtl %>%
  group_by(gene) %>%
  summarise(slope = coef(lm(b_smr ~ (-log10(p_smr)) + 0))[1])

sczeqtl_with_slopes <- left_join(sczeqtl, slopes_df, by = "gene")

ggplot(sczeqtl, aes(x = -log10(p_smr), y = b_smr, color = gene)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x - 1, se = FALSE, linewidth = 1)+
  theme_minimal(base_size = 12) +
  labs(
    x = "p_smr",
    y = "b_smr",
    title = "SMR Results by Gene"
  ) +
  theme(aspect.ratio = 1)


#################



library(dplyr)

smrbip <- smr_colocdb [grep('Bipolar disorder',smr_colocdb$trait_description),]
list<-c('NSUN2','TYW5','NSUN6','TRMT61A','MRM1',"QTRT1")
smrbip<-smrbip[smrbip$gene %in% list,]
smrbip<-smrbip[smrbip$qtl_type=='eQTL',]


smrscz <- smr_colocdb [grep('Schizophrenia',smr_colocdb$trait_description),]
list<-c('NSUN2','TYW5','NSUN6','TRMT61A','MRM1',"QTRT1")
smrscz<-smrscz[smrscz$gene %in% list,]
smrscz<-smrscz[smrscz$qtl_type=='eQTL',]

