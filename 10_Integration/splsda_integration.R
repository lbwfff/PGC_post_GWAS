
#对应splsda的结果，但是我还有另一个splsda的结果，最后用的哪个？
#两个结果应该内容是一样的，但是可视化有区别

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel) 
library(pcaMethods)
library(mixOmics)


p<-list()
cache<-data.frame()

##SCZ
positive_controls <- c("DRD2", "GRIN2A", "CACNA1C")  #"HTR2A","SYN2"
negative_controls <- c("TYR", "FABP4", "ALB", "TNNT3", "OR5M3")
Candidate<-unique(scoretable$RMP[scoretable$total.score>1])


filesave_wide <- filesave[filesave$Disorder %in% c('SCZ'),] %>% #
  pivot_wider(
    id_cols = Gene,
    names_from = c(Disorder, Method),
    names_sep = "_",
    values_from = Score)

filesave_wide<-filesave_wide[filesave_wide$Gene %in% c(positive_controls, negative_controls, Candidate),]
filesave_wide <- filesave_wide %>%
  select_if(~ !is.numeric(.) || (sum(!is.na(.)) > 4 && sum(. != 0, na.rm = TRUE) > 0))

pca_matrix<-as.data.frame(filesave_wide)
sczh2<-read.csv('sczh2_result.csv')
sczh2<-sczh2[match(pca_matrix$Gene,sczh2$GENE),]
sczh2$ratio<-sczh2$H2/sczh2$length
pca_matrix$sch2<-scale(sczh2$pve)
pca_matrix$sch2ratio<-scale(sczh2$ratio)

rownames(pca_matrix)<-pca_matrix$Gene
pca_matrix<-pca_matrix[,-1]
pca_matrix[is.na(pca_matrix)]<-0
pca_matrix<-abs(pca_matrix)

ppca_res <- pcaMethods::pca(pca_matrix,
                            method = "ppca", nPcs = 2,
                            scale = c("pareto"), center = TRUE,)

completed_data <- as.data.frame(ppca_res@completeObs) %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(gene_annotation, by = "Gene")

control_genes <- c(positive_controls, negative_controls) 

X_train <- completed_data[completed_data$Gene %in% control_genes, ]
rownames(X_train) <- X_train$Gene

X_train <- as.matrix(X_train %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

train_labels <- ifelse(rownames(X_train) %in% positive_controls, "Positive", "Negative")
Y_train <- as.factor(train_labels)

X_test <- completed_data[!completed_data$Gene %in% control_genes, ]
rownames(X_test) <- X_test$Gene

X_test <- as.matrix(X_test %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

splsda_model <- splsda(X_train, Y_train, ncomp = 2, keepX = c(5, 5))

splsda_pred <- predict(splsda_model, newdata = X_test)

train_coords <- as.data.frame(splsda_model$variates$X)
train_coords$Gene <- rownames(train_coords)
train_coords$Group <- train_labels

test_coords <- as.data.frame(splsda_pred$variates)
test_coords$Gene <- rownames(X_test)
test_coords$Group <- "Candidate Gene"

colnames(test_coords)<-c("comp1", "comp2", "Gene", "Group")
plsda_plot_df <- rbind(train_coords, test_coords)

expl_vars <- round(splsda_model$prop_expl_var$X * 100, 2)
comp1_label <- paste0("sPLS-DA Component 1 (", expl_vars[1], "%)")
comp2_label <- paste0("sPLS-DA Component 2 (", expl_vars[2], "%)")

#中心点计算
centroids <- plsda_plot_df %>%
  filter(Group %in% c("Positive", "Negative")) %>%
  group_by(Group) %>%
  summarize(Comp1 = mean(comp1), Comp2 = mean(comp2))
pos_center <- c(centroids$Comp1[centroids$Group == "Positive"], centroids$Comp2[centroids$Group == "Positive"])
neg_center <- c(centroids$Comp1[centroids$Group == "Negative"], centroids$Comp2[centroids$Group == "Negative"])

p[[1]]<-
  ggplot() +
  geom_point(data = filter(plsda_plot_df, Group == "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = filter(plsda_plot_df, Group != "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = centroids, 
             aes(x = Comp1, y = Comp2, fill = Group), 
             shape = 23, size = 5, color = NULL, stroke = 0) +
  scale_color_manual(values = c("Candidate Gene" = "#8c8c8c", 
                                "Positive" = "#d95f02", 
                                "Negative" = "#7570b3")) +
  scale_fill_manual(values = c("Positive" = "#d95f02", 
                               "Negative" = "#7570b3"), guide = "none") +
  stat_ellipse(data = filter(plsda_plot_df, Group %in% c("Positive", "Negative")),
               aes(x = comp1, y = comp2, fill = Group), 
               geom = "polygon", 
               alpha = 0, level = 0.95, 
               show.legend = FALSE) +
  geom_text_repel(data = plsda_plot_df, 
                  aes(x = comp1, y = comp2, label = Gene), 
                  size = 3.8, max.overlaps = 5, 
                  box.padding = 0.5,
                  point.padding = 0.3,
                  show.legend = FALSE) +
  labs(x = comp1_label, y = comp2_label, 
       title = "SCZ") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )+
  theme(aspect.ratio=1)

#距离分数

W1 <- splsda_model$prop_expl_var$X[1] 
W2 <- splsda_model$prop_expl_var$X[2]  
candidates <- plsda_plot_df %>% filter(Group == "Candidate Gene")
candidate_scores <- candidates %>%
  mutate(
    D_Pos = sqrt( W1 * (comp1 - pos_center[1])^2 + W2 * (comp2 - pos_center[2])^2 ),
    D_Neg = sqrt( W1 * (comp1 - neg_center[1])^2 + W2 * (comp2 - neg_center[2])^2 ),
    Score = D_Neg / (D_Pos + D_Neg)
  ) %>%
  dplyr::select(Gene, comp1, comp2, D_Pos, D_Neg, Score) %>%
  arrange(desc(Score))

candidate_scores$disorder<-c('SCZ')
cache<-rbind(cache,candidate_scores)


##


positive_controls <- c("CACNA1C", "ANK3", "TRANK1","SYNE1","TENM4") 
negative_controls <- c("TYR", "FABP4", "ALB", "TNNT3", "OR5M3")
Candidate<-unique(scoretable$RMP[scoretable$total.score>1])


filesave_wide <- filesave[filesave$Disorder %in% c('BIPI'),] %>% #
  pivot_wider(
    id_cols = Gene,
    names_from = c(Disorder, Method),
    names_sep = "_",
    values_from = Score)

filesave_wide<-filesave_wide[filesave_wide$Gene %in% c(positive_controls, negative_controls, Candidate),]
filesave_wide <- filesave_wide %>%
  select_if(~ !is.numeric(.) || (sum(!is.na(.)) > 4 && sum(. != 0, na.rm = TRUE) > 0))

pca_matrix<-as.data.frame(filesave_wide)
sczh2<-read.csv('BIPIh2_result.csv')
sczh2<-sczh2[match(pca_matrix$Gene,sczh2$GENE),]
sczh2$ratio<-sczh2$H2/sczh2$length
pca_matrix$sch2<-scale(sczh2$pve)
pca_matrix$sch2ratio<-scale(sczh2$ratio)

rownames(pca_matrix)<-pca_matrix$Gene
pca_matrix<-pca_matrix[,-1]
pca_matrix[is.na(pca_matrix)]<-0
pca_matrix<-abs(pca_matrix)

ppca_res <- pcaMethods::pca(pca_matrix,
                            method = "ppca", nPcs = 2,
                            scale = c("pareto"), center = TRUE,)

completed_data <- as.data.frame(ppca_res@completeObs) %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(gene_annotation, by = "Gene")

control_genes <- c(positive_controls, negative_controls) 

X_train <- completed_data[completed_data$Gene %in% control_genes, ]
rownames(X_train) <- X_train$Gene

X_train <- as.matrix(X_train %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

train_labels <- ifelse(rownames(X_train) %in% positive_controls, "Positive", "Negative")
Y_train <- as.factor(train_labels)

X_test <- completed_data[!completed_data$Gene %in% control_genes, ]
rownames(X_test) <- X_test$Gene

X_test <- as.matrix(X_test %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

splsda_model <- splsda(X_train, Y_train, ncomp = 2, keepX = c(5, 5))

splsda_pred <- predict(splsda_model, newdata = X_test)

train_coords <- as.data.frame(splsda_model$variates$X)
train_coords$Gene <- rownames(train_coords)
train_coords$Group <- train_labels

test_coords <- as.data.frame(splsda_pred$variates)
test_coords$Gene <- rownames(X_test)
test_coords$Group <- "Candidate Gene"

colnames(test_coords)<-c("comp1", "comp2", "Gene", "Group")
plsda_plot_df <- rbind(train_coords, test_coords)

expl_vars <- round(splsda_model$prop_expl_var$X * 100, 2)
comp1_label <- paste0("sPLS-DA Component 1 (", expl_vars[1], "%)")
comp2_label <- paste0("sPLS-DA Component 2 (", expl_vars[2], "%)")

#中心点计算
centroids <- plsda_plot_df %>%
  filter(Group %in% c("Positive", "Negative")) %>%
  group_by(Group) %>%
  summarize(Comp1 = mean(comp1), Comp2 = mean(comp2))
pos_center <- c(centroids$Comp1[centroids$Group == "Positive"], centroids$Comp2[centroids$Group == "Positive"])
neg_center <- c(centroids$Comp1[centroids$Group == "Negative"], centroids$Comp2[centroids$Group == "Negative"])

p[[2]]<-
  ggplot() +
  geom_point(data = filter(plsda_plot_df, Group == "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = filter(plsda_plot_df, Group != "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = centroids, 
             aes(x = Comp1, y = Comp2, fill = Group), 
             shape = 23, size = 5, color = NULL, stroke = 0) +
  scale_color_manual(values = c("Candidate Gene" = "#8c8c8c", 
                                "Positive" = "#d95f02", 
                                "Negative" = "#7570b3")) +
  scale_fill_manual(values = c("Positive" = "#d95f02", 
                               "Negative" = "#7570b3"), guide = "none") +
  stat_ellipse(data = filter(plsda_plot_df, Group %in% c("Positive", "Negative")),
               aes(x = comp1, y = comp2, fill = Group), 
               geom = "polygon", 
               alpha = 0, level = 0.95, 
               show.legend = FALSE) +
  geom_text_repel(data = plsda_plot_df, 
                  aes(x = comp1, y = comp2, label = Gene), 
                  size = 3.8, max.overlaps = 5, 
                  box.padding = 0.5,
                  point.padding = 0.3,
                  show.legend = FALSE) +
  labs(x = comp1_label, y = comp2_label, 
       title = "BIPI") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )+
  theme(aspect.ratio=1)

#距离分数

W1 <- splsda_model$prop_expl_var$X[1] 
W2 <- splsda_model$prop_expl_var$X[2]  
candidates <- plsda_plot_df %>% filter(Group == "Candidate Gene")
candidate_scores <- candidates %>%
  mutate(
    D_Pos = sqrt( W1 * (comp1 - pos_center[1])^2 + W2 * (comp2 - pos_center[2])^2 ),
    D_Neg = sqrt( W1 * (comp1 - neg_center[1])^2 + W2 * (comp2 - neg_center[2])^2 ),
    Score = D_Neg / (D_Pos + D_Neg)
  ) %>%
  dplyr::select(Gene, comp1, comp2, D_Pos, D_Neg, Score) %>%
  arrange(desc(Score))

candidate_scores$disorder<-c('BIPI')
cache<-rbind(cache,candidate_scores)


##AD
#AD这个结果是很符合直觉的

positive_controls <- c("APOE","BIN1","TREM2") 
negative_controls <- c("TYR", "FABP4", "ALB", "TNNT3", "OR5M3")
Candidate<-unique(scoretable$RMP[scoretable$total.score>1])


filesave_wide <- filesave[filesave$Disorder %in% c('AD'),] %>% #
  pivot_wider(
    id_cols = Gene,
    names_from = c(Disorder, Method),
    names_sep = "_",
    values_from = Score)

filesave_wide<-filesave_wide[filesave_wide$Gene %in% c(positive_controls, negative_controls, Candidate),]
filesave_wide <- filesave_wide %>%
  select_if(~ !is.numeric(.) || (sum(!is.na(.)) > 4 && sum(. != 0, na.rm = TRUE) > 0))

pca_matrix<-as.data.frame(filesave_wide)
sczh2<-read.csv('adh2_result.csv')
sczh2<-sczh2[match(pca_matrix$Gene,sczh2$GENE),]
sczh2$ratio<-sczh2$H2/sczh2$length
pca_matrix$sch2<-scale(sczh2$pve)
pca_matrix$sch2ratio<-scale(sczh2$ratio)

rownames(pca_matrix)<-pca_matrix$Gene
pca_matrix<-pca_matrix[,-1]
pca_matrix[is.na(pca_matrix)]<-0
pca_matrix[pca_matrix=='Inf']<-16
pca_matrix<-abs(pca_matrix)


ppca_res <- pcaMethods::pca(pca_matrix,
                            method = "ppca", nPcs = 2,
                            scale = c("pareto"), center = TRUE,)

completed_data <- as.data.frame(ppca_res@completeObs) %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(gene_annotation, by = "Gene")

control_genes <- c(positive_controls, negative_controls) 

X_train <- completed_data[completed_data$Gene %in% control_genes, ]
rownames(X_train) <- X_train$Gene

X_train <- as.matrix(X_train %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

train_labels <- ifelse(rownames(X_train) %in% positive_controls, "Positive", "Negative")
Y_train <- as.factor(train_labels)

X_test <- completed_data[!completed_data$Gene %in% control_genes, ]
rownames(X_test) <- X_test$Gene

X_test <- as.matrix(X_test %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

splsda_model <- splsda(X_train, Y_train, ncomp = 2, keepX = c(5, 5))

splsda_pred <- predict(splsda_model, newdata = X_test)

train_coords <- as.data.frame(splsda_model$variates$X)
train_coords$Gene <- rownames(train_coords)
train_coords$Group <- train_labels

test_coords <- as.data.frame(splsda_pred$variates)
test_coords$Gene <- rownames(X_test)
test_coords$Group <- "Candidate Gene"

colnames(test_coords)<-c("comp1", "comp2", "Gene", "Group")
plsda_plot_df <- rbind(train_coords, test_coords)

expl_vars <- round(splsda_model$prop_expl_var$X * 100, 2)
comp1_label <- paste0("sPLS-DA Component 1 (", expl_vars[1], "%)")
comp2_label <- paste0("sPLS-DA Component 2 (", expl_vars[2], "%)")

#中心点计算
centroids <- plsda_plot_df %>%
  filter(Group %in% c("Positive", "Negative")) %>%
  group_by(Group) %>%
  summarize(Comp1 = mean(comp1), Comp2 = mean(comp2))
pos_center <- c(centroids$Comp1[centroids$Group == "Positive"], centroids$Comp2[centroids$Group == "Positive"])
neg_center <- c(centroids$Comp1[centroids$Group == "Negative"], centroids$Comp2[centroids$Group == "Negative"])

p[[3]]<-
  ggplot() +
  geom_point(data = filter(plsda_plot_df, Group == "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = filter(plsda_plot_df, Group != "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = centroids, 
             aes(x = Comp1, y = Comp2, fill = Group), 
             shape = 23, size = 5, color = NULL, stroke = 0) +
  scale_color_manual(values = c("Candidate Gene" = "#8c8c8c", 
                                "Positive" = "#d95f02", 
                                "Negative" = "#7570b3")) +
  scale_fill_manual(values = c("Positive" = "#d95f02", 
                               "Negative" = "#7570b3"), guide = "none") +
  stat_ellipse(data = filter(plsda_plot_df, Group %in% c("Positive", "Negative")),
               aes(x = comp1, y = comp2, fill = Group), 
               geom = "polygon", 
               alpha = 0, level = 0.95, 
               show.legend = FALSE) +
  geom_text_repel(data = plsda_plot_df, 
                  aes(x = comp1, y = comp2, label = Gene), 
                  size = 3.8, max.overlaps = 5, 
                  box.padding = 0.5,
                  point.padding = 0.3,
                  show.legend = FALSE) +
  labs(x = comp1_label, y = comp2_label, 
       title = "AD") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )+
  theme(aspect.ratio=1)

#距离分数

W1 <- splsda_model$prop_expl_var$X[1] 
W2 <- splsda_model$prop_expl_var$X[2]  
candidates <- plsda_plot_df %>% filter(Group == "Candidate Gene")
candidate_scores <- candidates %>%
  mutate(
    D_Pos = sqrt( W1 * (comp1 - pos_center[1])^2 + W2 * (comp2 - pos_center[2])^2 ),
    D_Neg = sqrt( W1 * (comp1 - neg_center[1])^2 + W2 * (comp2 - neg_center[2])^2 ),
    Score = D_Neg / (D_Pos + D_Neg)
  ) %>%
  dplyr::select(Gene, comp1, comp2, D_Pos, D_Neg, Score) %>%
  arrange(desc(Score))

candidate_scores$disorder<-c('AD')
cache<-rbind(cache,candidate_scores)

##ADHD


positive_controls <- c('FOXP2','ST3GAL3')  #DRD4
negative_controls <- c("TYR", "FABP4", "ALB", "TNNT3", "OR5M3")
Candidate<-unique(scoretable$RMP[scoretable$total.score>1])


filesave_wide <- filesave[filesave$Disorder %in% c('ADHD'),] %>% #
  pivot_wider(
    id_cols = Gene,
    names_from = c(Disorder, Method),
    names_sep = "_",
    values_from = Score)

filesave_wide<-filesave_wide[filesave_wide$Gene %in% c(positive_controls, negative_controls, Candidate),]
filesave_wide <- filesave_wide %>%
  select_if(~ !is.numeric(.) || (sum(!is.na(.)) > 4 && sum(. != 0, na.rm = TRUE) > 0))

pca_matrix<-as.data.frame(filesave_wide)
sczh2<-read.csv('adhdh2_result.csv')
sczh2<-sczh2[match(pca_matrix$Gene,sczh2$GENE),]
sczh2$ratio<-sczh2$H2/sczh2$length
pca_matrix$sch2<-scale(sczh2$pve)
pca_matrix$sch2ratio<-scale(sczh2$ratio)

rownames(pca_matrix)<-pca_matrix$Gene
pca_matrix<-pca_matrix[,-1]
pca_matrix[is.na(pca_matrix)]<-0
pca_matrix[pca_matrix=='Inf']<-16
pca_matrix<-abs(pca_matrix)


ppca_res <- pcaMethods::pca(pca_matrix,
                            method = "ppca", nPcs = 2,
                            scale = c("pareto"), center = TRUE,)

completed_data <- as.data.frame(ppca_res@completeObs) %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(gene_annotation, by = "Gene")

control_genes <- c(positive_controls, negative_controls) 

X_train <- completed_data[completed_data$Gene %in% control_genes, ]
rownames(X_train) <- X_train$Gene

X_train <- as.matrix(X_train %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

train_labels <- ifelse(rownames(X_train) %in% positive_controls, "Positive", "Negative")
Y_train <- as.factor(train_labels)

X_test <- completed_data[!completed_data$Gene %in% control_genes, ]
rownames(X_test) <- X_test$Gene

X_test <- as.matrix(X_test %>% dplyr::select(-any_of(c("Gene", "hgnc_symbol", "ensembl_gene_id","Group"))))

splsda_model <- splsda(X_train, Y_train, ncomp = 2, keepX = c(5, 5))

splsda_pred <- predict(splsda_model, newdata = X_test)

train_coords <- as.data.frame(splsda_model$variates$X)
train_coords$Gene <- rownames(train_coords)
train_coords$Group <- train_labels

test_coords <- as.data.frame(splsda_pred$variates)
test_coords$Gene <- rownames(X_test)
test_coords$Group <- "Candidate Gene"

colnames(test_coords)<-c("comp1", "comp2", "Gene", "Group")
plsda_plot_df <- rbind(train_coords, test_coords)

expl_vars <- round(splsda_model$prop_expl_var$X * 100, 2)
comp1_label <- paste0("sPLS-DA Component 1 (", expl_vars[1], "%)")
comp2_label <- paste0("sPLS-DA Component 2 (", expl_vars[2], "%)")

#中心点计算
centroids <- plsda_plot_df %>%
  filter(Group %in% c("Positive", "Negative")) %>%
  group_by(Group) %>%
  summarize(Comp1 = mean(comp1), Comp2 = mean(comp2))
pos_center <- c(centroids$Comp1[centroids$Group == "Positive"], centroids$Comp2[centroids$Group == "Positive"])
neg_center <- c(centroids$Comp1[centroids$Group == "Negative"], centroids$Comp2[centroids$Group == "Negative"])

p[[4]]<-
  ggplot() +
  geom_point(data = filter(plsda_plot_df, Group == "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = filter(plsda_plot_df, Group != "Candidate Gene"), 
             aes(x = comp1, y = comp2, color = Group, shape = Group), 
             size = 3.5, alpha = 0.85) +
  geom_point(data = centroids, 
             aes(x = Comp1, y = Comp2, fill = Group), 
             shape = 23, size = 5, color = NULL, stroke = 0) +
  scale_color_manual(values = c("Candidate Gene" = "#8c8c8c", 
                                "Positive" = "#d95f02", 
                                "Negative" = "#7570b3")) +
  scale_fill_manual(values = c("Positive" = "#d95f02", 
                               "Negative" = "#7570b3"), guide = "none") +
  stat_ellipse(data = filter(plsda_plot_df, Group %in% c("Positive", "Negative")),
               aes(x = comp1, y = comp2, fill = Group), 
               geom = "polygon", 
               alpha = 0, level = 0.95, 
               show.legend = FALSE) +
  geom_text_repel(data = plsda_plot_df, 
                  aes(x = comp1, y = comp2, label = Gene), 
                  size = 3.8, max.overlaps = 5, 
                  box.padding = 0.5,
                  point.padding = 0.3,
                  show.legend = FALSE) +
  labs(x = comp1_label, y = comp2_label, 
       title = "ADHD") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )+
  theme(aspect.ratio=1)

#距离分数

W1 <- splsda_model$prop_expl_var$X[1] 
W2 <- splsda_model$prop_expl_var$X[2]  
candidates <- plsda_plot_df %>% filter(Group == "Candidate Gene")
candidate_scores <- candidates %>%
  mutate(
    D_Pos = sqrt( W1 * (comp1 - pos_center[1])^2 + W2 * (comp2 - pos_center[2])^2 ),
    D_Neg = sqrt( W1 * (comp1 - neg_center[1])^2 + W2 * (comp2 - neg_center[2])^2 ),
    Score = D_Neg / (D_Pos + D_Neg)
  ) %>%
  dplyr::select(Gene, comp1, comp2, D_Pos, D_Neg, Score) %>%
  arrange(desc(Score))

candidate_scores$disorder<-c('ADHD')
cache<-rbind(cache,candidate_scores)

write.csv(cache,'./property/riskdistance.csv')

pdf("./property/PLSDA_new.pdf", width = 12, height = 12)
patchwork::wrap_plots(p, nrow = 2)
dev.off()

plot<-merge(cache[cache$disorder=='SCZ',],cache[cache$disorder=='BIPI',],by=c('Gene'))
plot<-plot[,c(1,6,12)]
rownames(plot)<-plot$Gene
plot<-plot[,-1]
plot<-plot[rownames(plot) %in% c('NSUN2','NSUN6','QTRT1','TYW5','TRMT61A','THUMPD3'),]
col_fun1 = colorRamp2(c(0.3,0.5,0.8), c('#22B5AF','white','#F57F17'))

pdf('./property/weigth_distance.pdf',width = 8,height = 8)

Heatmap(
  plot,col = col_fun1,
  name = "Distance",na_col = "gray90",
  cluster_columns = FALSE,  
  cluster_rows = F,    
  show_column_names = T,
  width = ncol(plot)*unit(6, "mm"),
  height = nrow(plot)*unit(4, "mm"))

dev.off()



