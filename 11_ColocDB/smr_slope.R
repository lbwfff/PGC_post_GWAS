
#对应Figure 6a

library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ---------- 处理第一个数据集 (SCZ) ----------
scz_clean <- sczeqtl %>%
  mutate(log10_p = -log10(p_smr)) %>%
  filter(!is.na(log10_p) & !is.infinite(log10_p))

scz_global_max_x <- max(scz_clean$log10_p, na.rm = TRUE)

scz_slopes <- scz_clean %>%
  group_by(gene) %>%
  summarise(
    slope = if (n() == 1) {
      b_smr / log10_p
    } else if (var(log10_p) == 0) {
      mean(b_smr) / mean(log10_p)
    } else {
      coef(lm(b_smr ~ log10_p + 0))[1]},
    .groups = 'drop')

scz_plot <- left_join(scz_clean, scz_slopes, by = "gene")

# ---------- 处理第二个数据集 (BIP) ----------
bip_clean <- bipeqtl %>%
  mutate(log10_p = -log10(p_smr)) %>%
  filter(!is.na(log10_p) & !is.infinite(log10_p))

bip_global_max_x <- max(bip_clean$log10_p, na.rm = TRUE)

bip_slopes <- bip_clean %>%
  group_by(gene) %>%
  summarise(
    slope = if (n() == 1) {
      b_smr / log10_p
    } else if (var(log10_p) == 0) {
      mean(b_smr) / mean(log10_p)
    } else {
      coef(lm(b_smr ~ log10_p + 0))[1]},
    .groups = 'drop')

bip_plot <- left_join(bip_clean, bip_slopes, by = "gene")

# ---------- 构建统一的颜色映射 ----------
all_genes <- unique(c(scz_plot$gene, bip_plot$gene))
n_genes <- length(all_genes)

# 生成足够多的颜色（从 Set1 扩展）
if (n_genes <= 9) {
  gene_colors <- brewer.pal(max(3, n_genes), "Set1")
} else {
  gene_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_genes)
}
names(gene_colors) <- all_genes

# 将 gene 转为 factor 以保证两个图 levels 一致
scz_plot$gene <- factor(scz_plot$gene, levels = all_genes)
bip_plot$gene <- factor(bip_plot$gene, levels = all_genes)

# ---------- 绘图 ----------
p <- list()

p[[2]] <-
  ggplot(scz_plot, aes(x = log10_p, y = b_smr, color = gene)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_segment(
    aes(x = 0, y = 0, 
        xend = scz_global_max_x, 
        yend = slope * scz_global_max_x, 
        color = gene),
    linewidth = 0.6, alpha = 0.8, linetype = "dashed") +
  expand_limits(x = 0, y = 0) + 
  coord_cartesian(ylim = c(-0.5, 0)) +
  scale_color_manual(values = gene_colors) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 14) +
  labs(
    x = "-log10(P value)",
    y = "Beta value",
    title = "eQTL x SCZ GWAS") +
  theme(panel.border = element_rect(linewidth = 0.8, fill = 'transparent'), 
        axis.ticks = element_line(linewidth = 0.6),  
        legend.background = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.2, "cm")) + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))

p[[1]] <-
  ggplot(bip_plot, aes(x = log10_p, y = b_smr, color = gene)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_segment(
    aes(x = 0, y = 0, 
        xend = bip_global_max_x, 
        yend = slope * bip_global_max_x, 
        color = gene),
    linewidth = 0.6, alpha = 0.8, linetype = "dashed") +
  expand_limits(x = 0, y = 0) + 
  coord_cartesian(ylim = c(-1.25, 0)) +
  scale_color_manual(values = gene_colors) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 14) +
  labs(
    x = "-log10(P value)",
    y = "Beta value",
    title = "eQTL x BIP GWAS") +
  theme(panel.border = element_rect(linewidth = 0.8, fill = 'transparent'), 
        axis.ticks = element_line(linewidth = 0.6),  
        legend.background = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.2, "cm")) + 
  theme(aspect.ratio = 1) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))

#

pdf("./property/SMR_eQTL.pdf",width = 12,height = 8)
wrap_plots(p,nrow=1) 
dev.off()
