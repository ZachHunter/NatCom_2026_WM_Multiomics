library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ggrastr)
library(rstatix)
library(bvt)


#============================#
# The Evolution and Subtypes of Waldenstrom Macroglobulinemia:
# Findings from a Multi-omic Analysis of 249 Treatment Naive MYD88L265P Mutated Patients
#
# This R Script includes code used to generate:
#   Figure panel 9
#   Supplemental Figure 5B
#
# Single Cell analysis courtesy of Dr. Hao Sun
# Data from Sun H,  et al. Evolution of tumor subclones and T-cell dynamics underlie variable ibrutinib responses in Waldenström macroglobulinemia.
# Blood. Published online February 11, 2026. doi:10.1182/blood.2025032268
# Data archived available at the European Genome-Phenome Archive - Data set: EGAD50000002279
# https://ega-archive.org/dacs/EGAC50000000863
#============================#



##############################
# File dependencies
# Loading required files for analysis
##############################

# Set Input/Output directories
dataDir<-"~/Path/to/Directory/Data/"
outputDir<-"~/Path/to/Directory/"

# table of precalculated EScore, BCL, PCL and Subtype Score values by patient sample
df_score <- read.csv(file.path(dataDir,"BPsubtype_ESscore.csv"))
# pseudobulk expression of genes of interest by cell type with EScore annotation
df_wide <- read.csv(file.path(dataDir,"validation_genes_log2CPM.csv"), check.names = FALSE)
# Percentage of cells in group with detectable expression by cell type
df_wideP <- read.csv(file.path(dataDir,"validation_genes_PercentPositiveCells.csv"), check.names = FALSE)
# UMAP positions for single cell analysis of clonal WM LPCs, HD Memory B-cells, and Non-clonal patient Memory B-cells
df <- read.csv(file.path(dataDir,"UMAP_HD_early_late_WMcell_memoryB.csv"))

WMSubtypeColor<-c(3,1,2,4:9)
SubtypeTheme<-newNPTheme(npDefaultTheme, plotColors=list(points=npDefaultTheme$plotColors$points[WMSubtypeColor],fill=npDefaultTheme$plotColors$fill[WMSubtypeColor] ))


#-------------------------------
# Figure 9.1 kmeans of subtype and escore signatures
# K-means clustering on (subtype_score, Escore)
# ------------------------------
# 1.1 Prepare matrix & scale
# ------------------------------
X_scaled <- df_score %>%
  select(subtype_score, Escore) %>%
  as.matrix() %>%
  scale()

# ------------------------------
# 1.2 K-means clustering (K = 3)
# ------------------------------
set.seed(123)
km <- kmeans(X_scaled, centers = 3, nstart = 50)

df_score <- df_score %>%
  mutate(cluster = factor(paste0("cluster ", km$cluster)))  # cluster 1/2/3

# ------------------------------
# 1.3 Plot
# ------------------------------

pdf(file=file.path(outputDir,"Figures/Figure9/F9A_SC_EScore_Subtype_Kmeans.pdf"), width = 6.5, height = 5)
ggplot(df_score, aes(x = subtype_score, y = -Escore, color = cluster)) +
  scale_color_manual(values = SubtypeTheme$plotColors$fill) +
  scale_fill_manual(values = SubtypeTheme$plotColors$fill) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf) +
  stat_ellipse(
    aes(fill = cluster, color = cluster),
    geom = "polygon",
    alpha = 0.2,
    level = 0.95,
    linewidth = 1,
    show.legend = FALSE
  ) +
  labs(
    x = "PCL \u2190  subtype_score  \u2192  BCL",
    y = "early \u2190  Escore  \u2192  late",
    color = "K-means"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = -60, y = 15, label = "PCL",
           size = 5, fontface = "plain") +

  annotate("text", x = 70, y = 12, label = "BCL",
           size = 5, fontface = "plain") +

  annotate("text", x = -20, y = -35, label = "Early EScore",
           size = 5, fontface = "plain")
dev.off()


#-------------------------------
# Figure 9.2 log2 cpm 4 groups
# EScore genes in Early/Late + Controls
# ------------------------------

# ==============================
# 2.1 Wide -> long, set group order
# ==============================
group_levels <- c("HD-MemoryB", "Patient-MemoryB", "early-EScore-WM", "late-EScore-WM")

#Renaming Levels
df_wide$sample_cell_group[df_wide$sample_cell_group=="WM-MemoryB"]<-group_levels[2]
df_wide$sample_cell_group[df_wide$sample_cell_group=="earlyWM-WMcells"]<-group_levels[3]
df_wide$sample_cell_group[df_wide$sample_cell_group=="lateWM-WMcells"]<-group_levels[4]

df_wide <- df_wide %>%
  mutate(sample_cell_group = factor(sample_cell_group, levels = group_levels))

# ==============================
# 2.2 Wide -> long, pivot ALL numeric columns
# ==============================
df_long <- df_wide %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "gene",
    values_to = "log2CPM"
  ) %>%
  filter(
    !is.na(gene),
    !is.na(sample_cell_group),
    !is.na(log2CPM),
    is.finite(log2CPM)
  ) %>%
  mutate(
    gene = as.character(gene),
    sample_cell_group = droplevels(sample_cell_group)
  )

# ==============================
# 2.3 QC for statistical tests
#    Require per gene x group:
#    - n >= 2
#    - >=2 unique values
# ==============================
group_ok <- df_long %>%
  group_by(gene, sample_cell_group) %>%
  summarise(
    n = n(),
    n_uniq = n_distinct(log2CPM),
    .groups = "drop"
  ) %>%
  mutate(ok = (n >= 2 & n_uniq >= 2))

df_ok <- df_long %>%
  inner_join(
    group_ok %>% filter(ok) %>% select(gene, sample_cell_group),
    by = c("gene", "sample_cell_group")
  )

# ==============================
# 2.4 Wilcoxon (BH), keep Early vs Late only
# ==============================
stat_test <- df_ok %>%
  group_by(gene) %>%
  pairwise_wilcox_test(
    log2CPM ~ sample_cell_group,
    p.adjust.method = "BH",
    paired = FALSE,
    exact = FALSE
  ) %>%
  ungroup() %>%
  filter(group1 == "early-EScore-WM", group2 == "late-EScore-WM") %>%
  filter(p.adj.signif != "ns")

stat_test2 <- df_ok %>%
  group_by(gene) %>%
  pairwise_wilcox_test(
    log2CPM ~ sample_cell_group,
    p.adjust.method = "BH",
    paired = FALSE,
    exact = FALSE
  ) %>%
  ungroup() %>%
  filter(group1 == "Patient-MemoryB", group2 == "early-EScore-WM") %>%
  filter(p.adj.signif != "ns")

# Auto y-position per gene panel
y_pos <- df_long %>%
  group_by(gene) %>%
  summarise(y.position = max(log2CPM, na.rm = TRUE) * 0.9, .groups = "drop")

stat_test <- stat_test %>%
  left_join(y_pos, by = "gene") %>%
  mutate(
    xmin = "early-EScore-WM",
    xmax = "late-EScore-WM"
  )
stat_test2 <- stat_test2 %>%
  left_join(y_pos, by = "gene") %>%
  mutate(
    xmin = "Patient-MemoryB",
    xmax = "early-EScore-WM"
  )
stat_test<- stat_test %>%
  bind_rows(stat_test2)

# ==============================
# 2.5 Plot
# ==============================
pal <- c("HD-MemoryB"       = "#b5dffd",
         "Patient-MemoryB"       = "#d098ee",
         "early-EScore-WM"  = "#acd98d",
         "late-EScore-WM"   = "#f4737a")

pdf(file=file.path(outputDir,"Figures/Figure9/F9B_SC_EScore_CpM_4way.pdf"), width = 6, height = 4.5)
ggplot(df_long, aes(x = sample_cell_group, y = log2CPM)) +
  geom_boxplot(
    aes(fill = sample_cell_group),
    width = 0.55,
    outlier.shape = NA,
    alpha = 1
  ) +
  geom_point(
    aes(color = sample_cell_group),
    position = position_jitter(width = 0.18),
    size = 2, alpha = 0.3
  ) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.01,
    inherit.aes = FALSE,
    hide.ns = TRUE
  ) +
  labs(x = "", y = "Gene expression (log2 CPM)", fill = 'Cell group') +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )+
  guides(color = "none")
dev.off()

#-------------------------------
# Sup Figure 5.2 %Positive Cells 4 groups
# EScore genes in Early/Late + Controls
# ------------------------------

# ==============================
# S5.2.1 Wide -> long, pivot ALL numeric columns
# ==============================


#Renaming Levels
df_wideP$sample_cell_group[df_wideP$sample_cell_group=="WM-MemoryB"]<-group_levels[2]
df_wideP$sample_cell_group[df_wideP$sample_cell_group=="earlyWM-WMcells"]<-group_levels[3]
df_wideP$sample_cell_group[df_wideP$sample_cell_group=="lateWM-WMcells"]<-group_levels[4]
df_wideP$sample_cell_group<-factor(df_wideP$sample_cell_group,levels=group_levels)

df_long <- df_wideP %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "gene",
    values_to = "percentPositive"
  ) %>%
  filter(
    !is.na(gene),
    !is.na(sample_cell_group),
    !is.na(percentPositive),
    is.finite(percentPositive)
  ) %>%
  mutate(
    gene = as.character(gene),
    sample_cell_group = droplevels(sample_cell_group)
  )

# ==============================
# S5.2.2 QC for statistical tests
#    Require per gene x group:
#    - n >= 2
#    - >=2 unique values
# ==============================
group_ok <- df_long %>%
  group_by(gene, sample_cell_group) %>%
  summarise(
    n = n(),
    n_uniq = n_distinct(percentPositive),
    .groups = "drop"
  ) %>%
  mutate(ok = (n >= 2 & n_uniq >= 2))

df_ok <- df_long %>%
  inner_join(
    group_ok %>% filter(ok) %>% select(gene, sample_cell_group),
    by = c("gene", "sample_cell_group")
  )

# ==============================
# S5.2.3 Wilcoxon (BH), keep Early vs Late only
# ==============================
stat_test <- df_ok %>%
  group_by(gene) %>%
  pairwise_wilcox_test(
    percentPositive ~ sample_cell_group,
    p.adjust.method = "BH",
    paired = FALSE,
    exact = FALSE
  ) %>%
  ungroup() %>%
  filter(group1 == "early-EScore-WM", group2 == "late-EScore-WM") %>%
  filter(p.adj.signif != "ns")

stat_test2 <- df_ok %>%
  group_by(gene) %>%
  pairwise_wilcox_test(
    percentPositive ~ sample_cell_group,
    p.adjust.method = "BH",
    paired = FALSE,
    exact = FALSE
  ) %>%
  ungroup() %>%
  filter(group1 == "Patient-MemoryB", group2 == "early-EScore-WM") %>%
  filter(p.adj.signif != "ns")

# Auto y-position per gene panel
y_pos <- df_long %>%
  group_by(gene) %>%
  summarise(y.position = max(percentPositive, na.rm = TRUE) * 0.9, .groups = "drop")

stat_test <- stat_test %>%
  left_join(y_pos, by = "gene") %>%
  mutate(
    xmin = "early-EScore-WM",
    xmax = "late-EScore-WM"
  )
stat_test2 <- stat_test2 %>%
  left_join(y_pos, by = "gene") %>%
  mutate(
    xmin = "Patient-MemoryB",
    xmax = "early-EScore-WM"
  )
stat_test<- stat_test %>% bind_rows(stat_test2)

# ==============================
# S5.2.4 Plot
# ==============================
pal <- c("HD-MemoryB"       = "#b5dffd",
         "Patient-MemoryB"       = "#d098ee",
         "early-EScore-WM"  = "#acd98d",
         "late-EScore-WM"   = "#f4737a")

pdf(file=file.path(outputDir,"Figures/SFigure5/SF5B_SC_EScore_PercPos_4way.pdf"), width = 6, height = 4.5)
ggplot(df_long, aes(x = sample_cell_group, y = percentPositive)) +
  geom_boxplot(
    aes(fill = sample_cell_group),
    width = 0.55,
    outlier.shape = NA,
    alpha = 1
  ) +
  geom_point(
    aes(color = sample_cell_group),
    position = position_jitter(width = 0.18),
    size = 2, alpha = 0.3
  ) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    tip.length = 0.01,
    inherit.aes = FALSE,
    hide.ns = TRUE
  ) +
  labs(x = "", y = "Percent of positive cells", fill = 'Cell group') +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )+
  guides(color = "none")
dev.off()



#-------------------------------
# Figure 9.3 Single Cell UMAP
# Colored by imputed WM Subtype
# ------------------------------

# set cell_group
tempFact<-as.character(df$HD_early_late_WMcell_memoryB)
tempFact[tempFact=="Early WM-WM cells"] <- "Early EScore-WM cells"
df$HD_early_late_WMcell_memoryB <- factor(
  tempFact,
  levels = c("HD-Memory B cells", "WM-Memory B cells",
             "Early EScore-WM cells", "BCL-WM cells", "PCL-WM cells")
)

# ==============================
# 2) Plot
# ==============================

pdf(file=file.path(outputDir,"Figures/Figure9/F9C_SC_UMAP.pdf"), width = 6, height = 4.5)
ggplot(df, aes(x = UMAP1, y = UMAP2, color = HD_early_late_WMcell_memoryB)) +
  geom_point_rast(size = 0.01, alpha = 0.5) +
  labs(
    title = "Cell annotation group",
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = c('#b5dffd','#d098ee','#acd98d','#f4737a', '#2c69b0'))+
  guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()
