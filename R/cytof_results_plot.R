#' ---
#' title: "Visualization of mass cytometry analysis results: marker intensity in tSNE, barplot per patient, etc
#' author: "Fan Zhang"
#' date: "2018-05-30"
#' ---


setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

load("../data/synData.Fibro.downsample.SNE.RData")
load("../data/synData.Bcell.downsample.SNE.RData")
load("../data/synData.Mono.downsample.SNE.RData")
load("../data/synData.Tcell.downsample.SNE.RData")

# 300-0486C, 300-0511A
synData.Fibro.downsample$sampleID[which(synData.Fibro.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Fibro.downsample$sampleID[which(synData.Fibro.downsample$sampleID == "300-0511A")] <- "300-0511"
synData.Bcell.downsample$sampleID[which(synData.Bcell.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Bcell.downsample$sampleID[which(synData.Bcell.downsample$sampleID == "300-0511A")] <- "300-0511"
synData.Mono.downsample$sampleID[which(synData.Mono.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Mono.downsample$sampleID[which(synData.Mono.downsample$sampleID == "300-0511A")] <- "300-0511"
synData.Tcell.downsample$sampleID[which(synData.Tcell.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Tcell.downsample$sampleID[which(synData.Tcell.downsample$sampleID == "300-0511A")] <- "300-0511"

pacman::p_load(
  ggplot2,
  patchwork,
  dplyr,
  magrittr,
  cetcolor,
  seriation,
  data.table,
  pheatmap,
  ggrepel
)

# Plot tSNE for merged clusters for each cell type
# B cell
cluster_center <- synData.Bcell.downsample %>%
                  group_by(markers) %>%
                  summarise_at(vars(SNE1, SNE2), funs(median(., na.rm=TRUE)))
cluster_center <- as.data.frame(cluster_center)


ggplot() +
  geom_point(
    data = synData.Bcell.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "markers"),
    size = 1.5, stroke = 0.05, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$cytof_cluster, name = "") +
  geom_label_repel(
    data = cluster_center,
    aes(x = SNE1, y = SNE2, label = markers),
    size = 5, color = "black",
    fontface = 'bold',
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("bcell_cytof_tSNE", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()

# ---
# T cell
cluster_center <- synData.Tcell.downsample %>% group_by(markers) %>% summarise(
  mean_SNE1 = mean(SNE1),
  mean_SNE2 = mean(SNE2)
) 
cluster_center <- as.data.frame(cluster_center)

ggplot() +
  geom_point(
    data = synData.Tcell.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "markers"),
    size = 1.2, stroke = 0.05, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$cytof_cluster, name = "") +
  geom_label_repel(
    data = cluster_center,
    aes(x = mean_SNE1, y = mean_SNE2, label = markers),
    size = 5, color = "black",
    fontface = 'bold',
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("tcell_cytof_tSNE", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()

# ---
# Monocyte
cluster_center <- synData.Mono.downsample %>% group_by(markers) %>% summarise(
  mean_SNE1 = mean(SNE1),
  mean_SNE2 = mean(SNE2)
) 
cluster_center <- as.data.frame(cluster_center)

ggplot() +
  geom_point(
    data = synData.Mono.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "markers"),
    size = 1, stroke = 0.05, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$cytof_cluster, name = "") +
  geom_label_repel(
    data = cluster_center,
    aes(x = mean_SNE1, y = mean_SNE2, label = markers),
    size = 5, color = "black",
    fontface = 'bold',
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("mono_cytof_tSNE", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()


# ---
# Fibroblast
cluster_center <- synData.Fibro.downsample %>% group_by(markers) %>% summarise(
  mean_SNE1 = mean(SNE1),
  mean_SNE2 = mean(SNE2)
) 
cluster_center <- as.data.frame(cluster_center)

ggplot() +
  geom_point(
    data = synData.Fibro.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "markers"),
    size = 1.1, stroke = 0.05, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$cytof_cluster, name = "") +
  # scale_color_manual(values = meta_colors$cytof_cluster, name = "") +
  # geom_label_repel(
  #   data = cluster_center,
  #   aes(x = mean_SNE1, y = mean_SNE2, label = markers, color = markers),
  #   size = 5, #color = "black",
  #   fontface = 'bold',
  #   box.padding = unit(0.8, "lines"),
  #   point.padding = unit(0.3, "lines")
  # ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("fibro_cytof_tSNE", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()

# -------
# Plot barplots per patient
d <- readRDS("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/data/cytof_plot.rds")
d %<>% group_by(donor) %>% mutate(percent_donor = 100 * cells / sum(cells))
d %<>% group_by(donor, cell_type) %>% mutate(percent_donor_cell_type = 100 * cells / sum(cells))
d %<>% group_by(cluster, cell_type) %>% mutate(percent_cluster = 100 * cells / sum(cells))
head(d)

seriate_cols <- function(d, col1, col2, value.var = "percent") {
  mat <- dcast(
    data = d,
    formula = as.formula(sprintf("%s ~ %s", col1, col2)),
    value.var = value.var
  )
  rownames(mat) <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat_order <- seriate(mat, method = "BEA_TSP")
  d[[col1]] <- factor(as.character(d[[col1]]), rownames(mat)[mat_order[[1]]])
  d[[col2]] <- factor(as.character(d[[col2]]), colnames(mat)[mat_order[[2]]])
  return(d)
}

options(repr.plot.width = 12, repr.plot.height = 4)

plot_celltype_bars <- function(this_celltype = "Fibroblast") {
  this_d <- subset(d, cell_type == this_celltype)
  this_d <- seriate_cols(this_d, "donor", "cluster", "cells")
  this_d %<>% arrange(cluster)
  
  ggplot(
    data = this_d,
    mapping = aes(x = donor, y = cells, fill = cluster)
    ) +
    geom_col() +
    facet_grid( ~ Mahalanobis_20, scales = "free", space = "free_x") +
    scale_fill_manual(
      name = "",
      #values = alpha(meta_colors$cytof_cluster, 0.9),
      values = meta_colors$cytof_cluster,
      breaks = levels(this_d$cluster),
      labels = levels(this_d$cluster)
    ) +
    labs(
      x = "Donors",
      y = "Number of cells",
      title = this_celltype
    ) +
    theme_bw(base_size = 20) +
    theme(panel.border = element_rect(size = 0.2),
          axis.ticks = element_line(size = 0.2),
          strip.background = element_blank(),
          # axis.ticks = element_blank(), 
          panel.grid = element_blank(),
          axis.text = element_text(size = 15),
          axis.text.x=element_text(angle = 35, hjust = 1),
          axis.text.y = element_text(size = 15)
    ) 
}

plot_celltype_bars("Fibroblast")
ggsave(file = paste("fibro_cytof_bplot", ".pdf", sep = ""), width = 14, height = 4.5) # , useDingbats = FALSE
dev.off()

plot_celltype_bars("B cell")
ggsave(file = paste("bcell_cytof_bplot", ".pdf", sep = ""), width = 14, height = 4.5) # , useDingbats = FALSE
dev.off()

plot_celltype_bars("T cell")
ggsave(file = paste("tcell_cytof_bplot", ".pdf", sep = ""), width = 13.5, height = 4.5) # , useDingbats = FALSE
dev.off()

plot_celltype_bars("Monocyte")
ggsave(file = paste("mono_cytof_bplot", ".pdf", sep = ""), width = 13.5, height = 4.5) # , useDingbats = FALSE
dev.off()


# ------
# Fibroblast
# Plot intensity of a marker
ggplot() +
  geom_point(
    data = synData.Fibro.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "Cadherin.11"), 
    shape = 21, size = 1.1, stroke = 0.011
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9),
    name = "Normalized intensity"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("fibro_cytof_tSNE_Cadherin.11", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()


# B cell
# Plot intensity of a marker
ggplot() +
  geom_point(
    data = synData.Bcell.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "FcRL4"), 
    shape = 21, size = 1.1, stroke = 0.011
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9),
    name = "Normalized intensity"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("bcell_cytof_tSNE_FcRL4", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()

# T cell
# Plot intensity of a marker
ggplot() +
  geom_point(
    data = synData.Tcell.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "FoxP3"), 
    shape = 21, size = 1.1, stroke = 0.011
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9),
    name = "Normalized intensity"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("tcell_cytof_tSNE_FoxP3", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()

# Monocyte
# Plot intensity of a marker
ggplot() +
  geom_point(
    data = synData.Mono.downsample,
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "CCR2"), 
    shape = 21, size = 1.1, stroke = 0.011
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9),
    name = "Normalized intensity"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("mono_cytof_tSNE_CCR2", ".png", sep = ""), width = 5, height = 5) # , useDingbats = FALSE
dev.off()

# ------
# Plot disease status
cells_type <- synData.Tcell.downsample

inflam_label <- read.xls("../data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label$Patient <- as.character(inflam_label$Patient)
inter <- intersect(cells_type$sampleID, inflam_label$Patient)

cells_type <- cells_type[which(cells_type$sampleID %in% inter), ]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter), ]
inflam_label <- inflam_label[, c("Patient", "Mahalanobis_20")]
colnames(inflam_label)[1] <- "sampleID"
cells_type <- merge(cells_type, inflam_label, by = "sampleID")

ggplot() +
  geom_point(
    data = cells_type[sample(nrow(cells_type)),],
    mapping = aes_string(x = "SNE1", y = "SNE2", fill = "Mahalanobis_20"),
    size = 1.4, stroke = 0.0001, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$Case.Control, name = "") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 22) +
  theme(
    legend.position = "none",
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste(gsub(".*[.]([^.]+)[.].*", "\\1", deparse(substitute(synData.Tcell.downsample))), 
                    "_cytof_disease", ".png", sep = ""), width = 5, height = 5) 
dev.off()


