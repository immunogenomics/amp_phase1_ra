# Adress reviewers' comments for AMP Phase I RA paper
# Fan Zhang
# 2018-08-14

setwd("Documents/GitHub/amp_phase1_ra/")

library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
require(gdata)
library(Seurat)
library(pbapply)
library(parallel)
library(forcats)
library(pheatmap)
library(RColorBrewer)

source("R/meta_colors.R")

meta_colors$disease = c(
  "Leukocyte-rich RA" = "#FF7F00",
  "Leukocyte-poor RA" = "#FFD8B2",
  "OA" = "#6A3D9A"
)

# -------------------------
# Load single-cell data
log2cpm <- readRDS("../amp_phase1_ra_viewer/data/celseq_synovium_log2_5265cells_paper.rds")
meta <- readRDS("../amp_phase1_ra_viewer/data/celseq_synovium_meta_5265cells_paper.rds")
all(colnames(log2cpm) == meta$cell_name)

# Load bulk RNA-seq data
log2tpm <- readRDS("data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)

# Add manhalanobis labels
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label$Patient <- as.character(inflam_label$Patient)
table(inflam_label$Mahalanobis_20)

inter <- intersect(intersect(meta$sample, inflam_label$Patient), bulk_meta$Donor.ID)
# inter <- intersect(meta$sample, inflam_label$Patient)
meta <- meta[which(meta$sample %in% inter),]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
inflam_label <- inflam_label[, c(1,42)]

colnames(inflam_label)[1] <- "sample"
meta <- merge(meta, inflam_label, by ="sample")
dim(meta)
table(meta$Mahalanobis_20)

log2cpm <- log2cpm[, which(colnames(log2cpm) %in% meta$cell_name)]

meta <- meta[order(match(meta$cell_name, colnames(log2cpm))), ]
all(colnames(log2cpm) == meta$cell_name)


# -------------------------
# Use intersected bulk samples
# bulk_meta <- bulk_meta[which(bulk_meta$Donor.ID %in% inter),]
# log2tpm <- log2tpm[, which(colnames(log2tpm) %in% bulk_meta$Sample.ID)]
# all(colnames(log2tpm) == bulk_meta$Sample.ID)

cell_type <- "B cell"
log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_fibro) == bulk_meta_fibro$Sample.ID)
dim(bulk_meta_fibro)
dim(log2tpm_fibro)


# Add manhalanobis labels
# inflam_label <- read.xls("data/postQC_all_samples.xlsx")
# inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
# inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
# inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
# inflam_label$Patient <- as.character(inflam_label$Patient)
# table(inflam_label$Mahalanobis_20)
# inter <- intersect(bulk_meta_fibro$Donor.ID, inflam_label$Patient)

bulk_meta_fibro <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% inter),]
inflam_label <- inflam_label[which(inflam_label$sample %in% inter),]

bulk_meta_fibro <- bulk_meta_fibro[order(match(bulk_meta_fibro$Donor.ID, inflam_label$sample)), ]
all(bulk_meta_fibro$Donor.ID == inflam_label$sample)

bulk_meta_fibro$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(bulk_meta_fibro$Mahalanobis_20)

log2tpm_fibro <- log2tpm_fibro[, which(colnames(log2tpm_fibro) %in% bulk_meta_fibro$Sample.ID)]
log2tpm_fibro <- log2tpm_fibro[, order(match(colnames(log2tpm_fibro), bulk_meta_fibro$Sample.ID))]
all(bulk_meta_fibro$Sample.ID == colnames(log2tpm_fibro))
dim(log2tpm_fibro)


# ------------------------------
# For each cell type data
# scRNA-seq
log2cpm_fibro <- log2cpm[, which(meta$cell_type == "Fibroblast")]
meta_fibro <- meta[which(meta$cell_type == "Fibroblast"),]
all(colnames(log2cpm_fibro) == meta_fibro$cell_name)

log2cpm_tcell <- log2cpm[, which(meta$cell_type == "T cell")]
meta_tcell <- meta[which(meta$cell_type == "T cell"),]
all(colnames(log2cpm_tcell) == meta_tcell$cell_name)


log2cpm_bcell <- log2cpm[, which(meta$cell_type == "B cell")]
meta_bcell <- meta[which(meta$cell_type == "B cell"),]
all(colnames(log2cpm_bcell) == meta_bcell$cell_name)

log2cpm_monocyte <- log2cpm[, which(meta$cell_type == "Monocyte")]
meta_monocyte <- meta[which(meta$cell_type == "Monocyte"),]
all(colnames(log2cpm_monocyte) == meta_monocyte$cell_name)



# ------------------------------
# For one cluster
subset <- "SC-F2"
ind <- which(meta_fibro$cluster %in% subset)
log2cpm_fibro <- log2cpm_fibro[, ind]
meta_fibro <- meta_fibro[ind, ]
# log2cpm_fibro <- log2cpm_fibro[, order(match(colnames(log2cpm_fibro), meta_fibro$cell_name))]
meta_fibro <- meta_fibro[order(meta_fibro$cell_name),]
log2cpm_fibro <- log2cpm_fibro[, order(colnames(log2cpm_fibro))]
all(colnames(log2cpm_fibro) == meta_fibro$cell_name)

gene <- "HLA-DRA"
table(meta_fibro$cluster)
exp <- log2cpm_fibro[which(rownames(log2cpm_fibro) == gene),]
meta_fibro$gene <- as.numeric(exp)


dat_percent <- meta_fibro %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(percent = sum(gene > 0) / length(gene) * 100,
                   ave_mean = mean(gene),
                   ave_median = median(gene))
dat_percent <- as.data.frame(dat_percent)
dat_percent <- dat_percent[order(dat_percent$sample),]
dat_percent[1:4,]

# Add number percent of cells in SC-F2 for each donor
temp = data.frame(
  sample = as.character(names(table(meta_fibro$sample))),
  num_cell = as.numeric(table(meta_fibro$sample))
)
temp <- temp[-which(temp$num_cell == 0),]
temp$perc_cell <- 100 * temp$num_cell / sum(temp$num_cell)
all(as.character(dat_percent$sample) == as.character(temp$sample))
dat_percent$perc_cell <- temp$perc_cell

# Load bulk data
log2tpm <- readRDS("../data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("../data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)

log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == "Fibro")]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == "Fibro"),]
all(colnames(log2tpm_fibro) == bulk_meta_fibro$cell_name)

gene <- "HLA.DRA"
exp <- log2tpm_fibro[which(rownames(log2tpm_fibro) == gene),]
bulk_meta_fibro$gene_bulk <- as.numeric(exp) 

# ---------------------------------
# Intersect
inter <- intersect(meta_fibro$sample, bulk_meta_fibro$Donor.ID)

meta_1 <- dat_percent[which(dat_percent$sample %in% inter),]
meta_2 <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% inter),]
meta_1$sample <- as.character(meta_1$sample)

# meta_2 <- meta_2[-which(duplicated(meta_2$Donor.ID)),]

meta_1 <- meta_1[ order(match(meta_1$sample, meta_2$Donor.ID)), ]
all(meta_1$sample == meta_2$Donor.ID)
meta_1$gene_bulk <- meta_2$gene_bulk
meta_1$Sublining.fibroblasts <- meta_2$Sublining.fibroblasts
meta_1$Lining.fibroblasts <- meta_2$Lining.fibroblasts

inflam_label <- read.xls("../data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label <- inflam_label[which(inflam_label$Patient %in% meta_1$sample),]
inflam_label$Patient <- as.character(inflam_label$Patient)
meta_1 <- meta_1[ order(match(meta_1$sample, inflam_label$Patient)), ]
all(meta_1$sample == inflam_label$Patient)
meta_1$Mahalanobis <- inflam_label$Mahalanobis_20
meta_1[1:4, ]

ggplot() +
  geom_point(
    data = meta_1,
    mapping = aes_string(x = "gene_bulk", y = "perc_cell", fill = "Mahalanobis"),
    size = 5, stroke = 0.2, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "HLA-DRA expression by bulk RNA-seq",
    y = "Percent of fibroblasts in SC-F2 (%)",
    # y = "HLA-DRA mean expression for\n SC-F2 fibroblasts by scRNA-seq",
    # y = "Percent of HLA-DRA+ cells in \nall fibroblasts by scRNA-seq",
    title = "HLA-DRA"
    # subtitle = ""
  ) +
  theme_bw(base_size = 16) +
  theme(
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("bulk_perc_SC-F2_HLA-DRA", ".png", sep = ""),
       width = 6.2, height = 4, dpi = 300)
dev.off()





# ----------------------------
# Bulk vs flow


ggplot() +
  geom_point(
    data = bulk_meta_fibro,
    mapping = aes_string(x = "THY1", y = "Sublining.fibroblasts", fill = "Mahalanobis"),
    size = 5, stroke = 0.2, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "Gene expression by bulk RNA-seq",
    y = "THY1+/live cells by flow",
    title = ""
    # subtitle = ""
  ) +
  theme_bw(base_size = 16) +
  theme(
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("bulk_flow_THY1+", ".png", sep = ""),
       width = 6.2, height = 4, dpi = 300)
dev.off()


# ---------
# Update Fig 6a T cell heatmap

fan_genes <- c("IL7R", "CD4", "SELL","CD40LG", "AQP3", "ANK3", "TCF7",  "CCR7", "SC5D", "NFKBIZ", "TNFRSF25", "LEF1",
               "FOXP3", "TOMM7", "IKZF2", "LAYN", "FCRL3", "STAM", "CTLA4", "TIGIT", "DUSP4",
               "CXCL13", "PDCD1", "CD200",  # "MAF", "CXCR5",
               "GZMA", "CCL5", "CCL4", "CD8A", "NKG7", "CST7", "SLAMF7", "CRTAM", "GZMK",
               "GNLY", "FGFBP2", "CX3CR1", "TGFBR3", "GZMB", "ZNF683", "SPON2", "PRF1",
               # TpH are cells: PDCD1+, CXCL13+, ICOS+, CXCR5-
               "HLA-DQA1", "HLA-DRB5", "HLA-DRA", "APOBEC3G")

exp <- log2cpm_tcell[fan_genes,]
dim(exp)

mat_breaks <- seq(min(exp), max(exp), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(exp, n = 11)

annotation_col <- meta_tcell[, c("plate", "Mahalanobis_20", "cluster")]
colnames(annotation_col)[3] <- "fine_cluster"
colnames(annotation_col)[2] <- "disease"
rownames(annotation_col) <- meta_tcell$cell_name
rownames(meta_tcell) <- meta_tcell$cell_name
exp <- exp[,order(annotation_col$fine_cluster)]
scale_rows <- function(x) t(scale(t(x)))
exp <- scale_rows(exp) # Z-score
exp[exp > 2] <- 2
exp[exp < -2] <- -2

pdf("heatmap_markers_tcell_7clusters_v3.pdf", width=6.5, height=4, onefile = FALSE, bg = "white")
pheatmap(
  mat = exp,
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(6),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = meta_colors,
  fontsize = 7,
  fontsize_row = 7,
  scale = "none"
)
dev.off()

# ---------------------------------------------
# Update Fig 4a
sc_genes <- c("PTGFR", "PLA2G2A", "RPS19", "FOS", "SFRP1", "C3", "F3", "GAS6", "SFRP2", "TMEM150C", "FBLN5", "ABCA6",
               "CD34", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "IFI30", "RSPO3", "PLAU", "RARRES3", "JAK2",
               "IL6", "DKK3", "CADM1", "CAPG", "AKR1C2", "COL8A2", "ITGA11", "SFRP4", "LDLRAD4", "HBEGF", "CLIC5",
               "C10orf105", "HTRA4", "PCSK6", "DEFB1", "ITGA6", "ERRFI1", "NTN4", "MET")

sc_exp <- log2cpm_fibro[sc_genes,]
dim(sc_exp)
mat_breaks <- seq(min(sc_exp), max(sc_exp), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(sc_exp, n = 11)

annotation_col <- meta_fibro[, c("plate", "Mahalanobis_20", "cluster")]
colnames(annotation_col)[3] <- "fine_cluster"
colnames(annotation_col)[2] <- "disease"
rownames(annotation_col) <- meta_fibro$cell_name
rownames(meta_fibro) <- meta_fibro$cell_name
sc_exp <- sc_exp[,order(annotation_col$fine_cluster)]
scale_rows <- function(x) t(scale(t(x)))
sc_exp <- scale_rows(sc_exp) # Z-score
sc_exp[sc_exp > 2] <- 2
sc_exp[sc_exp < -2] <- -2

pdf("heatmap_markers_fibro.pdf", width=6.5, height=4, onefile = FALSE, bg = "white")
pheatmap(
  mat = exp,
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(6),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = meta_colors,
  fontsize = 7,
  fontsize_row = 7,
  scale = "none"
)
dev.off()

# ------------------------------------
# Both sc + bulk heatmap: not good, less columns for bulk data
bulk_genes <- c("PTGFR", "PLA2G2A", "RPS19", "FOS", "SFRP1", "C3", "F3", "GAS6", "SFRP2", "TMEM150C", "FBLN5", "ABCA6",
               "CD34", "HLA.DRA", "HLA.DPA1", "HLA.DPB1", "HLA.DRB1", "IFI30", "RSPO3", "PLAU", "RARRES3", "JAK2",
               "IL6", "DKK3", "CADM1", "CAPG", "AKR1C2", "COL8A2", "ITGA11", "SFRP4", "LDLRAD4", "HBEGF", "CLIC5",
               "C10orf105", "HTRA4", "PCSK6", "DEFB1", "ITGA6", "ERRFI1", "NTN4", "MET")
bulk_exp <- log2tpm_fibro[bulk_genes,]
dim(bulk_exp)

all(rownames(sc_exp) == rownames(bulk_exp))
rownames(bulk_exp) <- rownames(sc_exp)
both_exp <- cbind.data.frame(sc_exp, bulk_exp)
colnames(bulk_meta_fibro)[1] <- "sample"
both_meta <- rbind.data.frame(meta_fibro[, c("sample", "Mahalanobis_20")], bulk_meta_fibro[, c("sample", "Mahalanobis_20")])
both_meta$data_type <- c(rep("scRNA-seq", nrow(meta_fibro)), rep("bulk RNA-seq", nrow(bulk_meta_fibro)))
both_meta$cell_sample <- c(meta_fibro$cell_name, bulk_meta_fibro$Sample.ID)
both_meta$fine_cluster <- c(meta_fibro$cluster, rep("bulk RNA-seq", nrow(bulk_meta_fibro)))
all(colnames(both_exp) == both_meta$cell_sample)
dim(both_meta)
dim(both_exp)

mat_breaks <- seq(min(both_exp), max(both_exp), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
# mat_breaks <- quantile_breaks(both_exp, n = 11)

annotation_col <- both_meta[, c("sample", "Mahalanobis_20", "data_type", "fine_cluster")]
colnames(annotation_col)[2] <- "disease"
rownames(annotation_col) <- both_meta$cell_sample
rownames(both_meta) <- both_meta$cell_sample
both_exp <- both_exp[, order(annotation_col$fine_cluster)]
scale_rows <- function(x) t(scale(t(x)))
both_exp <- scale_rows(both_exp) # Z-score
both_exp[both_exp > 2] <- 2
both_exp[both_exp < -2] <- -2

pdf("heatmap_markers_fibro.pdf", width=6.5, height=4, onefile = FALSE, bg = "white")
pheatmap(
  mat = both_exp,
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(8),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  # annotation_colors = meta_colors,
  fontsize = 7,
  fontsize_row = 7,
  scale = "none"
)
dev.off()

# ----------------------------------
# Only bulk heatmap
# Fibro
bulk_genes <- c("PTGFR", "PLA2G2A", "RPS19", "FOS", "SFRP1", "C3", "F3", "GAS6", "SFRP2", "TMEM150C", "FBLN5", "ABCA6",
                "CD34", "HLA.DRA", "HLA.DPA1", "HLA.DPB1", "HLA.DRB1", "IFI30", "RSPO3", "PLAU", "RARRES3", "JAK2",
                "IL6", "DKK3", "CADM1", "CAPG", "AKR1C2", "COL8A2", "ITGA11", "SFRP4", "LDLRAD4", "HBEGF", "CLIC5",
                "C10orf105", "HTRA4", "PCSK6", "DEFB1", "ITGA6", "ERRFI1", "NTN4", "MET")
# Mono
bulk_genes <- c("IL1B", "NR4A2", "ATF3", "PLAUR", "CCR2", "HBEGF", "FOSB", "RGS2", "CD83", "DUSP1",
               "HTRA1", "NUPR1", "GPNMB", "ITGB5", "LGMN",
               "C1QA", "C1QB", "TMEM176B", "MARCO", "CD14", "HLA.DRA",
               "SPP1", "IFITM3", "IFI6", "LY6E", "S100A10")
# T cell
bulk_genes <- c("IL7R", "CD4", "SELL","CD40LG", "AQP3", "ANK3", "TCF7",  "CCR7", "SC5D", "NFKBIZ", "TNFRSF25", "LEF1",
               "FOXP3", "TOMM7", "IKZF2", "LAYN", "FCRL3", "STAM", "CTLA4", "TIGIT", "DUSP4",
               "CXCL13", "PDCD1", "CD200",  # "MAF", "CXCR5",
               "GZMA", "CCL5", "CCL4", "CD8A", "NKG7", "CST7", "SLAMF7", "CRTAM", "GZMK",
               "GNLY", "FGFBP2", "CX3CR1", "TGFBR3", "GZMB", "ZNF683", "SPON2", "PRF1",
               # TpH are cells: PDCD1+, CXCL13+, ICOS+, CXCR5-
               "HLA.DQA1", "HLA.DRB5", "HLA.DRA", "APOBEC3G")
# B cell
bulk_genes <- c("CD83", "CXCR4", "CD69", "IGHD", "IGHM", "BACH2", "IL4R", "IL6",
                "CD74", "MS4A1", "HLA.DPB1", "HLA.DRA", "ITGAX", "ZEB2", "CD52",
                "ACTB", "TBX21", "AICDA",
                "MZB1", "XBP1", "FKBP11", "SSR4", "DERL3")


bulk_exp <- log2tpm_fibro[bulk_genes,]
dim(bulk_exp)

mat_breaks <- seq(min(bulk_exp), max(bulk_exp), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
# mat_breaks <- quantile_breaks(both_exp, n = 11)

annotation_col <- bulk_meta_fibro[, c("Sample.ID", "Mahalanobis_20")]
colnames(annotation_col)[2] <- "disease"
rownames(annotation_col) <- bulk_meta_fibro$Sample.ID
rownames(bulk_meta_fibro) <- bulk_meta_fibro$Sample.ID
bulk_exp <- bulk_exp[, order(annotation_col$disease)]
# bulk_exp <- bulk_exp[, order(rev(p$tree_col$labels))]
scale_rows <- function(x) t(scale(t(x)))
bulk_exp <- scale_rows(bulk_exp) # Z-score
bulk_exp[bulk_exp > 2] <- 2
bulk_exp[bulk_exp < -2] <- -2

pdf("heatmap_markers_bulk_Bcell_clusCol_clusRow.pdf", width=4, height=3, onefile = FALSE, bg = "white")
pheatmap(
  mat = bulk_exp,
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(8),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = meta_colors,
  fontsize = 5,
  fontsize_row = 5,
  scale = "none"
)
dev.off()

p$tree_col$labels

