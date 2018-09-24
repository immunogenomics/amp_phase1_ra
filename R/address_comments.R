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
library(ggrepel)
source("../amp_phase1_ra_viewer/R/install-packages.R")
source("R/meta_colors.R")

meta_colors$disease = c(
  "Leukocyte-rich RA" = "#FF7F00",
  "Leukocyte-poor RA" = "#FFD8B2",
  "OA" = "#6A3D9A"
)

# -------------------------
# Load single-cell data
log2cpm <- readRDS("data/celseq_synovium_log2_5265cells_paper.rds")
meta <- readRDS("data/celseq_synovium_meta_5265cells_paper.rds")
all(colnames(log2cpm) == meta$cell_name)

# Load bulk RNA-seq data
log2tpm <- readRDS("data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)

# correct mislabeled samples
# bulk_meta$Cell.type[which(bulk_meta$Sample.ID == "S163")] <- "B cell"
# bulk_meta$Cell.type[which(bulk_meta$Sample.ID == "S164")] <- "T cell"
# bulk_meta$Cell.type[which(bulk_meta$Sample.ID == "S81")] <- "Mono"
# bulk_meta$Cell.type[which(bulk_meta$Sample.ID == "S82")] <- "Fibro"
# saveRDS(bulk_meta, "data/filtered_meta_lowinput_phase_1.rds")

# # Add manhalanobis labels
# inflam_label <- read.xls("data/postQC_all_samples.xlsx")
# inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
# inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
# inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
# inflam_label$Patient <- as.character(inflam_label$Patient)
# table(inflam_label$Mahalanobis_20)
# 
# # inter <- intersect(intersect(meta$sample, inflam_label$Patient), bulk_meta$Donor.ID)
# # inter <- intersect(inflam_label$Patient, bulk_meta$Donor.ID)
# inter <- intersect(meta$sample, inflam_label$Patient)
# meta <- meta[which(meta$sample %in% inter),]
# inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
# inflam_label <- inflam_label[, c(1,42)]
# 
# colnames(inflam_label)[1] <- "sample"
# meta <- merge(meta, inflam_label, by ="sample")
# dim(meta)
# table(meta$Mahalanobis_20)
# 
# log2cpm <- log2cpm[, which(colnames(log2cpm) %in% meta$cell_name)]
# 
# meta <- meta[order(match(meta$cell_name, colnames(log2cpm))), ]
# all(colnames(log2cpm) == meta$cell_name)


# -------------------------
# Use intersected bulk samples
bulk_meta <- bulk_meta[which(bulk_meta$Donor.ID %in% inter),]
log2tpm <- log2tpm[, which(colnames(log2tpm) %in% bulk_meta$Sample.ID)]
all(colnames(log2tpm) == bulk_meta$Sample.ID)


cell_type <- "Fibroblast"
log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_fibro) == bulk_meta_fibro$Sample.ID)
dim(bulk_meta_fibro)
dim(log2tpm_fibro)


# Add manhalanobis labels
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label$Patient <- as.character(inflam_label$Patient)
table(inflam_label$Mahalanobis_20)
inter <- intersect(bulk_meta_fibro$Donor.ID, inflam_label$Patient)


bulk_meta_fibro <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% inter),]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
bulk_meta_fibro <- bulk_meta_fibro[order(match(bulk_meta_fibro$Donor.ID, inflam_label$sample)), ]
all(bulk_meta_fibro$Donor.ID == inflam_label$sample)

bulk_meta_fibro$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(bulk_meta_fibro$Mahalanobis_20)

log2tpm_fibro <- log2tpm_fibro[, which(colnames(log2tpm_fibro) %in% bulk_meta_fibro$Sample.ID)]
log2tpm_fibro <- log2tpm_fibro[, order(match(colnames(log2tpm_fibro), bulk_meta_fibro$Sample.ID))]
all(bulk_meta_fibro$Sample.ID == colnames(log2tpm_fibro))
dim(log2tpm_fibro)

# ------------------------------------------------------------------------------------------------------------------------
# Reviewer 1 comment 1 (Option 1)
# ------------------------------------------------------------------------------------------------------------------------
# For each cell type data
# scRNA-seq
# "Fibroblast", "Monocyte", "T cell", "B cell"
type <- "Fibroblast"
log2cpm_fibro <- log2cpm[, which(meta$cell_type == type)]
meta_fibro <- meta[which(meta$cell_type == type),]
all(colnames(log2cpm_fibro) == meta_fibro$cell_name)

# For one cluster
subset <- "SC-F2"
ind <- which(meta_fibro$cluster %in% subset)
log2cpm_fibro <- log2cpm_fibro[, ind]
meta_fibro <- meta_fibro[ind, ]
# log2cpm_fibro <- log2cpm_fibro[, order(match(colnames(log2cpm_fibro), meta_fibro$cell_name))]
meta_fibro <- meta_fibro[order(meta_fibro$cell_name),]
log2cpm_fibro <- log2cpm_fibro[, order(colnames(log2cpm_fibro))]
all(colnames(log2cpm_fibro) == meta_fibro$cell_name)

# Two marker genes per cluster
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
log2tpm <- readRDS("data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)

log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == type)]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == type),]
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

inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label <- inflam_label[which(inflam_label$Patient %in% meta_1$sample),]
inflam_label$Patient <- as.character(inflam_label$Patient)
meta_1 <- meta_1[ order(match(meta_1$sample, inflam_label$Patient)), ]
all(meta_1$sample == inflam_label$Patient)
meta_1$Mahalanobis <- inflam_label$Mahalanobis_20
meta_1[1:4, ]
scaleFUN <- function(x) sprintf("%.1f", x)

ggplot() +
  geom_point(
    data = meta_1,
    mapping = aes_string(x = "gene_bulk", y = "perc_cell", fill = "Mahalanobis"),
    size = 4, stroke = 0.2, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "Bulk RNA-seq expression",
    y = paste0(type, " in ", subset,  " (%)"),
    # y = "HLA-DRA mean expression for\n SC-F2 fibroblasts by scRNA-seq",
    # y = "Percent of HLA-DRA+ cells in \nall fibroblasts by scRNA-seq",
    title = gene
  ) +
  scale_x_continuous(labels=scaleFUN) + 
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # panel.grid = element_blank(),
    plot.title = element_text(color="black", size=22, face = "italic")
  ) 

# pl <- grid.arrange(p1, p7, p5, p6, p8, p10, p11, p12, ncol=4, nrow = 2)
# pl <- grid.arrange(p22, p23, p24, p25, p26, p27, p28, p29, p31, p33, p34, p35, ncol=6, nrow = 2)
# pl <- grid.arrange(p42, p44, p45, p46, p48, p47, p51, p50, ncol = 4, nrow = 2)
pl <- grid.arrange(p61, p62, p63, p64, p65, p67, p68, p69, ncol = 4, nrow = 2)
ggsave(file = paste("bulk_perc_fibro", ".pdf", sep = ""), pl,
        width = 12, height = 6.5, dpi = 300)
dev.off()


# ------------------------------------------------------------------------------------------------------------------------
# Reviewer 1 comment 1 (Option 2) Percent of fibroblasts have >0 marker; ignore clusters
# ------------------------------------------------------------------------------------------------------------------------
# For each cell type data
# scRNA-seq
# "Fibroblast", "Monocyte", "T cell", "B cell"
type <- "T cell"
log2cpm_fibro <- log2cpm[, which(meta$cell_type == type)]
meta_fibro <- meta[which(meta$cell_type == type),]
all(colnames(log2cpm_fibro) == meta_fibro$cell_name)

# Two marker genes per cluster
gene <- "HLA-DRB1"
meta_fibro$gene <- as.numeric(log2cpm_fibro[which(rownames(log2cpm_fibro) == gene),])
dat_percent <- meta_fibro %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(percent = sum(gene > 0) / length(gene) * 100,
                   ave_mean = mean(gene),
                   ave_median = median(gene))
dat_percent <- as.data.frame(dat_percent)
dat_percent <- dat_percent[order(dat_percent$sample),]


# Load bulk data
log2tpm <- readRDS("data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)
log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == type)]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == type),]
all(colnames(log2tpm_fibro) == bulk_meta_fibro$cell_name)
gene <- "HLA.DRB1"
bulk_meta_fibro$gene_bulk <- as.numeric(log2tpm_fibro[which(rownames(log2tpm_fibro) == gene),]) 

# Intersect
inter <- intersect(meta_fibro$sample, bulk_meta_fibro$Donor.ID)
meta_1 <- dat_percent[which(dat_percent$sample %in% inter),]
meta_2 <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% inter),]
meta_1$sample <- as.character(meta_1$sample)
meta_1 <- meta_1[ order(match(meta_1$sample, meta_2$Donor.ID)), ]
all(meta_1$sample == meta_2$Donor.ID)
meta_1$gene_bulk <- meta_2$gene_bulk

# Add manhalonobis labels 
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label <- inflam_label[which(inflam_label$Patient %in% meta_1$sample),]
inflam_label$Patient <- as.character(inflam_label$Patient)
meta_1 <- meta_1[ order(match(meta_1$sample, inflam_label$Patient)), ]
all(meta_1$sample == inflam_label$Patient)
meta_1$Mahalanobis <- inflam_label$Mahalanobis_20
scaleFUN <- function(x) sprintf("%.1f", x)

fit <- lm(gene_bulk ~ percent, data = meta_1)
# summary(fit)
# summary(fit)$adj.r.squared
# summary(fit)$r.squared
# summary(fit)$coefficients[2,"Pr(>|t|)"]
fit$model$Mahalanobis <- meta_1$Mahalanobis
p42 <- ggplot(
  fit$model,
  aes(
    x=gene_bulk,
    y=percent,
    fill = Mahalanobis
    )
  ) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 1.5, linetype="dashed",
              col= "darkgrey", fill="lightgrey") +
  geom_point(
    shape=21, size = 4.5, stroke = 0.35
  ) +
  scale_fill_manual(values = meta_colors$Case.Control) +
  labs(
    x = "",
    y = "",
    title = gene
    # x = "Bulk RNA-seq expression",
    # y = "Percent of non-zero expressing cells",
    # title = paste0(gene, " in ", type)
  ) +
  scale_x_continuous(labels=scaleFUN) + 
  theme_bw(base_size = 25) +
  theme(
    legend.position = "none",
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(color="black", size=25)
  ) 


pf <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 4, nrow = 2)
ggsave(file = paste("bulk_perc_fibro_option2", ".pdf", sep = ""), pf,
       width = 15.5, height = 8, dpi = 300)

pm <- grid.arrange(p11, p12, p13, p14, p15, p16, p17, p18, ncol = 4, nrow = 2)
ggsave(file = paste("bulk_perc_mono_option2", ".pdf", sep = ""), pm,
       width = 15.5, height = 8, dpi = 300)

pb <- grid.arrange(p21, p22, p23, p24, p25, p26, p27, p28, ncol = 4, nrow = 2)
ggsave(file = paste("bulk_perc_bcell_option2", ".pdf", sep = ""), pb,
       width = 15.5, height = 8, dpi = 300)

pt <- grid.arrange(p31, p32, p33, p34, p35, p36, p37, p38, p39, p40, p41, p42, ncol = 4, nrow = 3)
ggsave(file = paste("bulk_perc_tcell_option2", ".pdf", sep = ""), pt,
       width = 15.5, height = 12, dpi = 300)


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
bulk_genes_fibro <- c("PTGFR", "PLA2G2A", "RPS19", "FOS", "SFRP1", "C3", "F3", "GAS6", "SFRP2", "TMEM150C", "FBLN5", "ABCA6",
                "CD34", "HLA.DRA", "HLA.DPA1", "HLA.DPB1", "HLA.DRB1", "IFI30", "RSPO3", "PLAU", "RARRES3", "JAK2",
                "IL6", "DKK3", "CADM1", "CAPG", "AKR1C2", "COL8A2", "ITGA11", "SFRP4", "LDLRAD4", "HBEGF", "CLIC5",
                "C10orf105", "HTRA4", "PCSK6", "DEFB1", "ITGA6", "DNASE1L3", "ERRFI1", "NTN4", "MET")
# Mono
bulk_genes_mono <- c("IL1B", "NR4A2", "ATF3", "PLAUR", "CCR2", "HBEGF", "FOSB", "RGS2", "CD83", "DUSP1",
               "HTRA1", "NUPR1", "GPNMB", "ITGB5", "LGMN",
               "C1QA", "C1QB", "TMEM176B", "MARCO", "CD14", "HLA.DRA",
               "SPP1", "IFITM3", "IFI6", "LY6E", "S100A10")
# T cell
bulk_genes_tcell <- c("IL7R", "CD4", "SELL","CD40LG", "AQP3", "ANK3", "TCF7",  "CCR7", "SC5D", "NFKBIZ", "TNFRSF25", "LEF1",
               "FOXP3", "TOMM7", "IKZF2", "LAYN", "FCRL3", "STAM", "CTLA4", "TIGIT", "DUSP4",
               "CXCL13", "PDCD1", "CD200",  # "MAF", "CXCR5",
               "GZMA", "CCL5", "CCL4", "CD8A", "NKG7", "CST7", "SLAMF7", "CRTAM", "GZMK",
               "GNLY", "FGFBP2", "CX3CR1", "TGFBR3", "GZMB", "ZNF683", "SPON2", "PRF1",
               # TpH are cells: PDCD1+, CXCL13+, ICOS+, CXCR5-
               "HLA.DQA1", "HLA.DRB5", "HLA.DRA", "APOBEC3G")
# B cell
bulk_genes_bcell <- c("CD83", "CXCR4", "CD69", "IGHD", "IGHM", "BACH2", "IL4R", "IL6",
                "CD74", "MS4A1", "HLA.DPB1", "HLA.DRA", "ITGAX", "ZEB2", "CD52",
                "ACTB", "TBX21", "AICDA",
                "MZB1", "XBP1", "FKBP11", "SSR4", "DERL3")


bulk_exp <- log2tpm_fibro[bulk_genes_fibro,]
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

pdf("heatmap_markers_bulk_fibro_allsamples.pdf", width=5, height=4, onefile = FALSE, bg = "white")
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

# -------
# Reviewer 1 comment 7
samp_clu <- table(meta$sample, meta$cluster)
samp_clu <- as.data.frame(samp_clu)
colnames(samp_clu) <- c("sample", "cluster", "freq")
samp_clu$sample <- as.character(samp_clu$sample)

non_samples <- names(which(table(meta$sample) == 0))
samp_clu <- samp_clu[-which(samp_clu$sample %in% non_samples),]
table(samp_clu$sample)
samp_clu$cluster = factor(samp_clu$cluster, levels=c("SC-F1", "SC-F2", "SC-F3", "SC-F4", "SC-M1", "SC-M2", "SC-M3", "SC-M4",
                                                     "SC-T1", "SC-T2", "SC-T3", "SC-T4", "SC-T5", "SC-T6", "SC-B1", "SC-B2", 
                                                     "SC-B3", "SC-B4"))
samp_clu$sample = factor(samp_clu$sample, levels = c( "300-0122" ,"300-0153", "300-0211", "300-0213", "300-0481", "300-0482",
                                                      "300-0483", "300-0485", "300-0486", "300-0487", "300-0511","300-0528",
                                                      "300-0546","300-2590", "301-0122", "301-0244","301-0151","301-0153",
                                                      "301-0155","301-0121", # lek-poor
                                                      "301-0163", # lek-poor
                                                      "301-0250", # lek-poor
                                                      "301-0132", # OA
                                                      "301-0159", # OA, 
                                                      "301-0161" # OA
                                                      )
                         )

# colfunc <- colorRampPalette(c("yellow", "orange", "red4"))
# plot(rep(1,15),col=colfunc(15),pch=19,cex=3)
# 
# colfunc <- colorRampPalette(c("lightgoldenrodyellow", "lightgoldenrod3"))
# plot(rep(1,3),col=colfunc(3),pch=19,cex=3)
# 
# colfunc <- colorRampPalette(c("plum1", "purple4"))
# plot(rep(1,3),col=colfunc(3),pch=19,cex=3)

ggplot(
  data=samp_clu,
  aes(x=cluster, y= freq, fill = sample)
  ) +
  geom_bar(stat="identity",
           position = "fill",
           # position = "stack",
           width = 0.8
  ) +
  # facet_grid(cluster ~ ., scales = "free", space = "free_x") +
  scale_fill_manual("", values = meta_colors$sample) +
  labs(
    x = "scRNA-seq clusters",
    y = "Proportion of cells per donor"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic(base_size = 12) +
  theme(    
    # legend.position="none",
    # axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x=element_text(angle=30, hjust=1),
    axis.text.y = element_text(size = 12)
    )
ggsave(file = paste("barplot_sc_cluster_per_patient_v2", ".pdf", sep = ""),
       width = 9, height = 6, dpi = 100)
dev.off()

# # Reviewer 1 comment 7
# sc
log2cpm <- as.data.frame(log2cpm)
markers <- gsub("[.]", "-", Reduce(union, list(bulk_genes_fibro, bulk_genes_mono, bulk_genes_tcell, bulk_genes_bcell)))
log2cpm_marker <- log2cpm[markers,]
all(colnames(log2cpm_marker) == meta$cell_name)
log2cpm_marker <- as.data.frame(t(log2cpm_marker))
log2cpm_marker$cluster <- meta$cluster

log2cpm_marker_ave <- log2cpm_marker %>%
                      group_by(cluster) %>%
                      summarise_at(vars(-cluster), funs(mean(., na.rm=TRUE)))
log2cpm_marker_ave <- as.data.frame(log2cpm_marker_ave)
rownames(log2cpm_marker_ave) <- log2cpm_marker_ave$cluster
log2cpm_marker_ave <- log2cpm_marker_ave[, -1]

# bulk
log2tpm_marker <- log2tpm[Reduce(union, list(bulk_genes_fibro, bulk_genes_mono, bulk_genes_tcell, bulk_genes_bcell)),]
log2tpm_marker <- t(log2tpm_marker)
colnames(log2tpm_marker) <- colnames(log2cpm_marker_ave)

# Combine
exp_marker <- rbind.data.frame(log2cpm_marker_ave, log2tpm_marker)
dim(exp_marker)

# Spearman correlation
cor_mat <- cor(t(exp_marker), method = "spearman")
cor_plot <- cor_mat[c(1:18), c(19:ncol(cor_mat))]
dim(cor_plot)
cor_plot <- cor_plot[c("SC-T1", "SC-T2", "SC-T3", "SC-T4", "SC-T5", "SC-T6", "SC-B1", "SC-B2", "SC-B3", "SC-B4", 
                       "SC-F1", "SC-F2", "SC-F3", "SC-F4", "SC-M1", "SC-M2", "SC-M3", "SC-M4"), ]

all(bulk_meta$Sample.ID == colnames(cor_plot))
annotation_col <- bulk_meta[, c("disease", "cluster")] 
rownames(annotation_col) <- ra_sle_meta$name
rownames(ra_sle_meta) <- ra_sle_meta$name

pheatmap(
  cor_plot,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = rev(colorRampPalette(brewer.pal(n=9, name = "RdYlBu"))(9)),
  fontsize_row = 10,
  scale = "none",
  legend = TRUE
  )
dev.copy(png,file = paste("heatmp_bulksamples_scclusters", ".png", sep = ""), width=3,height=4,units="in",res=300)
dev.off()

# ------
# Reviewer1 comment 6: Plot all viable cells vs CD45
# facs <- read.xls("data/171215_FACS_data_for_figure.xlsx")
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
facs <- read.csv("data/170829 AMP cell count.csv")
dim(facs)

# Plot the 21 scRNA-seq samples
meta_flow <- facs[which(facs$Sample %in% meta$sample),]
# Plot all the 51 samples
# meta_flow <- facs[which(facs$Sample %in% inflam_label$Patient),]
dim(meta_flow)
meta_flow$CD45 <- meta_flow$Tcell + meta_flow$Bcell
meta_flow$CD45_perc <- meta_flow$CD45 / meta_flow$All.viable.cells
meta_flow$leuk_perc <- (meta_flow$Tcell + meta_flow$Bcell + meta_flow$Monocyte) / meta_flow$All.viable.cells
meta_flow$All.viable.cells_2 <- meta_flow$All.viable.cells - meta_flow$Other.cells
meta_flow$CD45_perc_2 <- meta_flow$CD45 / meta_flow$All.viable.cells_2
meta_flow$leuk_perc_2 <- (meta_flow$Tcell + meta_flow$Bcell + meta_flow$Monocyte) / meta_flow$All.viable.cells_2


inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
table(inflam_label$Mahalanobis_20)

inflam_label <- inflam_label[which(inflam_label$Patient %in% meta_flow$Sample),]
inflam_label$Patient <- as.character(inflam_label$Patient)
inflam_label <- inflam_label[order(match(inflam_label$Patient, meta_flow$Sample)), ]
all(meta_flow$Sample == inflam_label$Patient)
meta_flow$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(meta_flow$Mahalanobis_20)
meta_flow$Sample <- as.character(meta_flow$Sample)

ggplot() +
  geom_point(
    data = meta_flow, 
    mapping = aes_string(x = "All.viable.cells", y = "CD45", fill = "Mahalanobis_20"),
    size = 2, stroke = 0.2, shape = 21
  ) +
  geom_text_repel(
    data = meta_flow[which(meta_flow$CD45 > 50000),],
    aes_string(x = "All.viable.cells", y = "CD45", label = "Sample"),
    size = 1.5, color = "black",
    box.padding = unit(0.15, "lines"),
    point.padding = unit(0.1, "lines")
  ) +
  scale_fill_manual("", values = meta_colors$disease) +
  scale_x_continuous(labels = scales::comma) + 
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "All viable cells",
    y = "Cell abundance of CD45",
    title = "Flow cytometry"
  ) +
  theme_classic(base_size = 7) +
  theme(
    # legend.position="none",
    # axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(size = 7),
    axis.text.x=element_text(angle= 30, hjust=1),
    axis.text.y = element_text(size = 7)
  )
ggsave(file = paste("allviablecells_vs_CD45_flow_cytometry", ".png", sep = ""),
       width = 3.5, height = 2, dpi = 300)
dev.off()


# Plot flow barplot ordering thes same way of Figure 3B
meta_flow <- meta_flow[, -7]
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
table(inflam_label$Mahalanobis_20)

inflam_label <- inflam_label[which(inflam_label$Patient %in% meta_flow$Sample),]
inflam_label$Patient <- as.character(inflam_label$Patient)
inflam_label <- inflam_label[order(match(inflam_label$Patient, meta_flow$Sample)), ]
all(meta_flow$Sample == inflam_label$Patient)
meta_flow$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(meta_flow$Mahalanobis_20)

t <- meta_flow[, c(1,2,7)]
t$cell_type <- rep("T cell", nrow(t))
colnames(t)[2] <- "freq"

b <- meta_flow[, c(1,3,7)]
b$cell_type <- rep("B cell", nrow(b))
colnames(b)[2] <- "freq"

m <- meta_flow[, c(1,4,7)]
m$cell_type <- rep("Monocyte", nrow(m))
colnames(m)[2] <- "freq"

f <- meta_flow[, c(1,5,7)]
f$cell_type <- rep("Fibroblast", nrow(f))
colnames(f)[2] <- "freq"

o <- meta_flow[, c(1,6,7)]
o$cell_type <- rep("Other cell", nrow(o))
colnames(o)[2] <- "freq"

flow_plot <- rbind.data.frame(t, b, m, f, o)
flow_plot[1:4,]
flow_plot$Sample <- as.character(flow_plot$Sample)
flow_plot$Mahalanobis_order = factor(flow_plot$Mahalanobis_20, levels=c('OA','Leukocyte-poor RA','Leukocyte-rich RA'))
flow_plot$Sample <- factor(flow_plot$Sample, 
                           levels = c("301-0161", "301-0159", "301-0132", "301-0121",
                                      "301-0250", "301-0163", "300-0481", "300-0153",
                                      "300-0485", "300-0483", "300-0211", "300-0213",
                                      "301-0122", "300-0528", "300-0122", "300-2590", 
                                      "300-0486", "300-0546", "300-0482", "300-0487", "300-0511"))

# Remove other cells
flow_plot <- flow_plot[-which(flow_plot$cell_type == "Other cell"),]

ggplot(
  data=flow_plot,
  aes(x=Sample, y= freq, fill = cell_type)
  ) +
  geom_bar(stat="identity",
           # position = "fill",
           position = "stack",
           width = 0.85
  ) +
  facet_grid(. ~ Mahalanobis_order, scales = "free", space = "free_x") +
  scale_fill_manual("", values = meta_colors$type) +
  labs(
    x = "Donors",
    y = "Number of cells",
    title = "Flow cytometry"
  ) +
  theme_bw(base_size = 12) +
  theme(    
    # legend.position="none",
    # axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x=element_text(angle=30, hjust=1),
    axis.text.y = element_text(size = 12)
  )
ggsave(file = paste("barplot_per_patient_flow_cytometry", ".pdf", sep = ""),
       width = 8.5, height = 6, dpi = 200)
dev.off()


# ----------------------------
# Reviewer 2 comment 5.4
# ----------------------------
# Correlation between cytof and bulk on the measured proteomic cytof markers

load("data/synData.Fibro.downsample.SNE.RData")
load("data/synData.Bcell.downsample.SNE.RData")
load("data/synData.Mono.downsample.SNE.RData")
load("data/synData.Tcell.downsample.SNE.RData")

# 300-0486C, 300-0511A
synData.Fibro.downsample$sampleID[which(synData.Fibro.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Fibro.downsample$sampleID[which(synData.Fibro.downsample$sampleID == "300-0511A")] <- "300-0511"
synData.Bcell.downsample$sampleID[which(synData.Bcell.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Bcell.downsample$sampleID[which(synData.Bcell.downsample$sampleID == "300-0511A")] <- "300-0511"
synData.Mono.downsample$sampleID[which(synData.Mono.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Mono.downsample$sampleID[which(synData.Mono.downsample$sampleID == "300-0511A")] <- "300-0511"
synData.Tcell.downsample$sampleID[which(synData.Tcell.downsample$sampleID == "300-0486C")] <- "300-0486"
synData.Tcell.downsample$sampleID[which(synData.Tcell.downsample$sampleID == "300-0511A")] <- "300-0511"


# ggplot() +
#   geom_point(
#     data = synData.Bcell.downsample,
#     mapping = aes_string(x = "SNE1", y = "SNE2", fill = "CD19"),
#     shape = 21, size = 1.1, stroke = 0.011
#   ) +
#   scale_fill_gradientn(
#     colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9),
#     name = "Normalized intensity"
#   ) +
#   labs(
#     x = NULL,
#     y = NULL
#   ) +
#   theme_bw(base_size = 22) +
#   theme(
#     legend.position = "none",
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.grid = element_blank()
#   )

# fibro
synData.Fibro_sampleID <- synData.Fibro.downsample[, c("sampleID", "CD90")]
synData.Fibro_mean <- synData.Fibro_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean = mean(CD90),
                   ave_median = median(CD90))
synData.Fibro_mean <- as.data.frame(synData.Fibro_mean)

# B cell
# synData.Bcell_sampleID <- synData.Bcell.downsample[, c("sampleID", "IgD")]
synData.Bcell_sampleID <- synData.Bcell.downsample[, c("sampleID", "CD11c")]
synData.Bcell_mean <- synData.Bcell_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean = mean(CD11c),
                   ave_median = median(CD11c))
synData.Bcell_mean <- as.data.frame(synData.Bcell_mean)

# Mono
# synData.Mono_sampleID <- synData.Mono.downsample[, c("sampleID", "CCR2")]
synData.Mono_sampleID <- synData.Mono.downsample[, c("sampleID", "CD38")]
synData.Mono_mean <- synData.Mono_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean = mean(CD38),
                   ave_median = median(CD38))
synData.Mono_mean <- as.data.frame(synData.Mono_mean)

# T cell
synData.Tcell_sampleID <- synData.Tcell.downsample[, c("sampleID", "CD8a")]
synData.Tcell_mean <- synData.Tcell_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean = mean(CD8a),
                   ave_median = median(CD8a))
synData.Tcell_mean <- as.data.frame(synData.Tcell_mean)

# bulk data
# fibro
bulk_meta$CD90_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "THY1"), ])
cell_type <- "Fibroblast"
log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_fibro) == bulk_meta_fibro$Sample.ID)
dim(bulk_meta_fibro)
dim(log2tpm_fibro)

# B cell
bulk_meta$MS4A1 <- as.numeric(log2tpm[which(rownames(log2tpm) == "MS4A1"), ])
bulk_meta$ITGAX <- as.numeric(log2tpm[which(rownames(log2tpm) == "ITGAX"), ])
bulk_meta$HLA.DRA_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "HLA.DRA"), ])
bulk_meta$CD19_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CD19"), ])
bulk_meta$IGHM_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "IGHM"), ])
bulk_meta$IGHD_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "IGHD"), ])
cell_type <- "B cell"
log2tpm_bcell <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_bcell <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_bcell) == bulk_meta_bcell$Sample.ID)
dim(bulk_meta_bcell)
dim(log2tpm_bcell)

# Mono
bulk_meta$CD38_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CD38"), ])
cell_type <- "Monocyte"
log2tpm_mono <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_mono <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
bulk_meta_mono <- bulk_meta_mono[-duplicated(bulk_meta_mono$Donor.ID),]
log2tpm_mono <- log2tpm_mono[, which(colnames(log2tpm_mono) %in% bulk_meta_mono$Sample.ID)]
all(colnames(log2tpm_mono) == bulk_meta_mono$Sample.ID)
dim(bulk_meta_mono)
dim(log2tpm_mono)

# T cell
bulk_meta$PDCD1_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "PDCD1"), ])
bulk_meta$CD8A_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CD8A"), ])
cell_type <- "T cell"
log2tpm_tcell <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_tcell <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_tcell) == bulk_meta_tcell$Sample.ID)
dim(bulk_meta_tcell)
dim(log2tpm_tcell)

# Get manhalonobis labels
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label$Patient <- as.character(inflam_label$Patient)
table(inflam_label$Mahalanobis_20)

# Merge manhalonobis labels to bulk meta
# Fibro
inter <- intersect(bulk_meta_fibro$Donor.ID, inflam_label$Patient)
bulk_meta_fibro <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% inter),]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
bulk_meta_fibro <- bulk_meta_fibro[order(match(bulk_meta_fibro$Donor.ID, inflam_label$Patient)), ]
all(bulk_meta_fibro$Donor.ID == inflam_label$Patient)
bulk_meta_fibro$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(bulk_meta_fibro$Mahalanobis_20)

# B cell
inter <- intersect(bulk_meta_bcell$Donor.ID, inflam_label$Patient)
bulk_meta_bcell <- bulk_meta_bcell[which(bulk_meta_bcell$Donor.ID %in% inter),]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
bulk_meta_bcell <- bulk_meta_bcell[order(match(bulk_meta_bcell$Donor.ID, inflam_label$Patient)), ]
all(bulk_meta_bcell$Donor.ID == inflam_label$Patient)
bulk_meta_bcell$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(bulk_meta_bcell$Mahalanobis_20)

# Mono
inter <- intersect(bulk_meta_mono$Donor.ID, inflam_label$Patient)
bulk_meta_mono <- bulk_meta_mono[which(bulk_meta_mono$Donor.ID %in% inter),]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
inflam_label <-  inflam_label[order(match(inflam_label$Patient, bulk_meta_mono$Donor.ID)),]
all(bulk_meta_mono$Donor.ID == inflam_label$Patient)
bulk_meta_mono$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(bulk_meta_mono$Mahalanobis_20)

# T cell
inter <- intersect(bulk_meta_tcell$Donor.ID, inflam_label$Patient)
bulk_meta_tcell <- bulk_meta_tcell[which(bulk_meta_tcell$Donor.ID %in% inter),]
inflam_label <- inflam_label[which(inflam_label$Patient %in% inter),]
inflam_label <-  inflam_label[order(match(inflam_label$Patient, bulk_meta_tcell$Donor.ID)),]
all(bulk_meta_tcell$Donor.ID == inflam_label$Patient)
bulk_meta_tcell$Mahalanobis_20 <- inflam_label$Mahalanobis_20
table(bulk_meta_tcell$Mahalanobis_20)


# Overlapped samples between bulk fibroblast samples with cytof fibroblast samples
# Fibro
over <- intersect(synData.Fibro_mean$sampleID, bulk_meta_fibro$Donor.ID)
synData.Fibro_mean <- synData.Fibro_mean[which(synData.Fibro_mean$sampleID %in% over), ]
bulk_meta_fibro <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% over), ]
bulk_meta_fibro <- bulk_meta_fibro[order(match(bulk_meta_fibro$Donor.ID, synData.Fibro_mean$sampleID)), ]
all(synData.Fibro_mean$sampleID == bulk_meta_fibro$Donor.ID)
bulk_meta_fibro$cytof_ave_mean <- synData.Fibro_mean$ave_mean
bulk_meta_fibro$cytof_ave_median <- synData.Fibro_mean$ave_median

# B cell
over <- intersect(synData.Bcell_mean$sampleID, bulk_meta_bcell$Donor.ID)
synData.Bcell_mean <- synData.Bcell_mean[which(synData.Bcell_mean$sampleID %in% over), ]
bulk_meta_bcell <- bulk_meta_bcell[which(bulk_meta_bcell$Donor.ID %in% over), ]
bulk_meta_bcell <- bulk_meta_bcell[order(match(bulk_meta_bcell$Donor.ID, synData.Bcell_mean$sampleID)), ]
all(synData.Bcell_mean$sampleID == bulk_meta_bcell$Donor.ID)
bulk_meta_bcell$cytof_ave_mean <- synData.Bcell_mean$ave_mean
bulk_meta_bcell$cytof_ave_median <- synData.Bcell_mean$ave_median

# Mono
over <- intersect(synData.Mono_mean$sampleID, bulk_meta_mono$Donor.ID)
synData.Mono_mean <- synData.Mono_mean[which(synData.Mono_mean$sampleID %in% over), ]
bulk_meta_mono <- bulk_meta_mono[which(bulk_meta_mono$Donor.ID %in% over), ]
bulk_meta_mono <- bulk_meta_mono[order(match(bulk_meta_mono$Donor.ID, synData.Mono_mean$sampleID)), ]
all(synData.Mono_mean$sampleID == bulk_meta_mono$Donor.ID)
bulk_meta_mono$cytof_ave_mean <- synData.Mono_mean$ave_mean
bulk_meta_mono$cytof_ave_median <- synData.Mono_mean$ave_median

# T cell
over <- intersect(synData.Tcell_mean$sampleID, bulk_meta_tcell$Donor.ID)
synData.Tcell_mean <- synData.Tcell_mean[which(synData.Tcell_mean$sampleID %in% over), ]
bulk_meta_tcell <- bulk_meta_tcell[which(bulk_meta_tcell$Donor.ID %in% over), ]
bulk_meta_tcell <- bulk_meta_tcell[order(match(bulk_meta_tcell$Donor.ID, synData.Tcell_mean$sampleID)), ]
all(synData.Tcell_mean$sampleID == bulk_meta_tcell$Donor.ID)
bulk_meta_tcell$cytof_ave_mean <- synData.Tcell_mean$ave_mean
bulk_meta_tcell$cytof_ave_median <- synData.Tcell_mean$ave_median


# Fibro
ggplot() +
  geom_point(
    data = bulk_meta_fibro,
    # mapping = aes_string(x = "CD90", y = "cytof_ave_mean", fill = "Mahalanobis_20"),
    mapping = aes_string(x = "CD34_gene", y = "cytof_ave_mean"),
    size = 3, stroke = 0.2, shape = 20, color = "#08519C"
  ) +
  # scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "Bulk RNA-seq expression",
    y = "Averaged mass cytometry \n proteomic expression",
    title = "Fibroblast samples\nCD34"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CD34", ".pdf", sep = ""),
       width = 3, height = 3, dpi = 300)
dev.off()

# B cell
ggplot() +
  geom_point(
    data = bulk_meta_bcell,
    mapping = aes_string(x = "IGHD_gene", y = "cytof_ave_mean"),
    size = 3, stroke = 0.2, shape = 20, color = "#E31A1C"
  ) +
  # scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "Bulk RNA-seq expression",
    y = "Averaged mass cytometry \n proteomic expression",
    title = "B cell samples\nIgD"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_IGH", ".pdf", sep = ""),
       width = 3, height = 3, dpi = 300)
dev.off()

# Mono
ggplot() +
  geom_point(
    data = bulk_meta_mono,
    mapping = aes_string(x = "CCR2_gene", y = "cytof_ave_mean"),
    size = 3, stroke = 0.2, shape = 20, color = "#DE77AE"
  ) +
  # scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "Bulk RNA-seq expression",
    y = "Averaged mass cytometry \n proteomic expression",
    title = "Monocyte samples\nCCR2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CCR2", ".pdf", sep = ""),
       width = 3, height = 3, dpi = 300)
dev.off()


# T cell
ggplot() +
  geom_point(
    data = bulk_meta_tcell,
    mapping = aes_string(x = "PDCD1_gene", y = "cytof_ave_mean"),
    size = 3, stroke = 0.2, shape = 20, color = "#A65628"
  ) +
  # scale_fill_manual(values = meta_colors$disease, name = "") +
  labs(
    x = "Bulk RNA-seq expression",
    y = "Averaged mass cytometry \n proteomic expression",
    title = "T cell samples\nPD1"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_PDCD1", ".pdf", sep = ""),
       width = 3.2, height = 3, dpi = 300)
dev.off()


# Plot fitted line -------
# Fibro
fit <- lm(CD90_gene ~ cytof_ave_mean, data = bulk_meta_fibro)
ggplot(
  fit$model,
  aes(
    x=CD90_gene,
    y=cytof_ave_mean
   )
  ) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 1.5, linetype="dashed",
              col= "darkgrey", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#08519C"
  ) +
  labs(
    x = "",
    y = "",
    title = "Fibroblast samples\nCD90"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CD90_fit", ".pdf", sep = ""),
       width = 3, height = 3.2, dpi = 300)
dev.off()

# Mono
fit <- lm(CD38_gene ~ cytof_ave_mean, data = bulk_meta_mono)
ggplot(
  fit$model,
  aes(
    x=CD38_gene,
    y=cytof_ave_mean
    )
  ) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 1.5, linetype="dashed",
              col= "darkgrey", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#DE77AE"
  ) +
  labs(
    x = "",
    y = "",
    title = "Monocyte samples\nCD38"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CD38_fit", ".pdf", sep = ""),
       width = 3, height = 3.2, dpi = 300)
dev.off()


# B cell
fit <- lm(ITGAX ~ cytof_ave_mean, data = bulk_meta_bcell)
ggplot(
  fit$model,
  aes(
    x=ITGAX,
    y=cytof_ave_mean
   )
  ) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 1.5, linetype="dashed",
              col= "darkgrey", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#E31A1C"
  ) +
  labs(
    x = "",
    y = "",
    title = "B cell samples\nCD11c"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_ITGAX_fit", ".pdf", sep = ""),
       width = 3, height = 3.2, dpi = 300)
dev.off()


# T cell
fit <- lm(CD8A_gene ~ cytof_ave_mean, data = bulk_meta_tcell)
ggplot(
  fit$model,
  aes(
    x=CD8A_gene,
    y=cytof_ave_mean
    )
  ) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 1.5, linetype="dashed",
              col= "darkgrey", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#A65628"
  ) +
  labs(
    x = "",
    y = "",
    title = "T cell samples\nCD8"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CD8A_fit", ".pdf", sep = ""),
       width = 3, height = 3.2, dpi = 300)
dev.off()


# ----------------------------
# Correlation between cytof (downsample the same number of cells per donor with scRNA-seq) and scRNA-sq 
type <- "Fibroblast"
log2cpm_fibro <- as.data.frame(log2cpm[, which(meta$cell_type == type)])
meta_fibro <- meta[which(meta$cell_type == type),]
# log2cpm_fibro <- log2cpm_fibro[,order(match(colnames(log2cpm_fibro), meta_fibro$cell_name))]
all(colnames(log2cpm_fibro) == meta_fibro$cell_name)
table(meta_fibro$sample)
meta_fibro$CD90 <- as.numeric(log2cpm_fibro[which(rownames(log2cpm_fibro) == "THY1"), ])

# overlapped samples between scRNA-seq and cytof: only 7 in fibroblast 
over <- intersect(synData.Fibro.downsample$sampleID, meta_fibro$sample)

sc_fibro_over <- meta_fibro[which(meta_fibro$sample %in% over), c("sample", "Mahalanobis_20", "CD90")]
table(sc_fibro_over$sample)

cytof_fibro_over <- synData.Fibro.downsample[which(synData.Fibro.downsample$sampleID %in% over), c("CD90", "sampleID")]
table(cytof_fibro_over$sampleID)

sample(x, size, replace = FALSE, prob = NULL)


# ----------------------------
# Load Proteomic markers for each cell
meta <- meta[order(meta$sample), ]
# facs1 <- read.table(file = 'data/Phase1_facs.tsv', sep = '\t', header = TRUE)
facs2 <- read.table(file = 'data/celseq_flow.tsv', sep = '\t', header = TRUE) 
facs2 <- facs2[order(facs2$sample),]

facs2$test <- paste(facs2$sample, "_", facs2$plate, "_", facs2$well384, sep="")
meta$test <- paste(meta$sample, "_", substr(meta$cell_name, 1, 4), "_", substr(meta$cell_name, 11, 13), sep="")
facs2_protein <- facs2[which(facs2$test %in% meta$test),]
dim(facs2_protein)
meta_protein <- meta[which(meta$test %in% facs2_protein$test),]
meta_protein <- meta_protein[order(match(meta_protein$test, facs2_protein$test)), ]
all(facs2_protein$test == meta_protein$test)
meta_protein <- merge(meta_protein, facs2_protein[, c(10:27)], by = "test")

log2cpm_protein <- log2cpm[,which(colnames(log2cpm) %in% meta_protein$cell_name)]
meta_protein <- meta_protein[order(match(meta_protein$cell_name, colnames(log2cpm))),] 
all(meta_protein$cell_name == colnames(log2cpm_protein))
meta_protein$CD3D_gene <- log2cpm_protein[which(rownames(log2cpm_protein) == "CD3D"),]
# Use meta_protein and log2cpm_protein to compare mRNA and protein expression
# meta_protein$CD3_log <- log10(meta_protein$CD3)
# is.na(meta_protein$CD3_log) <- 0

# ggplot() +
#   geom_point(
#     data=meta_protein,
#     aes_string(x = "CD3D_gene", y = "CD3", fill = "cell_type"),
#     size = 4, stroke = 0.2, shape = 21
#   ) 

all_cells <- readRDS("data/all_cells_fine_cluster_label.rds")
all_cells <- all_cells[which(rownames(all_cells) %in% meta_protein$cell_name),]
all_cells <- all_cells[order(match(rownames(all_cells), meta_protein$cell_name)),] 
all(meta_protein$cell_name == rownames(all_cells))
meta_protein$T1 <- all_cells$T1
meta_protein$T2 <- all_cells$T2

# Remove the extreme values for some rare dots. Or we can take the log10
meta_protein$CD3[is.na(meta_protein$CD3)] <- 20000
meta_protein$CD3[which(meta_protein$CD3 > 30000)] <- 20000
meta_protein$CD14[is.na(meta_protein$CD14)] <- 20000
meta_protein$CD19[is.na(meta_protein$CD19)] <- 20000
meta_protein$CD19[which(meta_protein$CD19 > 40000)] <- 20000
meta_protein$PDPN[which(meta_protein$PDPN > 40000)] <- 20000
meta_protein$CD14[which(meta_protein$CD14 > 40000)] <- 20000
meta_protein$CD90[which(meta_protein$CD90 > 100000)] <- 50000
meta_protein$CD34[which(meta_protein$CD34 > 50000)] <- 40000
meta_protein$CD27[which(meta_protein$CD27 > 10000)] <- 7000
meta_protein$CD45[which(meta_protein$CD45 > 50000)] <- 2000

protein <- "CD45"
ggplot() +
  geom_point(
    data = meta_protein[order(meta_protein$CD45),],
    mapping = aes_string(x = "T1", y = "T2", fill = protein),
    size = 1.2, stroke = 0.001, shape = 21
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(10),
    name = "Protein fluorescence"
  ) +
  labs(
    x = "tSNE1",
    y = "tSNE2",
    title = protein
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    # panel.grid = element_blank(),
    plot.title = element_text(color="black")
  ) 
ggsave(file = paste("protein_", protein, ".png", sep = ""),
       width = 3.5, height = 3, dpi = 300)
dev.off()
