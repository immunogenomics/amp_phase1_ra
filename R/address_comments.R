# Adress reviewers' comments for AMP Phase I RA paper
# Fan Zhang
# 2018-08-14

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/")

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
library(cetcolor)
library(gridExtra)
require(viridis)
require(ggbeeswarm)
require(scales)
require(reshape2)
require(ggrepel)

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
# Change HLA.DRA to HLA-DRA
temp <- rownames(log2tpm)[grep("HLA.", rownames(log2tpm))]
rownames(log2tpm)[grep("HLA.", rownames(log2tpm))] <- gsub(".", "-", temp, fixed = TRUE)


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
# Reviewer 1 comment 1 (Option 2: much better correlation) Percent (%) of fibroblasts have >0 marker; ignore clusters
# ------------------------------------------------------------------------------------------------------------------------
# For each cell type single-cell RNA-seq data:
# Generate all the subplots all at once
# "Fibroblast", "Monocyte", "B cell", "T cell"
type <- "Monocyte"

# Two marker genes per cluster
# markers <- c("HLA-DRA", "IFI30", "PTGFR", "CD34", "DKK3", "COL8A2", "CLIC5", "HBEGF")
markers <- c("NR4A2", "ATF3", "NUPR1", "HTRA1", "CD14", "MARCO", "IFI6", "IFITM3")
# markers <- c("IGHM", "CXCR4", "HLA-DRA", "IGHG3", "ITGAX", "ZEB2", "MZB1", "XBP1")
# markers <- c("CCR7", "SELL", "FOXP3", "TIGIT", "CXCL13", "PDCD1", "GZMK", "NKG7", "GZMB", "PRF1", "HLA-DQA1", "HLA-DRB1")

log2cpm_type <- log2cpm[, which(meta$cell_type == type)]
meta_type <- meta[which(meta$cell_type == type),]
all(colnames(log2cpm_type) == meta_type$cell_name)

log2tpm_type <- log2tpm[, which(bulk_meta$Cell.type == type)]
bulk_meta_type <- bulk_meta[which(bulk_meta$Cell.type == type),]
all(colnames(log2tpm_type) == bulk_meta_type$cell_name)

myplots <- list()
for (i in 1:length(markers)) {
  gene <- markers[i]
  meta_type$gene <- as.numeric(log2cpm_type[which(rownames(log2cpm_type) == gene),])
  dat_percent <- meta_type %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(percent = sum(gene > 0) / length(gene),
                     ave_mean = mean(gene),
                     ave_median = median(gene))
  dat_percent <- as.data.frame(dat_percent)
  dat_percent <- dat_percent[order(dat_percent$sample),]

  bulk_meta_type$gene_bulk <- as.numeric(log2tpm_type[which(rownames(log2tpm_type) == gene),]) 
  
  # Intersect
  inter <- intersect(meta_type$sample, bulk_meta_type$Donor.ID)
  meta_1 <- dat_percent[which(dat_percent$sample %in% inter),]
  meta_2 <- bulk_meta_type[which(bulk_meta_type$Donor.ID %in% inter),]
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
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  # Linear regression
  fit <- lm(gene_bulk ~ percent, data = meta_1)
  fit$model$Mahalanobis <- meta_1$Mahalanobis
  lm_eqn <- function(fit){
    eq <- substitute(# italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
      ~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
      list(
        # a = format(coef(m)[1], digits = 2), 
        # b = format(coef(m)[2], digits = 2), 
        r2 = format(summary(fit)$r.squared, digits = 2),
        pvalue = format(summary(fit)$coefficients[2,'Pr(>|t|)'], digits=1)
      )
    )
    as.character(as.expression(eq));                 
  }
  ind <- paste("p", i, sep = "")
  ind <- ggplot(
    fit$model,
    aes(
      x=gene_bulk,
      y=percent,
      fill = Mahalanobis
      )
    ) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
                size = 0.7, linetype="dashed",
                col= "black", fill="lightgrey") +
    geom_point(
      shape=21, size = 5, stroke = 0.1
    ) +
    geom_text(x = -Inf, y = Inf, hjust = 0, vjust = 1.8, # display the text on the top left; for fibroblast: vjust = 1.8
              label=lm_eqn(fit), 
              color='black', fontface="plain", size=7, parse=T) +
    scale_fill_manual(values = meta_colors$Case.Control) +
    labs(
      x = NULL,
      y = NULL,
      title = gene
      # x = "Bulk RNA-seq expression",
      # y = "Percent of non-zero expressing cells",
      # title = paste0(gene, " in ", type)
    ) +
    scale_x_continuous(labels=scaleFUN) +
    # scale_x_continuous(breaks= pretty_breaks()) + # display integer values in axis
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_classic(base_size = 30) +
    theme(
      legend.position = "none",
      # axis.text = element_blank(),
      # axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(color="black", size=27)
    ) 
  myplots[[i]] <- ind
}

all <- do.call("grid.arrange", c(myplots, ncol = 4))
# ggsave(file = paste("bulk_percSC_", type, ".pdf", sep = ""), all, width = 27, height = 9, dpi = 300)
ggsave(file = paste("bulk_percSC_", type, ".pdf", sep = ""), all, width = 18, height = 9, dpi = 300)
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
               # Tph are cells: PDCD1+, CXCL13+, ICOS+, CXCR5-
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


# ---------------------------------------------
# Update Fig 5a
markers_4cluster <- read.xls("../../HMS/amp/results/2017_11_05_Run_knn_community_detection_on_cca_mono/markers_4clusters_mono.xlsx")
markers_4cluster$gene <- as.character(markers_4cluster$gene)
sc_genes <- c(markers_4cluster$gene[c(1:20, 25, 38)],
               markers_4cluster$gene[101:120],
               arkers_4cluster$gene[200:220],
               markers_4cluster$gene[c(300:315,384)])

log2cpm_mono <- log2cpm[, which(meta$cell_type == "Monocyte")]
meta_mono <- meta[which(meta$cell_type == type),]
all(colnames(log2cpm_mono) == meta_mono$cell_name)

sc_exp <- log2cpm_mono[sc_genes,]
dim(sc_exp)
mat_breaks <- seq(min(sc_exp), max(sc_exp), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(sc_exp, n = 11)

annotation_col <- meta_mono[, c("plate", "Mahalanobis_20", "cluster")]
colnames(annotation_col)[3] <- "fine_cluster"
colnames(annotation_col)[2] <- "disease"
rownames(annotation_col) <- meta_mono$cell_name
rownames(meta_mono) <- meta_mono$cell_name
sc_exp <- sc_exp[,order(annotation_col$fine_cluster)]
scale_rows <- function(x) t(scale(t(x)))
sc_exp <- scale_rows(sc_exp) # Z-score
sc_exp[sc_exp > 2] <- 2
sc_exp[sc_exp < -2] <- -2

pdf("heatmap_markers_mono.pdf", width=8, height=4, onefile = FALSE, bg = "white")
pheatmap(
  mat = sc_exp,
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(6),
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = meta_colors,
  fontsize = 10,
  fontsize_row = 7,
  scale = "none"
)
dev.off()

# ---------------------------------------------
# Update Fig 7a
sc_genes <- c("IGHD", "CXCR4", "IL6", "IGHM", "CD83", "BACH2", "KDM6B", "NR4A2", "CD69",
              "CD74", "HLA-DPB1", "MS4A1", "HLA-DRA", "HLA-DRB1", "SELL",
              "ITGAX", "ZEB2", "ACTB", "CD52", "TBX21", "IFI44L", "OAS3", "GBP1", "ISG15", # "RSAD2", "IFIT1",
              "SSR4", "MZB1", "FKBP11", "XBP1", "DERL3", "SLAMF7", "CD27")

log2cpm_bcell <- log2cpm[, which(meta$cell_type == "B cell")]
meta_bcell <- meta[which(meta$cell_type == "B cell"),]
all(colnames(log2cpm_bcell) == meta_bcell$cell_name)

sc_exp <- log2cpm_bcell[sc_genes,]
dim(sc_exp)
mat_breaks <- seq(min(sc_exp), max(sc_exp), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(sc_exp, n = 11)

annotation_col <- meta_bcell[, c("plate", "Mahalanobis_20", "cluster")]
colnames(annotation_col)[3] <- "fine_cluster"
colnames(annotation_col)[2] <- "disease"
rownames(annotation_col) <- meta_bcell$cell_name
rownames(meta_bcell) <- meta_bcell$cell_name
sc_exp <- sc_exp[,order(annotation_col$fine_cluster)]
scale_rows <- function(x) t(scale(t(x)))
sc_exp <- scale_rows(sc_exp) # Z-score
sc_exp[sc_exp > 2] <- 2
sc_exp[sc_exp < -2] <- -2

pdf("heatmap_markers_bcell.pdf", width=6, height=3, onefile = FALSE, bg = "white")
pheatmap(
  mat = sc_exp,
  border_color = NA,
  color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(8),
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

# -------------------------------------
# Plot heatmap of donors and clusters
samp_clu <- table(meta$sample, meta$cluster)
samp_clu <- as.data.frame(samp_clu)
colnames(samp_clu) <- c("sample", "cluster", "cells")
# Add leukocyte-rich, ..., labels
inflam_label <- read.xls("data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "Leukocyte-rich RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "Leukocyte-poor RA"
inflam_label <- inflam_label[which(inflam_label$Patient %in% samp_clu$sample),]
inflam_label <- inflam_label[, c("Patient", "Mahalanobis_20")]
colnames(inflam_label)[1] <- "sample"
samp_clu_merge <- merge(samp_clu, inflam_label, by = "sample")

# 'percent_donor': is the percent of the donor's cells in this cluster,
# where all four cell types are pooled together.
# samp_clu_merge %<>% group_by(sample) %>% mutate(percent_donor = 100 * cells / sum(cells))
# mat <- dcast(data = samp_clu, formula = sample ~ cluster, value.var = "percent_donor")
mat <- dcast(data = samp_clu_merge, formula = sample ~ cluster, value.var = "cells")
rownames(mat) <- mat[[1]]
mat[[1]] <- NULL
mat <- as.matrix(mat)
# All in one heatmap: number of cells in one donor and one cluster
pheatmap(
  mat = mat,
  color = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10
)

# ---
# Saperate by cell type
this <- as.data.frame(samp_clu_merge)
title = "Cell count: T cell"
# this_celltype <- this[which(this$cluster %in% c("SC-F1", "SC-F2", "SC-F3", "SC-F4")),]
# this_celltype <- this[which(this$cluster %in% c("SC-B1", "SC-B2", "SC-B3", "SC-B4")),]
this_celltype <- this[which(this$cluster %in% c("SC-T1", "SC-T2", "SC-T3", "SC-T4", "SC-T5", "SC-T6")),]
# this_celltype <- this[which(this$cluster %in% c("SC-M1", "SC-M2", "SC-M3", "SC-M4")),]
this_celltype %<>% group_by(sample) %>% mutate(percent_donor = 100 * cells / sum(cells))
# ggplot() +
#   geom_tile(
#     data = this_celltype,
#     mapping = aes(x = cluster, y = sample, fill = percent_donor)
#   ) +
#   # facet_grid( ~ Mahalanobis_20) +
#   # scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(9)) +
#   scale_fill_gradientn(colours = viridis(10)) +
#   theme_minimal(base_size = 15) +
#   theme(
#     # strip.text = element_blank(),
#     strip.background = element_blank(),
#     # axis.ticks = element_blank(), 
#     # panel.grid = element_blank(),
#     axis.text.x = element_text(angle = 35, hjust = 1)
#   ) +
#   labs(x = NULL, y = NULL, title = "Fibroblast")


# One cell type heatmap: number of cells in one donor and one cluster
# this_celltype <- dcast(data = this_celltype, formula = sample ~ cluster, value.var = "percent_donor")
this_celltype <- dcast(data = this_celltype, formula = sample ~ cluster, value.var = "cells")
rownames(this_celltype) <- this_celltype[[1]]
this_celltype[[1]] <- NULL
this_celltype <- as.matrix(this_celltype)

# Only plot the 21 post-QC samples
this_celltype <- this_celltype[which(rownames(this_celltype) %in% names(which(table(meta$sample) > 0))), ]
temp <- inflam_label[which(inflam_label$Patient %in% rownames(this_celltype)),]
colnames(temp) <- c("sample", "disease")
annotation_row = data.frame(
  disease = temp$disease
)
rownames(annotation_row) <- temp$sample
annotation_col = data.frame(
  cluster = colnames(this_celltype)
)
rownames(annotation_col) <- annotation_col$cluster
meta_colors$cluster <- meta_colors$fine_cluster

this_celltype <- as.data.frame(this_celltype)
this_celltype$sample <- rownames(this_celltype)
# Use the order of OA, leukocyte-poor RA, leukocyte-rich RA
this_celltype <- this_celltype[order(match(this_celltype$sample, rownames(annotation_row))), ]
all(rownames(this_celltype) == rownames(annotation_row))
this_celltype <- this_celltype[, -which(colnames(this_celltype) == "sample")]

# Overwrite default draw_colnames in the pheatmap package.
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_50 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 50, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_50",
  ns = asNamespace("pheatmap")
)
pdf(paste(title, ".pdf", sep=""), width=4.8, height=4, onefile = FALSE, bg = "white")
pheatmap(
  mat = this_celltype,
  border_color = "white",
  # color = colorRampPalette(viridis(9))(9),
  color = cet_pal(100, name = "kgy"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = meta_colors,
  show_rownames = FALSE,
  fontsize = 7,
  main = title
)
dev.off()



# ------------------------------------
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

# -------
# fibro
synData.Fibro_sampleID <- synData.Fibro.downsample[, c("sampleID", "CD90", "CD34")]
synData.Fibro_mean <- synData.Fibro_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean_CD90 = mean(CD90),
                   ave_median_CD90 = median(CD90),
                   ave_mean_CD34 = mean(CD34),
                   ave_median_CD34 = median(CD34)
                   )
synData.Fibro_mean <- as.data.frame(synData.Fibro_mean)

# B cell
synData.Bcell_sampleID <- synData.Bcell.downsample[, c("sampleID", "CD11c", "IgD")]
synData.Bcell_mean <- synData.Bcell_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean_CD11c = mean(CD11c),
                   ave_median_CD11c = median(CD11c),
                   ave_mean_IgD = mean(IgD),
                   ave_medianIgD = median(IgD)
                   )
synData.Bcell_mean <- as.data.frame(synData.Bcell_mean)

# Mono
synData.Mono_sampleID <- synData.Mono.downsample[, c("sampleID", "CD38", "CCR2")]
synData.Mono_mean <- synData.Mono_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean_CD38 = mean(CD38),
                   ave_median_CD38 = median(CD38),
                   ave_mean_CCR2 = mean(CCR2),
                   ave_median_CCR2 = median(CCR2)
                   )
synData.Mono_mean <- as.data.frame(synData.Mono_mean)

# T cell
synData.Tcell_sampleID <- synData.Tcell.downsample[, c("sampleID", "CD8a", "PD.1")]
synData.Tcell_mean <- synData.Tcell_sampleID %>%
  dplyr::group_by(sampleID) %>%
  dplyr::summarise(ave_mean_CD8a = mean(CD8a),
                   ave_median_CD8a = median(CD8a),
                   ave_mean_PD1 = mean(PD.1),
                   ave_median_PD1 = median(PD.1)
                   )
synData.Tcell_mean <- as.data.frame(synData.Tcell_mean)

# -------
# bulk data
# fibro
bulk_meta$CD90_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "THY1"), ])
bulk_meta$CD34_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CD34"), ])
cell_type <- "Fibroblast"
log2tpm_fibro <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_fibro <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_fibro) == bulk_meta_fibro$Sample.ID)
dim(bulk_meta_fibro)
dim(log2tpm_fibro)

# B cell
# bulk_meta$MS4A1 <- as.numeric(log2tpm[which(rownames(log2tpm) == "MS4A1"), ])
bulk_meta$ITGAX <- as.numeric(log2tpm[which(rownames(log2tpm) == "ITGAX"), ])
# bulk_meta$HLA.DRA_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "HLA.DRA"), ])
# bulk_meta$CD19_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CD19"), ])
# bulk_meta$IGHM_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "IGHM"), ])
bulk_meta$IGHD_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "IGHD"), ])
cell_type <- "B cell"
log2tpm_bcell <- log2tpm[, which(bulk_meta$Cell.type == cell_type)]
bulk_meta_bcell <- bulk_meta[which(bulk_meta$Cell.type == cell_type),]
all(colnames(log2tpm_bcell) == bulk_meta_bcell$Sample.ID)
dim(bulk_meta_bcell)
dim(log2tpm_bcell)

# Mono
bulk_meta$CD38_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CD38"), ])
bulk_meta$CCR2_gene <- as.numeric(log2tpm[which(rownames(log2tpm) == "CCR2"), ])
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

# -------
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

# -------
# Overlapped samples between bulk fibroblast samples with cytof fibroblast samples
# Fibro
over <- intersect(synData.Fibro_mean$sampleID, bulk_meta_fibro$Donor.ID)
synData.Fibro_mean <- synData.Fibro_mean[which(synData.Fibro_mean$sampleID %in% over), ]
bulk_meta_fibro <- bulk_meta_fibro[which(bulk_meta_fibro$Donor.ID %in% over), ]
bulk_meta_fibro <- bulk_meta_fibro[order(match(bulk_meta_fibro$Donor.ID, synData.Fibro_mean$sampleID)), ]
all(synData.Fibro_mean$sampleID == bulk_meta_fibro$Donor.ID)
bulk_meta_fibro$cytof_ave_mean_CD90 <- synData.Fibro_mean$ave_mean_CD90
bulk_meta_fibro$cytof_ave_median_CD90 <- synData.Fibro_mean$ave_median_CD90
bulk_meta_fibro$cytof_ave_mean_CD34 <- synData.Fibro_mean$ave_mean_CD34
bulk_meta_fibro$cytof_ave_median_CD34 <- synData.Fibro_mean$ave_median_CD34

# B cell
over <- intersect(synData.Bcell_mean$sampleID, bulk_meta_bcell$Donor.ID)
synData.Bcell_mean <- synData.Bcell_mean[which(synData.Bcell_mean$sampleID %in% over), ]
bulk_meta_bcell <- bulk_meta_bcell[which(bulk_meta_bcell$Donor.ID %in% over), ]
bulk_meta_bcell <- bulk_meta_bcell[order(match(bulk_meta_bcell$Donor.ID, synData.Bcell_mean$sampleID)), ]
all(synData.Bcell_mean$sampleID == bulk_meta_bcell$Donor.ID)
bulk_meta_bcell$cytof_ave_mean_CD11c <- synData.Bcell_mean$ave_mean_CD11c
bulk_meta_bcell$cytof_ave_median_CD11c <- synData.Bcell_mean$ave_median_CD11c
bulk_meta_bcell$cytof_ave_mean_IgD <- synData.Bcell_mean$ave_mean_IgD
bulk_meta_bcell$cytof_ave_median_IgD <- synData.Bcell_mean$ave_medianIgD

# Mono
over <- intersect(synData.Mono_mean$sampleID, bulk_meta_mono$Donor.ID)
synData.Mono_mean <- synData.Mono_mean[which(synData.Mono_mean$sampleID %in% over), ]
bulk_meta_mono <- bulk_meta_mono[which(bulk_meta_mono$Donor.ID %in% over), ]
bulk_meta_mono <- bulk_meta_mono[order(match(bulk_meta_mono$Donor.ID, synData.Mono_mean$sampleID)), ]
all(synData.Mono_mean$sampleID == bulk_meta_mono$Donor.ID)
bulk_meta_mono$cytof_ave_mean_CD38 <- synData.Mono_mean$ave_mean_CD38
bulk_meta_mono$cytof_ave_median_CD38 <- synData.Mono_mean$ave_median_CD38
bulk_meta_mono$cytof_ave_mean_CCR2 <- synData.Mono_mean$ave_mean_CCR2
bulk_meta_mono$cytof_ave_median_CCR2 <- synData.Mono_mean$ave_median_CCR2

# T cell
over <- intersect(synData.Tcell_mean$sampleID, bulk_meta_tcell$Donor.ID)
synData.Tcell_mean <- synData.Tcell_mean[which(synData.Tcell_mean$sampleID %in% over), ]
bulk_meta_tcell <- bulk_meta_tcell[which(bulk_meta_tcell$Donor.ID %in% over), ]
bulk_meta_tcell <- bulk_meta_tcell[order(match(bulk_meta_tcell$Donor.ID, synData.Tcell_mean$sampleID)), ]
all(synData.Tcell_mean$sampleID == bulk_meta_tcell$Donor.ID)
bulk_meta_tcell$cytof_ave_mean_CD8a <- synData.Tcell_mean$ave_mean_CD8a
bulk_meta_tcell$cytof_ave_median_CD8a <- synData.Tcell_mean$ave_median_CD8a
bulk_meta_tcell$cytof_ave_mean_PD1 <- synData.Tcell_mean$ave_mean_PD1
bulk_meta_tcell$cytof_ave_median_PD1 <- synData.Tcell_mean$ave_median_PD1



# Plot fitted line -------
# Fibro
fit <- lm(CD34_gene ~ cytof_ave_mean_CD34, data = bulk_meta_fibro)
lm_eqn <- function(fit){
  eq <- substitute(# italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
    ~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
    list(
      # a = format(coef(m)[1], digits = 2), 
      # b = format(coef(m)[2], digits = 2), 
      r2 = format(summary(fit)$r.squared, digits = 3),
      pvalue = format(summary(fit)$coefficients[2,'Pr(>|t|)'], digits=1)
    )
  )
  as.character(as.expression(eq));                 
}
ggplot(
  fit$model,
  aes(
    x=CD34_gene,
    y=cytof_ave_mean_CD34
   )
  ) +
  geom_text(x=2.2, y=2.7, label=lm_eqn(fit), color='black', fontface="plain", size=4, parse=T) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 0.7, linetype="dashed",
              col= "black", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#08519C"
  ) +
  labs(
    x = "",
    y = "",
    title = "CD34"
    # title = "Fibroblast samples\nCD90"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CD34_fit", ".pdf", sep = ""),
       width = 3, height = 3, dpi = 300)
dev.off()


# Mono
fit <- lm(CCR2_gene ~ cytof_ave_mean_CCR2, data = bulk_meta_mono)
lm_eqn <- function(fit){
  eq <- substitute(# italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
    ~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
    list(
      # a = format(coef(m)[1], digits = 2), 
      # b = format(coef(m)[2], digits = 2), 
      r2 = format(summary(fit)$r.squared, digits = 3),
      pvalue = format(summary(fit)$coefficients[2,'Pr(>|t|)'], digits=1)
    )
  )
  as.character(as.expression(eq));                 
}
ggplot(
  fit$model,
  aes(
    x=CCR2_gene,
    y=cytof_ave_mean_CCR2
    )
  ) +
  geom_text(x=1.5, y=3, label=lm_eqn(fit), color='black', fontface="plain", size=4, parse=T) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 0.7, linetype="dashed",
              col= "black", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#DE77AE"
  ) +
  labs(
    x = "",
    y = "",
    title = "CCR2"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_CCR2_fit", ".pdf", sep = ""),
       width = 3, height = 3.2, dpi = 300)
dev.off()


# B cell
fit <- lm(IGHD_gene ~ cytof_ave_mean_IgD, data = bulk_meta_bcell)
lm_eqn <- function(fit){
  eq <- substitute(# italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
    ~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
    list(
      # a = format(coef(m)[1], digits = 2), 
      # b = format(coef(m)[2], digits = 2), 
      r2 = format(summary(fit)$r.squared, digits = 3),
      pvalue = format(summary(fit)$coefficients[2,'Pr(>|t|)'], digits=1)
    )
  )
  as.character(as.expression(eq));                 
}
ggplot( 
  fit$model,
  aes(
    x=IGHD_gene,
    y=cytof_ave_mean_IgD
   )
  ) +
  geom_text(x=5.5, y=2.1, label=lm_eqn(fit), color='black', fontface="plain", size=4, parse=T) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 0.7, linetype="dashed",
              col= "black", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#E31A1C"
  ) +
  labs(
    x = "",
    y = "",
    title = "IgD"
  ) +
  scale_y_continuous(breaks=c(0,1,2)) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_IgD_fit", ".pdf", sep = ""),
       width = 3, height = 3, dpi = 300)
dev.off()


# T cell
fit <- lm(PDCD1_gene ~ cytof_ave_mean_PD1, data = bulk_meta_tcell)
lm_eqn <- function(fit){
  eq <- substitute(# italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
    ~~italic(r)^2~"="~r2*","~~italic(p)~"="~pvalue,
    list(
      # a = format(coef(m)[1], digits = 2), 
      # b = format(coef(m)[2], digits = 2), 
      r2 = format(summary(fit)$r.squared, digits = 3),
      pvalue = format(summary(fit)$coefficients[2,'Pr(>|t|)'], digits=1)
    )
  )
  as.character(as.expression(eq));                 
}
ggplot(
  fit$model,
  aes(
    x=PDCD1_gene,
    y=cytof_ave_mean_PD1
    )
  ) +
  geom_text(x=4.95, y=2.2, label=lm_eqn(fit), color='black', fontface="plain", size=4, parse=T) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, 
              size = 0.7, linetype="dashed",
              col= "black", fill="lightgrey") +
  geom_point(
    shape=21, size = 2.5, stroke = 0.1, fill = "#A65628"
  ) +
  labs(
    x = "",
    y = "",
    title = "PD1"
  ) +
  scale_y_continuous(breaks=c(0,1,2)) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none"
    # axis.text = element_blank()
    # axis.ticks = element_blank()
    # panel.grid = element_blank()
  ) 
ggsave(file = paste("cytof_bulk_per_sample_PD1_fit", ".pdf", sep = ""),
       width = 3, height = 3.2, dpi = 300)
dev.off()


# # Fibro
# ggplot() +
#   geom_point(
#     data = bulk_meta_fibro,
#     # mapping = aes_string(x = "CD90", y = "cytof_ave_mean", fill = "Mahalanobis_20"),
#     mapping = aes_string(x = "CD34_gene", y = "cytof_ave_mean_CD90"),
#     size = 3, stroke = 0.2, shape = 20, color = "#08519C"
#   ) +
#   # scale_fill_manual(values = meta_colors$disease, name = "") +
#   labs(
#     x = "Bulk RNA-seq expression",
#     y = "Averaged mass cytometry \n proteomic expression",
#     title = "Fibroblast samples\nCD34"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none"
#     # axis.text = element_blank()
#     # axis.ticks = element_blank()
#     # panel.grid = element_blank()
#   ) 
# ggsave(file = paste("cytof_bulk_per_sample_CD34", ".pdf", sep = ""),
#        width = 3, height = 3, dpi = 300)
# dev.off()
# 
# # B cell
# ggplot() +
#   geom_point(
#     data = bulk_meta_bcell,
#     mapping = aes_string(x = "IGHD_gene", y = "cytof_ave_mean"),
#     size = 3, stroke = 0.2, shape = 20, color = "#E31A1C"
#   ) +
#   # scale_fill_manual(values = meta_colors$disease, name = "") +
#   labs(
#     x = "Bulk RNA-seq expression",
#     y = "Averaged mass cytometry \n proteomic expression",
#     title = "B cell samples\nIgD"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none"
#     # axis.text = element_blank()
#     # axis.ticks = element_blank()
#     # panel.grid = element_blank()
#   ) 
# ggsave(file = paste("cytof_bulk_per_sample_IGH", ".pdf", sep = ""),
#        width = 3, height = 3, dpi = 300)
# dev.off()
# 
# # Mono
# ggplot() +
#   geom_point(
#     data = bulk_meta_mono,
#     mapping = aes_string(x = "CCR2_gene", y = "cytof_ave_mean"),
#     size = 3, stroke = 0.2, shape = 20, color = "#DE77AE"
#   ) +
#   # scale_fill_manual(values = meta_colors$disease, name = "") +
#   labs(
#     x = "Bulk RNA-seq expression",
#     y = "Averaged mass cytometry \n proteomic expression",
#     title = "Monocyte samples\nCCR2"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none"
#     # axis.text = element_blank()
#     # axis.ticks = element_blank()
#     # panel.grid = element_blank()
#   ) 
# ggsave(file = paste("cytof_bulk_per_sample_CCR2", ".pdf", sep = ""),
#        width = 3, height = 3, dpi = 300)
# dev.off()
# 
# 
# # T cell
# ggplot() +
#   geom_point(
#     data = bulk_meta_tcell,
#     mapping = aes_string(x = "PDCD1_gene", y = "cytof_ave_mean"),
#     size = 3, stroke = 0.2, shape = 20, color = "#A65628"
#   ) +
#   # scale_fill_manual(values = meta_colors$disease, name = "") +
#   labs(
#     x = "Bulk RNA-seq expression",
#     y = "Averaged mass cytometry \n proteomic expression",
#     title = "T cell samples\nPD1"
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none"
#     # axis.text = element_blank()
#     # axis.ticks = element_blank()
#     # panel.grid = element_blank()
#   ) 
# ggsave(file = paste("cytof_bulk_per_sample_PDCD1", ".pdf", sep = ""),
#        width = 3.2, height = 3, dpi = 300)
# dev.off()


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

# --------------------------------------------
# Reviewer 2 Comment 5.1 
# Calculate entropy or KL-divergence in R 
meta$sample <- as.character(meta$sample)
table(meta$sample)
table(meta$cluster)

library(factoextra)
library(cluster)

file_mean_sd <- "celseq_synovium_log2cpm_mean_sd.rds"
if (!file.exists(file_mean_sd)) {
  # dat_mean_sd <- data.frame(
  #   mean  = Matrix::rowMeans(dat),
  #   sd    = apply(dat, 1, sd),
  #   count = apply(dat, 1, function(x) sum(x > 0))
  # )
  dat_mean_sd <- data.frame(
    mean  = Matrix::rowMeans(log2cpm_fibro),
    sd    = apply(log2cpm_fibro, 1, sd),
    count = apply(log2cpm_fibro, 1, function(x) sum(x > 0))
  )
  xs <- dat_mean_sd$mean
  xs[xs == 0] <- 1
  dat_mean_sd$cv <- dat_mean_sd$sd / xs
  rm(xs)
  dat_mean_sd$density_cv <- with(dat_mean_sd, get_density(mean, cv))
  dat_mean_sd$density_sd <- with(dat_mean_sd, get_density(mean, sd))
  saveRDS(dat_mean_sd, file_mean_sd)
} else {
  dat_mean_sd <- readRDS(file_mean_sd)
}

top_percent <- 0.9
mat_idx <- with(
  dat_mean_sd,
  mean > quantile(mean, top_percent) &
    sd > quantile(sd, top_percent)
)

test <- log2cpm_fibro[mat_idx,]
dim(test)

dat_dist <- dist(as.matrix(t(test)))
dat_dist <- as.matrix(dat_dist)
dim(dat_dist)
# per = 50
# tsne1 <- Rtsne(
#   X = log2cpm_fibro,
#   is_distance = TRUE,
#   dims = 2,
#   perplexity = per,
#   theta = 0.5,
#   pca = TRUE
# )

meta_fibro$clus_num <- substring(meta_fibro$cluster, 5)
sil = silhouette(x = as.numeric(meta_fibro$clus_num), dmatrix = 1-dat_dist^2)
plot(sil,
     col = meta_colors$fine_cluster, 
     xlab = "Silhouette width",
     ylab = "Cells",
     main = "Silhouette plot of fibroblast clusters and \ncell-to-cell similarity matrix",
     cex.lab = 1,
     cex = 1,
     cex.axis = 1
)


# # https://github.com/MarioniLab/MNN2017/blob/master/SomeFuncs/BatchMixingEntropy.R
# BatchEntropy <- function(dataset, batch0, L=100, M=100, k=500) {
#   # entropy of batch mixing
#   # L is the number bootstrapping times
#   # M is the number of randomly picked cells    
#   # k is the number of nearest neighbours of cell (from all batches) to check   
#   
#   require(RANN)  
#   nbatches<-length(unique(batch0))
#   
#   entropy<-matrix(0,L,1)
#   set.seed(0) 
#   for (boot in 1:L) {
#     bootsamples<-sample(1:nrow(dataset),M)
#     W21<-nn2(dataset,query=dataset[bootsamples,],k)
#     
#     for (i in 1:length(bootsamples)){
#       
#       for (j in 1:nbatches) {
#         xi<-max(1,sum(batch0[W21$nn.idx[i,]]==j))
#         entropy[boot]<-entropy[boot]+xi*log(xi)
#       }
#     }
#   }
#   
#   return( (-1)*entropy/length(bootsamples) )
# }
# 
# en <- BatchEntropy(meta[, c("T1_all", "T2_all")], meta$sample, L=100, M=100, k=500)
# boxplot(en)

find_entropy <- function(donor_list) {
  entropies <- c()
  for (i in 1:length(donor_list)) {
    donor <- donor_list[[i]]
    distributions <- c()
    entropy_i <- c()
    # d <- length(table(donor$cell_type))
    # for (j in 1:length(table(donor$cell_type))) {
    #   distributions[j] <- (table(donor$cell_type)[j])/(sum(table(donor$cell_type)))
    #   entropy_i[j] <- distributions[j]*log(distributions[j], d)
    # } 
    d <- length(table(donor$cluster))
    for (j in 1:length(table(donor$cluster))) {
      distributions[j] <- (table(donor$cluster)[j])/(sum(table(donor$cluster)))
      entropy_i[j] <- distributions[j]*log(distributions[j], d)
    } 
    entropy <- -sum(entropy_i)
    entropies[i] <- entropy
  }
  entropies
}
meta$sample <- as.character(meta$sample)
meta$cell_type <- as.character(meta$cell_type)
meta$cluster <- as.character(meta$cluster)
  
# Plot entropy of 21 donors based CCA results
donors <- list()
for (i in 1:length(table(meta$sample))){
  donor_name <- names(table(meta$sample))[i]
  donor_name <- meta[meta$sample == donor_name, ]
  donors[[i]] <- donor_name
}
entropies <- find_entropy(donors)

plot_entropy = data.frame(
  donor =names(table(meta$sample)),
  entropy = entropies
)
ggplot(
  data=plot_entropy,
  aes(x=donor, y= entropy)
  ) +
  geom_bar(stat="identity",
           position = "stack",
           width = 0.7
  ) +
  # scale_fill_manual("", values = meta_colors$sample) +
  labs(
    x = "donors",
    y = "Entropy"
  ) +
  theme_classic(base_size = 20) +
  theme(    
    # legend.position="none",
    # axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x=element_text(angle=30, hjust=1),
    axis.text.y = element_text(size = 20)
  )
ggsave(file = paste("sc_cluster_donor_entropy", ".pdf", sep = ""), width = 9, height = 5, dpi = 300)
dev.off()

# # # Plot entropy of 21 donors based PCA results
# donors_seurat <- list()
# for (i in 1:length(table(meta_seurat$donor))){
#   donor_name <- names(table(meta_seurat$donor))[i]
#   donor_name <- meta_seurat[meta_seurat$donor == donor_name, ]
#   donors_seurat[[i]] <- donor_name
# }
# seurat_entropies <- find_entropy(donors_seurat)
# barplot(seurat_entropies)


# Chamith's function on KL divergence
calcKL <- function(dataset, clusterCol, clusterName, markers, binSize = 100) {
  message(paste("Ranking markers for cluster", clusterName))
  divergence <- rep(NA, length(markers))
  names(divergence) <- markers
  for (i in seq_along(markers)) {
    stain <- markers[i]
    q <- dataset[, stain]
    p <- dataset[dataset[[clusterCol]] == clusterName, stain]
    # Create 100 evenly spaced bins from 0 to maximum value observed for stain
    bins <- seq(0, max(q, p), length.out = binSize)
    # Create vectors of counts at each bin for p and q
    q.counts <- hist(q, breaks = bins, include.lowest = T, plot = F)$counts
    p.counts <- hist(p, breaks = bins, include.lowest = T, plot = F)$counts
    # Normalize so counts sum to 1 (converting to probability vector)
    q.probs <- q.counts/sum(q.counts)
    p.probs <- p.counts/sum(p.counts)
    # Calculate divergence for each bin, setting NaN to 0 and sum
    divr <- p.probs * log(p.probs/q.probs)
    divr[is.na(divr)] <- 0
    divr <- sum(divr)
    divergence[stain] <- divr
  }
  return(divergence)
}

# ---
# Update silhouette plot in the supplemental figure
sil_fibro <- readRDS("/Users/fanzhang/Documents/HMS/amp/results/2017_10_18_Run_knn_community_detection_on_cca_results/fibro_silhouette_width_euclidean.rds")
sil_mono <- readRDS("/Users/fanzhang/Documents/HMS/amp/results/2017_11_05_Run_knn_community_detection_on_cca_mono/mono_silhouette_width_euclidean.rds")
sil_tcell <- readRDS("/Users/fanzhang/Documents/HMS/amp/results/2017_10_25_Run_knn_community_detection_on_cca_Tcell/tcell_silhouette_width_euclidean.rds")
sil_bcell <- readRDS("/Users/fanzhang/Documents/HMS/amp/results/2017_11_05_Run_knn_community_detection_on_cca_bcell/bcell_silhouette_width_euclidean.rds")
all <- rbind.data.frame(sil_fibro, sil_mono, sil_tcell, sil_bcell)

all$cluster_name <- factor(all$cluster_name, 
                           levels=rev(c("SC-F1", "SC-F2", "SC-F3", "SC-F4", "SC-M1", "SC-M2", "SC-M3", "SC-M4",
                                    "SC-T1", "SC-T2", "SC-T3", "SC-T4", "SC-T5", "SC-T6", "SC-B1", "SC-B2", 
                                    "SC-B3", "SC-B4"))
                           )
# all$cluster_numCells <- paste0(all$cluster_name, "(n=", all$cell, ")", sep="")

ggplot(all, 
       aes(y=sil_width, x= cluster_name, color=cluster_name)) + 
  geom_boxplot() +
  scale_color_manual(values = meta_colors$fine_cluster, name = "") +
  geom_hline(yintercept=0, linetype = "dashed", color = "grey") + 
  # coord_cartesian(ylim = c(-0.4, 0.8)) +
  coord_flip(ylim = c(-0.2, 0.9)) +
  # scale_x_discrete(position = "top") +
  labs(
    x = "",
    y = "Silhouette width"
  ) +
  theme_bw(base_size = 30) +
  theme(axis.ticks = element_blank(),
        panel.grid = element_blank()
  ) 
ggsave(file = paste("silhouette_width_eucliDist_allClusters", ".pdf", sep = ""), 
       width = 10, height = 10)
dev.off()


# ---
# Update Figure 2 cytof panel

synData.3ksample <- readRDS("R/synData.3ksample.postQC.SNE.rds")
pal <- inferno(40)[c(11, seq(12, 40, 2))]

mat <- synData.3ksample[sample(nrow(synData.3ksample)),]
mat <- mat[, which(colnames(mat) %in% c("SNE1", "SNE2", "group"))]

p <- ggplot(mat, aes(x = SNE1, y = SNE2))
p <- p + geom_point(color = "grey40", size = 0.0002)
p <- p + stat_density2d(geom = "polygon", aes(fill = ..level.., alpha = ..level..), n = 30)
p <- p + coord_fixed(xlim = c(-35, 35), ylim = c(-35, 35))  
p <- p + scale_fill_gradientn(colors = pal)
p <- p + scale_alpha_continuous(range = c(0.0, 0.75), guide = F)
p <- p + facet_grid(. ~ group)
p

ggsave(p, filename = "Figure2.png", width = 10, height = 8, dpi = 400)

# ---
# Update Figure 3 tSNE plots

dat <- readRDS("Documents/GitHub/amp_phase1_ra/data/celseq_synovium_log2_5265cells_paper.rds")
meta <- readRDS("Documents/GitHub/amp_phase1_ra/data/celseq_synovium_meta_5265cells_paper.rds")
all(colnames(dat) == meta$cell_name)

type <- "Fibroblast"
marker <- c("THY1", "HLA-DRA", "IL6", "CD34", "CD55", "DKK3")
# marker <- c("CD14", "IL1B", "NUPR1", "SPP1", "IFI6")
# marker <- c("CD4", "CD8A", "HLA-DRB1", "CCR7", "FOXP3", "CXCL13",
#             "GZMB", "GNLY", "GZMK")
# marker <- c("HLA-DRA", "XBP1", "IGHD", "IGHG3", "CD27", "ITGAX")

myplots <- list()
for (i in 1:length(marker)) {
  gene <- marker[i]
  type_meta <- meta[which(meta$cell_type == type), ]
  type_meta$gene <- as.numeric(dat[which(rownames(dat) == gene), which(meta$cell_type == type)])
  ind <- paste("p", i, sep = "")
  ind <- ggplot() +
    geom_point(
      data = type_meta,
      # data = type_meta[order(type_meta$gene),],
      mapping = aes_string(x = "T1", y = "T2", fill = "gene"),
      size = 1.4, stroke = 0.05, shape = 21
    ) +
    scale_fill_gradientn(
      colours = colorRampPalette(RColorBrewer::brewer.pal(8, "Greens"))(8),
      name = bquote("Log"[2]~"(CPM+1)")
    ) +
    guides(
      # fill = guide_colorbar(barwidth = 1, barheight = 10),
      fill = FALSE,
      alpha = "none"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = gene
    ) +
    theme_bw(base_size = 2) +
    theme(
      plot.title = element_text(color="black", size=25, face = "italic"), # face="bold.italic"
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      panel.grid = element_blank()
    ) 
  # print(ind)
  myplots[[i]] <- ind
}

all <- do.call("grid.arrange", c(myplots, ncol = 2, nrow = 3))
ggsave(file = paste("tsne_markers_", type, ".png", sep = ""), all, width = 4, height = 6.5, dpi = 400)
dev.off()

