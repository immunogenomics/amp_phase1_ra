#' ---
#' title: "Visualizing CCA results for integrating bulk RNA-seq with single-cell RNA-seq for each cell types"
#' author: "Fan Zhang"
#' 
#' Fine clustering of CCA-based pipleline is performed on each cell type
#' We are visualizing the linear projections of genes and cells in the canonical variates,
#' which are results by running `scRNAseq_bulkRNAseq_integrative_pipeline.R` one each cell type data
#' 
#' date: "2018-03-22"
#' ---


setwd("/Users/fanzhang/Documents/HMS/amp/results/2017_06_26_CCA_bulk_singlecell")

library(CCA)
library(glmnet)
library("gridExtra")
library(Matrix)
library(fastcluster)
library(MASS)
library(reshape2)
library(stringr)
library(dplyr)
library(readr)
library(parallel)
library(limma)
library(Rtsne)
library(scde)
library(org.Hs.eg.db)
library(GO.db)
library(ggplot2)
library(scales)
library(ggrepel)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(superheat)
library(pheatmap)
library(biomaRt)
require(gdata)

source("../2017_02_28_Phase1_cellseq_RA_single_cell_data/pure_functions.R")
source("../2017_02_28_Phase1_cellseq_RA_single_cell_data/meta_colors.R")


# Load CCA results for each cell type -----------------------------
# Fibroblasts results
# 7,016 genes; 45 samples; 1,910 cells
res <- readRDS("data/Fibro_cca_7016_genes.rds")

# Load meta data accordingly
bulk_m <- readRDS("data/Fibro_cca_bulk_m.rds")
cell_m <- readRDS("data/Fibro_cca_cell_m.rds")

# T cell results
# 7,003 genes; 47 samples; 1,990 cells
res <- readRDS("data/Tcell_cca_7003_genes.rds")

# B cell results
# 7,023 genes; 29 samples; 1,428 cells
res <- readRDS("data/Bcell_cca_7023_genes.rds")

# Monocyte results
# 7,016 genes; 47 samples; 940 cells
res <- readRDS("data/Monocyte_cca_7016_genes.rds")


# ----------------------------------------------------------------
# Visualization (Take fibroblasts as an example) 
# Plot the canonical variates (CCA dimensions)
barplot(res$cor, xlab = "Dimension", ylab = "Canonical correlations", ylim = c(0,1))

type_colors <- list(
  "data_type" = c(
    "single cell" = "#BC80BD",
    "bulk"   = "#66C2A5"
  )
)

# Plot genes projected in the top CCA dimensions
res_xscores <- data.frame(
  X = res$scores$xscores[, c(1:20)]
)
res_xscores$GENE_NAME <- as.factor(rownames(res_xscores))

# Only plot a few marker genes
marker_Fibro =  c("HTRA1", "PRG4", "APOE", "CXCL12", "MMP2", 
                  "THY1", "HLA-DQA1", "COL14A1", "CD55", "MMP3",
                  "SHANK2", "ENPP1", "PODN")
size_label = 5
ggplot() +
  geom_point(
    data = res_xscores,
    aes(x = X.1, y = X.2),
    size = 0.5,
    color = "dimgrey"
  ) +
  geom_label_repel(
    data = subset(res_xscores, GENE_NAME %in% marker_Fibro),
    aes(x = X.1, y = X.2, label = GENE_NAME),
    size = size_label, color = "black",
    fontface = 'bold',
    box.padding = unit(-0.15, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  #scale_y_continuous(limits = c(-12, 12)) +
  #scale_x_continuous(limits = c(-5, 5)) +
  xlab("Dimension 1 (CCA)") + ylab("Dimension 2 (CCA)")+
  theme(axis.text = element_text(size = 25),
        axis.text.y = element_text(size = 25)) +
  theme_bw(base_size = 25) 
ggsave(file = paste("Fibro_cca_marker_genes", ".png", sep = ""), width = 7, height = 7, dpi = 200)
dev.off()



# ----------------------------------------------------------
# Plot both samples and cells in the same dimensions 
m_bulk <- data.frame(
  res$scores$corr.X.xscores[, 1:20],
  disease = bulk_m$Tissue.type,
  donor = bulk_m$Donor.ID,
  data_type = rep("bulk", nrow(m_bulk))
)
m_sc <- data.frame(
  res$scores$corr.Y.xscores[, 1:20],
  disease = cell_m$disease,
  donor = cell_m$sample,
  data_type = rep("single cell", nrow(m_sc))
)
both <- rbind.data.frame(m_bulk, m_sc)
both$disease[which(both$disease == "OA-arthro")] <- "OA"
both$disease[which(both$disease == "RA-arthro")] <- "RA"
both$disease[which(both$disease == "RA-biopsy")] <- "RA"
both$disease <- as.character(both$disease)

ggplot() +
  # geom_point(
  #   # data = both[sample(nrow(both)),],
  #   data = both,
  #   mapping = aes(m_d1, m_d2, fill = cell_type),
  #   # mapping = aes(m_d1, m_d2, fill = disease), 
  #   shape = 21, size = 4.5, stroke = 0.1
  # ) +
  geom_point(
    data = both[which(both$data_type == "bulk"),],
    mapping = aes(X1, X2, fill = disease),
    shape = 21, size = 4.4, stroke = 0.1
    # alpha = 0.7
  ) +
  geom_point(
    data = both[which(both$data_type == "single cell"),],
    mapping = aes(X1, X2, fill = disease),
    shape = 21, size = 2, stroke = 0.1
    # alpha = 0.6
  ) +
  scale_fill_manual(values = meta_colors$Subject.type, name = "") +
  xlab("Dimension 1 (CCA)") + ylab("Dimension 2 (CCA)")+
  theme_bw(base_size = 25) +
  theme(
    axis.text = element_text(size = 25, colour = "black"),
    axis.text.y = element_text(size = 25)) 
ggsave(file = paste("allcelltypes_cca_celltype", ".png", sep = ""), width = 7, height = 5.5, dpi = 300)
dev.off()






