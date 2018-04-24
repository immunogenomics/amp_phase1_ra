#' ---
#' title: "Integrative pipeline of integrating bulk RNA-seq with identified mass cytometry clusters"
#' author: "Fan Zhang"
#' date: "2018-04-24"
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

library(CCA)
library(nscancor)
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
library(pheatmap)
library(biomaRt)
library("pspearman")
require(gdata)

source("pure_functions.R")
source("meta_colors.R")

# -------------------
# Load the low input samples
log2tpm <- readRDS("../data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("../data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)

file_mean_sd <- "lowinput_log2tpm_mean_sd_fibro.rds"
if (!file.exists(file_mean_sd)) {
  dat_mean_sd <- data.frame(
    mean  = Matrix::rowMeans(log2tpm),
    sd    = apply(log2tpm, 1, sd),
    count = apply(log2tpm, 1, function(x) sum(x > 0))
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

top_percent <- 0.8
mat_idx <- with(
  dat_mean_sd,
  mean > quantile(mean, top_percent) &
    sd > quantile(sd, top_percent)
)

# Take fibroblast as an example
cell_type <- "Fibro"
bulk_samples_type <- log2tpm[, which(bulk_meta$Cell.type==cell_type)]
bulk_meta <- bulk_meta[which(bulk_meta$Cell.type==cell_type),]
bulk_samples <- bulk_samples_type[mat_idx,]
dim(bulk_samples)
dim(bulk_meta)
all(bulk_meta$Sample.ID == colnames(bulk_samples))

# Convert to: samples x genes
bulk_samples <- t(bulk_samples)
dim(bulk_samples)

# -------------------------------
# Load cytof clusters 
cytof = read.xls ("../data/AMP.cluster.sample.frequency.xlsx", sheet = 5, header = TRUE) 
temp <- gsub('[(.)]','-',colnames(cytof))
colnames(cytof) <- gsub('[(X)]','',temp)
rownames(cytof) <- cytof$cluster
cytof <- cytof[, -1]

# -------------------------------
# Take the intersected donors
inter_donor <- intersect(names(table(bulk_meta$Donor.ID)), colnames(cytof))
bulk_samples <- bulk_samples[which(bulk_meta$Donor.ID %in% inter_donor),]
bulk_meta <- bulk_meta[which(bulk_meta$Donor.ID %in% inter_donor),]
dim(bulk_samples)
dim(bulk_meta)
all(bulk_meta$Sample.ID == rownames(bulk_samples))
rownames(bulk_samples) <- bulk_meta$Donor.ID

cytof_samples <- cytof[, which(colnames(cytof) %in% inter_donor)]
cytof_samples <- t(cytof_samples)
colnames(cytof_samples) <- seq(1, ncol(cytof_samples))
dim(cytof_samples)

bulk_samples <- bulk_samples[ order(match(rownames(bulk_samples), rownames(cytof_samples))), ]
all(rownames(bulk_samples) == rownames(cytof_samples))
bulk_meta <- bulk_meta[ order(match(bulk_meta$Donor.ID, rownames(bulk_samples))),]
all(bulk_meta$Donor.ID == rownames(bulk_samples))
dim(bulk_samples)
dim(cytof_samples)

# -------------------------------
# Run regularized CCA 
data_1 <- scale(bulk_samples)
data_2 <- scale(cytof_samples)

res.regul <- estim.regul(data_1, data_2, plt = TRUE)
saveRDS(res.regul, file = paste(cell_type, "_regul_rcca_bulk", ncol(bulk_samples), "genes_cytof", ncol(cytof_samples), "clusters", ".rds", sep = ""))

res.rcc <- rcc(data_1, data_2, res.regul$lambda1, res.regul$lambda2)
saveRDS(res.rcc, file = paste(cell_type, "_res_rcca_bulk", ncol(bulk_samples), "genes_cytof", ncol(cytof_samples), "clusters", ".rds", sep = ""))


# # Note: add a Gaussian noise to the cytof matrix in case if some matrix is not positive semidefinite
#  Noisify <- function(data) {
# 	    if (is.vector(data)) {
# 		         noise <- runif(length(data), -1e-7, 1e-7)
#      noisified <- data + noise
#         } else {
# 		     length <- dim(data)[1] * dim(data)[2]
#           noise <- matrix(runif(length, -1e-7, 1e-7), dim(data)[1])
# 	       noisified <- data + noise
# 	     }
#    return(noisified)
#     }
# bulk_samples_noise <- Noisify(bulk_samples)
# cytof_samples_noise <- Noisify(cytof_samples)

