#' ---
#' title: "Single-cell RNA-seq clustering"
#' Integrative pipeline of integrating bulk RNA-seq with single-cell RNA-seq and then unbiased clustering
#' author: "Fan Zhang"
#' date: "2018-03-19"
#' ---
#' 
setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

library(devtools)
library(igraph)
library(pheatmap)
library(Rtsne)
library(viridis)
library(RColorBrewer)
library(loe)
library(limma)
library(pbapply)
library(vegan)
library(cluster)
library(dplyr)
library(ggplot2)
library(CCA)
library(superheat)
library(biomaRt)
require(gdata)

source("meta_colors.R")
source("pure_functions.R")

# -------
# Main steps
# 1) We first selected the highly variable genes from both sides that we 
#    can keep the data set specific variation; 
# 2) Based on the shared highly variable genes, we integrate single-cell RNA-seq 
#    with bulk RNA-seq by CCA (González et al. 2008); 
# 3) We calculate the cell-to-cell similarity matrix on top 10 CCA canonical 
#    dimensions based on Euclidean distance;
# 4) We build up a K-nearest neighbors (KNN) graph based on the cell-to-cell similarity matrix; 
#    We then convert the KNN neighbor relation matrix into an adjacency matrix; 
# 5) We cluster the cells using Infomap algorithm (cluster_infomap function from igraph package) 
#    for community detection to decompose the cell-to-cell adjacency matrix into major modules; 
# 6) We then construct a low dimensional embedding using tSNE based on 
#    the cell-to-cell distance matrix with cell clusters labeled; 
# 7) We identify and prioritize significantly differential expressed (DE) genes for 
#    each distinct cluster using AUC and Wilcox test; 
# 8) For pathway analysis, we downloaded gene sets from Gene Ontology (GO) terms on April 2016. 
#    We also use the immunologic signatures from 4,872 hallmark gene sets from MSigDB (Liberzon et al. 2015)
#    to test enrichment of all the tested DE genes sorted by decreased AUC scores 
#    for each cluster using liger.
#
# NOTE: 
# 1) coarse clustering and fine clustering follow the same pipeline
#    For fine clustering, please run the piepeline on cells from one cell type 
# 2) Since graph-based clustering is based on random walk (a stochastic process), the results depand on nk (KNN step), 
#    and also the number dimensions ncc (CCA step), and the number of selected features as well.
# 
# -------


# Read the post-QC single-cell RNA-seq expression data and meta data
dat <- readRDS(file = paste("../data/celseq_synovium_log2_postQC", ".rds", sep = ""))
sc_meta <- readRDS(file = paste("../data/celseq_synovium_meta", ".rds", sep = ""))
all(colnames(dat) == sc_meta$cell_name)
dim(dat)

# all_fine <- readRDS(file = paste("data/all_cells_fine_cluster_label", ".rds", sep = "")) 
# dat <- dat[, which(colnames(dat) %in% rownames(all_fine))]
# sc_meta <- sc_meta[which(sc_meta$cell_name %in% rownames(all_fine)),]
# all(colnames(dat) == sc_meta$cell_name)

file_mean_sd <- "../data/celseq_synovium_log2cpm_mean_sd_test.rds"
if (!file.exists(file_mean_sd)) {
  dat_mean_sd <- data.frame(
    mean  = Matrix::rowMeans(dat),
    sd    = apply(dat, 1, sd),
    count = apply(dat, 1, function(x) sum(x > 0))
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
top_percent <- 0.7
mat_idx <- with(
  dat_mean_sd,
  mean > quantile(mean, top_percent) &
    sd > quantile(sd, top_percent)
)

# 1) Select the most highly variable genes based on mat_idx
sc_samples <- dat[mat_idx,]
colnames(sc_meta)[3] <- "Donor.ID"
all(sc_meta$cell_name == colnames(sc_samples))
dim(sc_samples)


# -----------------------------------------------------------------
# Load post-QC bulk RNA-seq data
log2tpm <- readRDS("../data/filtered_log2tpm_lowinput_phase_1.rds")
bulk_meta <- readRDS("../data/filtered_meta_lowinput_phase_1.rds")
all(colnames(log2tpm) == bulk_meta$Sample.ID)

file_mean_sd <- "data/celseq_synovium_log2tpm_mean_sd.rds"
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

bulk_samples <- log2tpm[mat_idx,]
all(bulk_meta$cell_name == colnames(bulk_samples))
dim(bulk_samples)

# 1) We first selected the highly variable genes from both sides 
# Take the union genes from both sides
genes <- union(rownames(bulk_samples), rownames(sc_samples))
length(genes)
bulk_samples_new <- log2tpm[which(rownames(log2tpm) %in% genes),]
sc_samples_new <- dat[which(rownames(dat) %in% genes),]

genes <- intersect(rownames(bulk_samples_new), rownames(sc_samples_new))
length(genes)
bulk_samples_new <- bulk_samples_new[which(rownames(bulk_samples_new) %in% genes),]
sc_samples_new <- sc_samples_new[which(rownames(sc_samples_new) %in% genes),]
sc_samples_new <- sc_samples_new[ order(match(rownames(sc_samples_new), rownames(bulk_samples_new))), ]
all(rownames(sc_samples_new) == rownames(bulk_samples_new))
dim(sc_samples_new)
dim(bulk_samples_new)


# -----------------------------------------------------------------
# 2) CCA: using the same genes for bulk and single-cell, different sizes of samples and cells
data_1 <- scale(bulk_samples_new)
data_2 <- scale(sc_samples_new)

# Do Classical CCA: bulk + single-cell (using the same genes)
res = cc(data_1, data_2)
dim(res$scores$corr.Y.xscores)

# A regularized step based on leave-one-out cross-validation process is needed when 
# the number of single cells is greater than the number of shared genes
# res.regul <- estim.regul(data_1, data_2, plt = TRUE) 
# res.rcc <- rcc(data_1, data_2, res.regul$lambda1, res.regul$lambda2)

saveRDS(res, file = paste("allcelltypes_cca_", nrow(bulk_samples_new), "_genes", ".rds", sep = ""))
# res <- readRDS("../data/allcelltypes_cca_7465_genes.rds")

# Plot the canonical variates
barplot(res$cor, xlab = "Dimension", ylab = "Canonical correlations", ylim = c(0,1))

# Combine top 30 CCA dimensions with the meta data
cc_sc <- data.frame(
  res$scores$corr.Y.xscores[, c(1:30)],
  disease = sc_meta$disease,
  plate = sc_meta$plate,
  type = sc_meta$type,
  cell_name = sc_meta$cell_name 
)
dim(cc_sc)
all(cc_sc$cell_name == sc_meta$cell_name)


# -----------------------------------------------------------------
# 3) Calculate Euclidean distance on top 10 CCs
ncc <- 10
top_cc <- cc_sc[, 1:ncc]
dat_dist <- dist(as.matrix(top_cc), method = "euclidean")
dat_dist <- as.matrix(dat_dist)
dim(dat_dist)


# -----------------------------------------------------------------
# 4) Make a KNN graph based on the Euclidean distance
# Large number of nk will result in less number of clusters
nk = 300
dat_knn <- make.kNNG(dat_dist, k=nk) 
dim(dat_knn)

# Convert Euclidean distance to adjacency matrix
g_adj  <- graph.adjacency(dat_knn, mode="undirected", weighted=TRUE)
# g <- graph_from_adjacency_matrix(g_adj, weighted=FALSE)


# -----------------------------------------------------------------
# 5) Community detection on adjacency matirx using cluster_infomap
# cluster_louvain is an alternative algorithm for community detection
imc <- cluster_infomap(g_adj)
# imc <- cluster_louvain(g_adj) 
membership(imc)
communities(imc)
cc_sc$cluster <- as.factor(membership(imc))

# Visualize the established graph (This step is too slow for plotting)
# plot(imc, g_adj, vertex.label=NA, vertex.size = 5)
# dev.copy(png,file = paste("plot_infomap", ".png", sep = ""),width=12,height=12,units="in",res=300)
# dev.off()


# -----------------------------------------------------------------
# 6) Visualize the identified clusters of cells in tSNE based on the distance matrix
per <- 40
tsne1 <- Rtsne(
  X = dat_dist,
  is_distance = TRUE,
  dims = 2,
  perplexity = per,
  theta = 0.8,
  pca = TRUE
)
cc_sc$T1 <- tsne1$Y[,1]
cc_sc$T2 <- tsne1$Y[,2]

# saveRDS(cc_sc, paste("all_infomap_", "cc_", nk, "nk_", per, "per", ".rds", sep = ""))


# tSNE plot
ggplot() +
  geom_point(
    data = cc_sc,
    mapping = aes_string(x = "T1", y = "T2", fill = "cluster"),
    # mapping = aes_string(x = "T1", y = "T2", fill = "type"),
    # mapping = aes_string(x = "T1", y = "T2", fill = "plate"),
    size = 2.5, stroke = 0.05, shape = 21
  ) +
  scale_fill_manual(values = meta_colors$cluster_all, name = "") +
  # scale_fill_manual(values = meta_colors$cluster_all_pair, name = "") +
  # scale_fill_manual(values = meta_colors$type, name = "") +
  # scale_fill_manual(values = meta_colors$plate, name = "") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 50) +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
    #axis.text = element_text(size = 35),
    #axis.text.y = element_text(size = 35)
  ) 
ggsave(file = paste("all_cluster_", "cc_", nk, "k_", per, "perp", ".png", sep = ""), width = 12.9, height = 9, dpi = 300)
dev.off()



# marker gene plot
gene <- "CD19"
gene_sc <- dat[which(rownames(dat) == gene),]
cc_sc$gene <- as.numeric(gene_sc)

ggplot() +
  geom_point(
    data = cc_sc,
    mapping = aes_string(x = "T1", y = "T2", fill = "gene"),
    size = 3, stroke = 0.05, shape = 21
  ) +
  # scale_fill_viridis(
  #   option = "viridis",
  #   name = bquote("Log"[2]~"(CPM)")
  # ) +
  scale_fill_gradientn(
    colours = colorRampPalette(RColorBrewer::brewer.pal(8, "Greens"))(10),
    name = ""
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
  theme_bw(base_size = 40) +
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(), 
    panel.grid = element_blank()
  ) 
ggsave(file = paste("all_", gene, ".png", sep = ""), width = 11, height = 9, dpi = 200)
dev.off()


# -----------------------------------------------------------------
# Feature selection for highly variable genes
file_mean_sd <- "celseq_synovium_log2cpm_mean_sd.rds"
dat_mean_sd <- readRDS(file_mean_sd)

top_percent <- 0.8
mat_idx <- with(
  dat_mean_sd,
  mean > quantile(mean, top_percent) &
    sd > quantile(sd, top_percent)
)
variable_genes <- rownames(dat[mat_idx,])
length(variable_genes)

# # Considering remove mitochondrial genes
# mito.genes_1 <- grep(pattern = "^MT-", x = variable_genes, value = TRUE)
# mito.genes_2 <- grep(pattern = "^MTRNR", x = variable_genes, value = TRUE)
# mito.genes <- c(mito.genes_1, mito.genes_2)
# variable_genes <- variable_genes[-which(variable_genes %in% mito.genes)]



# -----------------------------------------------------------------
# 7) Find differential expressed genes based on AUC score
DifferentialAUC <- function(x, y) {
  prediction.use <- ROCR::prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- ROCR::performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}

# To save the time, here we only test the most highly variable genes
gene_res <- do.call(rbind, pbapply::pblapply(X = variable_genes, cl = 5, FUN = function(gene_name) {
  gene <- as.numeric(dat[gene_name,])
  do.call(rbind, lapply(1:length(table(cc_sc$cluster)), function(i) {
    x <- gene[cc_sc$cluster == i]
    y <- gene[cc_sc$cluster != i]
    data.frame(
      "gene" = gene_name,
      "cluster" = i,
      "wilcox" = wilcox.test(x, y, alternative = "greater")$p.value,
      "auc" = DifferentialAUC(x, y),
      # "auc2" = DifferentialAUC(x > 0, y > 0),
      "mean" = mean(x)
    )
  }))
}))
saveRDS(gene_res, file = paste("cluster_markers_res", "_cc_",  nk, "knn_", per, "per", ".rds", sep = ""))


# Sorted genes based on AUC score
gene_res_best <- do.call(rbind, lapply(1:length(table(cc_sc$cluster)), function(i) {
  x <- gene_res %>% group_by(gene) %>%
    summarise(
      wilcox = wilcox[cluster == i],
      auc = auc[cluster == i],
      # Sort by mean diff between cluster and next biggest cluster
      mean_diff = mean[cluster == i] - sort(mean[cluster != i], decreasing = TRUE)[1]
      # mean_diff_1 (log2FC): mean_diff = mean[cluster == i] - sort(mean[cluster != i], decreasing = TRUE)[1]
      # mean_diff_2 (log2FC)
      # Effect size: Z-score
      # t.test p-value
    )
  x <- x[order(x$auc, decreasing = TRUE),]
  x <- head(x, 200)
  x$cluster <- i
  x
}))
View(gene_res_best)

# Export the cluster marker genes
saveRDS(gene_res_best, file = paste("cluster_markers_top200", "_cc_",  nk, "knn_", per, "per", ".rds", sep = ""))


# -----------------------------------------------------------------
# 8) Gene enrichment for distinct clusters using liger based on MSIGDB

# Load Gene sets (translate GENE NAME to ENSEMBLE ID)
match <- read.table(file = "../data/gene_ID_name_match_gencodeV24.txt.gz", header = T, stringsAsFactors = F )
ids <- sapply(match$GENE_ID, function(gene){
  strsplit(gene, "\\.")[[1]][1]
})
rownames(match) <- ids
match <- subset(match, select = -c(GENE_ID) )
match$GENE_ID <- rownames(match)

load("data/gene_sets.rda")
symbol_to_ensembl <- unlist(split(match$GENE_ID, match$GENE_NAME))

# Take all the variable_genes (sorted on the auc value) for each cluster
c <- gene_res_best[which(gene_res_best$cluster == "1"),]$auc
names(c) <- gene_res_best[which(gene_res_best$cluster == "1"),]$gene
names(c)[1:4]

names(c) <- symbol_to_ensembl[names(c)]
c <- c[!is.na(names(c))]
length(c)


lig1 <- liger::bulk.gsea(
  values = c,
  set.list = MSIGDB_C7_ENSEMBL,
  n.rand = 1e5,
  mc.cores = 4,
  rank = TRUE
)
lig1 <- lig1[order(lig1$q.val),]
head(lig1, 50)
testName = "C7_cluster1"
write.table(lig1, file = paste("liger_", testName, ".txt", sep = ""), quote = F, sep = "\t")



# -----------------------------------------------------------------
# Estimate the optimal number of clusters (k) using Silhouette 
# dat_dist <- dat_dist[which(rownames(dat_dist) %in% rownames(all_fine)),which(rownames(dat_dist) %in% rownames(all_fine))]
# all(rownames(dat_dist) == colnames(dat_dist))
# cc_sc <- cc_sc[which(cc_sc$cell_name %in% rownames(all_fine)),]

dat_dist <- dat_dist[ order(match(rownames(dat_dist), rownames(all_fine))), ]
all(rownames(dat_dist) == rownames(all_fine))
table(all_fine$clus_number)

sil = silhouette(as.numeric(all_fine$clus_number), dat_dist)
pdf('silhouette_plot_all_fine_k19_newcolor.pdf', height = 10, width = 8)
plot(sil,
     col = meta_colors$fine_cluster, 
     xlab = "Silhouette width",
     ylab = "Cells",
     main = "Silhouette plot of scRNA-seq clusters and \ncell-to-cell similarity matrix",
     cex.lab = 1,
     cex = 1,
     cex.axis = 1
     )
# abline(v = 0.22, lty = 2)
dev.off()


# Calculate silhouette width for many k using PAM
# PAM it is more robust (compared to Kmeans) because it minimizes a sum of dissimilarities 
# instead of a sum of squared euclidean distances;
# Use Gene expression data matrix to caluculate euclidean distance is too slow to get results
# Here I used CCs (the same one that used for graph-based clustering) to calculate distance.

# Estimate the optimal value for k by run PMA on the distance matrix with multiple k
sil_width <- c(NA)
for(i in 2:15){
  pam_fit <- pam(dat_dist,
                 diss = TRUE,
                 k = i)
  sil_width[i] <- pam_fit$silinfo$avg.width
}
# saveRDS(sil_width, "sil_width.rds")
# sil_width <- readRDS("sil_width.rds")

# Plot sihouette width (higher is better)
plot(sil_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")


