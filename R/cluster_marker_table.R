# title: "Generating single-cell RNA-seq cluster marker table"

library(pbapply)
library(parallel)
library(forcats)
library(stringr)
library(data.table)

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/")

auroc <- function(score, cls) {
  n1 <- sum(!cls)
  n2 <- sum(cls)
  U <- sum(rank(score)[!cls]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

# cca <- readRDS("data/allcelltypes_cca_7465_genes.rds")

log2cpm <- readRDS("data/celseq_synovium_log2_postQC.rds")

meta <- readRDS("data/celseq_synovium_meta.rds")
meta <- meta[order(meta$cell_name),]
meta$cell_name <- as.character(meta$cell_name)

x <- readRDS("/data/srlab/slowikow/github.com/immunogenomics/ampviewer/data/all_cells_fine_cluster_label.rds")
x <- x[order(rownames(x)),]

meta <- meta[meta$cell_name %in% rownames(x),]

stopifnot(all(rownames(x) == meta$cell_name))

log2cpm <- log2cpm[,meta$cell_name]

stopifnot(all(colnames(log2cpm) == meta$cell_name))

# tSNE coordinates
meta$T1 <- x$T1
meta$T2 <- x$T2

# cluster assignment
meta$cluster <- x$fine_cluster
meta$cluster <- sprintf("C-%s", str_replace(meta$cluster, "-", ""))

# Count how many cells express each gene.
log2cpm_ncells <- apply(log2cpm, 1, function(row) {
  sum(row > 0)
})

# Exclude genes detected in fewer than 10 cells.
log2cpm <- log2cpm[log2cpm_ncells >= 10,]

cell_clusters <- sort(unique(meta$cluster))

get_markers <- function(log2cpm, cell_clusters) {
  # Compute statistics for each cluster.
  dat_marker <- rbindlist(pblapply(
    X = rownames(log2cpm),
    cl = 20,
    FUN = function(gene_name) {
      gene      <- as.numeric(log2cpm[gene_name,])
      gene_mean <- mean(gene)
      gene_sd   <- sd(gene)
      rbindlist(lapply(unique(cell_clusters), function(cell_cluster) {
        ix <- cell_clusters == cell_cluster
        x <- gene[ix]
        x_mean <- mean(x)
        pct_nonzero <- sum(x > 0) / length(x)
        y <- gene[!ix]
        pct_nonzero_other <- sum(y > 0) / length(y)
        data.frame(
          "gene"          = gene_name,
          "cluster"       = cell_cluster,
          "wilcox_pvalue" = wilcox.test(x, y, alternative = "greater")$p.value,
          "auc"           = auroc(gene, ix),
          "pct_nonzero"   = pct_nonzero,
          "pct_nonzero_other" = pct_nonzero_other,
          "cluster_mean"  = x_mean,
          "log2FC"        = x_mean - mean(y),
          "zscore"        = (x_mean - gene_mean) / gene_sd
        )
      }))
    }
  ))
}

dat_marker <- get_markers(
  log2cpm,
  meta$cluster
)
dat_marker$cell_type <- "All cells"

cell_types <- unique(as.character(meta$type))

cell_type_markers <- rbindlist(lapply(cell_types, function(cell_type) {
  ix <- meta$type == cell_type
  d <- get_markers(
    log2cpm[,ix],
    meta$cluster[ix]
  )
  d$cell_type <- cell_type
  return(d)
}))

dat_marker <- rbindlist(list(dat_marker, cell_type_markers))

x <- dat_marker[
  cell_type != "All cells",
  .SD[order(-.SD$auc),][1:3],
  by = cluster
]

saveRDS(
  object = y,
  file = "data/cluster_marker_table.rds"
)

