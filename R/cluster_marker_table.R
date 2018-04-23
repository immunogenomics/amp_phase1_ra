library(pbapply)
library(parallel)
library(forcats)
library(stringr)
library(data.table)

auroc <- function(score, cls) {
  n1 <- sum(!cls)
  n2 <- sum(cls)
  U <- sum(rank(score)[!cls]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

# cca <- readRDS("data/allcelltypes_cca_7465_genes.rds")

log2cpm <- readRDS("data/celseq_synovium_log2_5452cells_paper.rds")

meta <- readRDS("data/celseq_synovium_meta_5452cells_paper.rds")
# meta <- meta[order(meta$cell_name),]
meta$cell_name <- as.character(meta$cell_name)

stopifnot(all(colnames(log2cpm) == meta$cell_name))

# cluster assignment
meta$cluster <- sprintf("C-%s", str_replace(meta$fine_cluster, "-", ""))

# Count how many cells express each gene.
log2cpm_ncells <- apply(log2cpm, 1, function(row) {
  sum(row > 0)
})

# Exclude genes detected in fewer than 10 cells.
log2cpm <- log2cpm[log2cpm_ncells >= 10,]
dim(log2cpm)

unique_cell_clusters <- sort(unique(meta$cluster))

cell_clusters <- meta$cluster

get_markers <- function(log2cpm, cell_clusters) {
  # Compute statistics for each cluster.
  dat_marker <- rbindlist(pblapply(
    X = rownames(log2cpm),
    cl = 20,
    FUN = function(gene_name) {
      gene      <- as.numeric(log2cpm[gene_name,])
      rbindlist(lapply(unique(cell_clusters), function(cell_cluster) {
        ix <- cell_clusters == cell_cluster
        x <- gene[ix]
        x_mean <- mean(x)
        x_sd   <- sd(x)
        x_pct_nonzero <- sum(x > 0) / length(x)
        y <- gene[!ix]
        y_mean <- mean(y)
        y_sd   <- sd(y)
        y_pct_nonzero <- sum(y > 0) / length(y)
        test_w <- wilcox.test(x, y, alternative = "two.sided")
        test_t <- t.test(x, y, alternative = "two.sided")
        data.frame(
          "gene"              = gene_name,
          "cluster"           = cell_cluster,
          "wilcox_pvalue"     = test_w$p.value,
          "ttest_pvalue"      = test_t$p.value,
          "auc"               = auroc(gene, ix),
          "pct_nonzero"       = x_pct_nonzero,
          "pct_nonzero_other" = y_pct_nonzero,
          "log2FC"            = x_mean - y_mean,
          "mean"              = x_mean,
          "sd"                = x_sd,
          "mean_other"        = y_mean,
          "sd_other"          = y_sd
        )
      }))
    }
  ))
  # Check if the mean is highest in this cluster.
  dat_marker[
    ,
    mean_highest := mean >= max(mean),
    by = gene
  ]
  dat_marker[
    ,
    pct_nonzero_highest := pct_nonzero >= max(pct_nonzero),
    by = gene
  ]
  return(dat_marker)
}

# Test each cluster against all cells.
dat_marker <- get_markers(
  log2cpm,
  meta$cluster
)
dat_marker$cell_type <- "All cells"

# Test each cluster against other clusters of the same cell type.
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

# Put the results together.
dat_marker <- rbindlist(list(dat_marker, cell_type_markers))

x <- dat_marker[
  cell_type != "All cells",
  .SD[order(-.SD$auc),][1:10],
  by = cluster
]

saveRDS(
  object = dat_marker,
  file = "data/cluster_marker_table.rds"
)

