setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R")

library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)

risk_genes_okada <- read.csv("../../../HMS/amp/results/2017_01_17_Downstream_analysis/Okada2014_riskGenesRA.csv")
ra_genes_tiffany <- load('../data/curated_regions_OR_RA_RAneg_SLE_PsA_AS_Gout.RData')

exp <- readRDS("../data/celseq_synovium_log2_5452cells_paper.rds")
meta <- readRDS("../data/celseq_synovium_meta_5452cells_paper.rds")

meta <- meta[-which(meta$fine_cluster == "T-1"), ]
exp <- exp[, which(colnames(exp) %in% meta$cell_name)]
all(meta$cell_name == colnames(exp))

# RA risk genes from okada
exp <- exp[which(rownames(exp) %in% risk_genes_okada$Gene),]

# RA risk genes from Tiffany summarized GWAS and genes from HM and Rachel
as.character(curated_regions_OR_RA_RAneg_SLE_PsA_AS_Gout$Genes)
exp <- exp[which(rownames(exp) %in% as.character(curated_regions_OR_RA_RAneg_SLE_PsA_AS_Gout$Genes)),]

dim(exp)

combine <- cbind.data.frame(t(exp), meta)
combine <- combine[, c(1:nrow(exp), which(colnames(combine) == "fine_cluster"))]

# Take the average of expression for each gene each cluster

mean_gene <- combine %>% 
  group_by(fine_cluster) %>%
  summarise_all("mean")

mean_gene <- as.data.frame(mean_gene)
rownames(mean_gene) <- mean_gene$fine_cluster
mean_gene <- mean_gene[, -1]
mean_gene <- t(mean_gene)
# mean_gene <- as.data.frame(mean_gene)
dim(mean_gene)

# Check the mean and sd of each gene across all the clusters
plot(apply( mean_gene, 1, mean ), apply( mean_gene, 1, sd )) 
gene_plot <- mean_gene[which(apply( mean_gene, 1, max ) > 1.5 &
                                apply( mean_gene, 1, sd ) > 1 ), ]
dim(gene_plot)

pdf("RA_riskgenes_sc_clusters.pdf", width=5, height=10)
pheatmap(
  mat = gene_plot,
  border_color = "white",
  # color  = viridis(10),
  color = colorRampPalette(brewer.pal(n=9, name = "Blues"))(9),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  # cellwidth = 5,
  # cellheight = 5,
  fontsize = 10,
  fontsize_row = 8,
  scale = 'none',
  legend = TRUE
)
dev.off()

