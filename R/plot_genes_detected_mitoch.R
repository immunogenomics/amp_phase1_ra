setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

library(scales)
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
library(stringr)
library(data.table)

source("R/meta_colors.R")

# dat <- readRDS(file = paste("data/celseq_synovium_log2_postQC", ".rds", sep = ""))
# sc_meta <- readRDS(file = paste("data/celseq_synovium_meta", ".rds", sep = ""))
# all(colnames(dat) == sc_meta$cell_name)
# dim(dat)
# 
# all_fine <- readRDS(file = paste("data/all_cells_fine_cluster_label", ".rds", sep = ""))
# dat <- dat[, which(colnames(dat) %in% rownames(all_fine))]
# sc_meta <- sc_meta[which(sc_meta$cell_name %in% rownames(all_fine)),]
# all(colnames(dat) == sc_meta$cell_name)
# all_fine <- all_fine[order(match(rownames(all_fine), sc_meta$cell_name)), ]
# all(sc_meta$cell_name == rownames(all_fine))
# sc_meta$fine_cluster <- all_fine$fine_cluster
# dim(dat)

dat <- readRDS(file = paste("../data/celseq_synovium_log2_5452cells_paper", ".rds", sep = ""))
sc_meta <- readRDS(file = paste("../data/celseq_synovium_meta_5452cells_paper", ".rds", sep = ""))
all(colnames(dat) == sc_meta$cell_name)
sc_meta$fine_cluster <- sprintf("C-%s", str_replace(sc_meta$fine_cluster, "-", ""))

genes_mt  <- grep("^MT-", rownames(dat), value = TRUE)
genes_rp  <- grep("^RP[SL]", rownames(dat), value = TRUE)

counts_rowsums <- Matrix::rowSums(dat)
sum(counts_rowsums[genes_mt]) / sum(counts_rowsums)
sum(counts_rowsums[genes_rp]) / sum(counts_rowsums)

cutoff_mt <- 0.25
cutoff_genes <- 1000

gd <- apply(dat, 2, function(col) sum(col > 0))
sc_meta$gd <- gd


sc_meta$fine_cluster = factor(sc_meta$fine_cluster, 
                              levels=c("C-F1", "C-F2", "C-F3", "C-F4",
                                       "C-M1", "C-M2", "C-M3", "C-M4",
                                       "C-B1", "C-B2", "C-B3", "C-B4",
                                       "C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7"
                              ))
ggplot(
  data = sc_meta,
  mapping = aes(x = gd, y = percent_mt_molecules)
  ) +
  geom_point(size = 0.5) +
  geom_density_2d(bins = 4, color = "red") +
  geom_hline(yintercept = cutoff_mt) +
  geom_vline(xintercept = cutoff_genes) +
  scale_y_continuous(labels = percent) +
  scale_x_continuous(labels = function(x) ifelse(x > 0, paste(x / 1e3, "k"), 0)) +
  facet_wrap(~ fine_cluster, ncol = 4) +
  labs(
    x = "Genes Detected",
    y = "Percent of Molecules from MT",
    title = sprintf(
      "Percent of molecules assigned to %s mitochondrial genes",
      comma(length(genes_mt))
    )
  ) + 
  theme_bw(base_size = 20) +
  theme(    
    # axis.ticks = element_blank(), 
    # panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    legend.text = element_text(size = 19),
    plot.title = element_text(size=20)
  ) +
  # scale_y_discrete(breaks=NULL)
  coord_cartesian(
    xlim=c(0, 8000),
    ylim=c(0, 0.3)
  ) 
ggsave(file = paste("goodcells_per_cluster", ".pdf", sep = ""), width = 13, height = 12)
dev.off()


sc_meta$good <- rep(1, nrow(sc_meta))
sc_meta$good[which(sc_meta$gd < 1000 | sc_meta$percent_mt_molecules > 0.25)] <- 0
sum_bad <- table(sc_meta$fine_cluster, sc_meta$good)[,1]
sum_good <- table(sc_meta$fine_cluster, sc_meta$good)[,2]
per <- percent(sum_bad/(sum_bad + sum_good))
clus_good <- data.frame(
  sum_good = sum_good,
  sum_bad = sum_bad,
  per_bad = per
)


# ---------------------------------
dat_percent <- sc_meta %>%
  group_by(fine_cluster) %>%
  summarise(sum_good = sum(good == 0),
            sum_all = sum(good)) 
  # summarise(percent = sum(good == 0) / length(good) * 100)

ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = sum_all, fill = fine_cluster),
    color = "black", size = 0.15, width = 0.8
    # fill = "grey60"
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(dat_percent$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  labs(y = "Cells") +
  coord_flip() 


