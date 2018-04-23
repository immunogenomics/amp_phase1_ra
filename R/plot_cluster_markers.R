#' ---
#' title: "Plot DE markers for each single-cell RNA-seq cluster"
#' author: "Fan Zhang"
#' date: "2018-03-30"
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R")

library(ggplot2)
library(ggrepel)
library(gdata) 
library(pbapply)
library(dplyr)

source("meta_colors.R")

dat_table <- readRDS("../data/cluster_marker_table.rds")
table(dat_table$cluster)
test <- dat_table$cell_type
test[which(test != "All cells")] <- "One cell type cells"
dat_table$test <- test

# "All cells" are used for comparing one cluster versus all (18) the other clusters
# dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
# dim(dat_table)

# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2)), ] 

dat_high_pct <- dat_table[which(dat_table$pct_nonzero > 0.6),]
table(dat_high_pct$cluster)
table(dat_high_pct$cell_type)

dat_high_pct$cluster = factor(dat_high_pct$cluster, 
                              levels=c("C-F1", "C-F2", "C-F3", "C-F4",
                                       "C-M1", "C-M2", "C-M3", "C-M4",
                                       "C-B1", "C-B2", "C-B3", "C-B4",
                                       "C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7"
                                       ))

# Plot AUC vs log2 FC
dat_plot <- dat_high_pct[which(dat_high_pct$test == "All cells"),]
dat_plot$sig <- rep(0, nrow(dat_plot))
dat_plot$sig[which(dat_plot$auc > 0.7 & dat_plot$log2FC > 1)] <- 1
dat_plot$sig <- as.factor(dat_plot$sig)

ggplot() + 
  geom_point(
    data=dat_plot,
    aes(x=auc, 
        y=log2FC,
        color = sig
    ),
    size = 0.5
  ) +
  scale_color_manual(values = c('grey', '#D95F02'), # c('#1B9E77', '#D95F02'),
                     name="",
                     # labels= c("One cluster vs all the other clusters",
                     #          "One cluster vs the other clusters within the same cell type")
                     labels= c("",
                               "Genes with AUC > 0.7, Log2(FC) > 1, and non-zero expression > 60%")
                     ) +
  geom_vline(xintercept = 0.7, linetype="solid", color = "black", size=0.2) +
  geom_hline(yintercept = 1, linetype="solid", color = "black", size=0.2) +
  facet_wrap(~ cluster, ncol = 4) +
  labs(
    x = "AUC",
    y = "Log2 (FC)",
    title = "Genes that percent of non-zero expression > 60% for each cluster"
  ) +
  theme_bw(base_size = 22) +
  theme(    
    # axis.ticks = element_blank(), 
    # panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 30),
    legend.text = element_text(size = 19),
    plot.title = element_text(size=22)
  ) 
  # coord_cartesian(
  #   xlim=c(0, 1),
  #   ylim=c(-5, 10)
  # ) 
ggsave(file = paste("auc_log2FC_test2_goodcells", ".pdf", sep = ""), width = 13, height = 12)
dev.off()


# Plot AUC vs wilcox p-value
ggplot() + 
  geom_point(
    data=dat_high_pct,
    aes(x=auc, 
        y= -log10(wilcox_pvalue),
        # y= -log10(ttest_pvalue),
        color = test
    ),
    size = 0.5
  ) +
  scale_color_manual(values = c('#1B9E77', '#D95F02'),
                     name="",
                     labels=c("One cluster vs all the other clusters",
                              "One cluster vs the other clusters within the same cell type")) +
  geom_vline(xintercept = 0.7, linetype="solid", color = "black", size=0.2) +
  geom_hline(yintercept = 20, linetype="solid", color = "black", size=0.2) +
  facet_wrap(~ cluster, ncol = 4) +
  labs(
    x = "AUC",
    y = "-Log10 (wilcox P)",
    title = "Genes that percent of non-zero expression > 60% for each cluster"
  ) +
  theme_bw(base_size = 20) +
  theme(    
    # axis.ticks = element_blank(), 
    # panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 30),
    legend.text = element_text(size = 19),
    plot.title = element_text(size=22)
  ) +
  # scale_y_discrete(breaks=NULL)
  coord_cartesian(
    xlim=c(0, 1),
    ylim=c(0, 300)
  ) 
ggsave(file = paste("auc_wilcox", ".pdf", sep = ""), width = 13, height = 12)
dev.off()



# Report # genes that pass four criteria:
# 1) pct_nonzero > 60%, 2) AUC > 0.7, 3) -log10(wilcox p) > 20, and 4) log2FC > 1

dat_table$wilcox_log <- -log10(dat_table$wilcox_pvalue)
dat_table$wilcox_log[which(dat_table$wilcox_log == "Inf")] <- 300

# good_markers <- dat_table[which(dat_table$pct_nonzero > 0.6 & dat_table$auc > 0.7 & dat_table$wilcox_log > 20 & dat_table$log2FC > 1),]
good_markers <- dat_table[which(dat_table$pct_nonzero > 0.6 & dat_table$auc > 0.7 & dat_table$log2FC > 1),]
good_markers$cluster = factor(good_markers$cluster, 
                              levels=c("C-F1", "C-F2", "C-F3", "C-F4",
                                       "C-M1", "C-M2", "C-M3", "C-M4",
                                       "C-B1", "C-B2", "C-B3", "C-B4",
                                       "C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7"
                              ))
table(good_markers$cluster)


# Export top 10 - 20 cluster marker genes 
x <- dat_table %>%
  group_by(cluster) %>%
  top_n(20, wt = auc)
dim(x)
x$cluster = factor(x$cluster, 
                  levels=c("C-F1", "C-F2", "C-F3", "C-F4",
                           "C-M1", "C-M2", "C-M3", "C-M4",
                           "C-B1", "C-B2", "C-B3", "C-B4",
                           "C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7"
                              ))
x <- x[order(x$auc, decreasing = TRUE),]
x <- x[order(x$cluster, decreasing = FALSE),]

x_save <- x[, c("gene", "cluster", "cell_type", "auc", "wilcox_pvalue", "ttest_pvalue", "pct_nonzero", "pct_nonzero_other", "log2FC")]

write.table(x_save, file = "top20_singlecell_cluster_markers.txt", quote = F, sep = "\t")


