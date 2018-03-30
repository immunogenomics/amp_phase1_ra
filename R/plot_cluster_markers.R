#' ---
#' title: "Plot DE markers for each single-cell RNA-seq cluster"
#' author: "Fan Zhang"
#' date: "2018-03-30"
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R")

library(ggplot2)
library(ggrepel)
library(gdata) 

source("meta_colors.R")

dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table[1:4,]
table(dat_table$cluster)

dat_high_pct <- dat_table[which(dat_table$pct_nonzero > 0.66),]
table(dat_high_pct$cluster)
table(dat_high_pct$cell_type)

dat_high_pct$cluster = factor(dat_high_pct$cluster, 
                              levels=c("C-F1", "C-F2", "C-F3", "C-F4",
                                       "C-M1", "C-M2", "C-M3", "C-M4",
                                       "C-B1", "C-B2", "C-B3", "C-B4",
                                       "C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7"
                                       ))


ggplot() + 
  geom_point(
    data=dat_high_pct,
    aes(x=auc, 
        y=log2FC,
        fill = cluster
        ),
    size = 1
  ) +
  geom_vline(xintercept = 0.6, linetype="solid", 
             color = "blue", size=0.2) +
  geom_hline(yintercept = 2, linetype="solid", 
             color = "blue", size=0.2) +
  facet_wrap(~ cluster, ncol = 4) +
  labs(
    x = "AUROC",
    y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 20) +
  theme(    
    legend.position = "none",
    # axis.ticks = element_blank(), 
    # panel.grid = element_blank(),
    axis.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 30)
    ) +
  # scale_y_discrete(breaks=NULL)
  coord_cartesian(
    xlim=c(0, 1), 
    ylim=c(-5, 10)
  ) 
ggsave(file = paste("auc_log2FC", ".pdf", sep = ""),
       width = 11, height = 10, dpi = 150)
dev.off()


ggplot() + 
  geom_point(
    data=dat_high_pct,
    aes(x=auc, 
        y= -log10(wilcox_pvalue)
        # y= -log10(ttest_pvalue)
    ),
    size = 1
  ) +
  geom_vline(xintercept = 0.65, linetype="solid", 
             color = "blue", size=0.2) +
  geom_hline(yintercept = 20, linetype="solid", 
              color = "blue", size=0.2) +
  facet_wrap(~ cluster, ncol = 4) +
  labs(
    x = "AUROC",
    y = "-Log10(wilcox P)"
  ) +
  theme_bw(base_size = 20) +
  theme(    
    legend.position = "none",
    # axis.ticks = element_blank(), 
    # panel.grid = element_blank(),
    axis.text = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20, angle = 30)
  ) +
  # scale_y_discrete(breaks=NULL)
  coord_cartesian(
    xlim=c(0, 1),
    ylim=c(0, 300)
  ) 
ggsave(file = paste("auc_wilcox", ".pdf", sep = ""),
       width = 12, height = 11, dpi = 150)
dev.off()

