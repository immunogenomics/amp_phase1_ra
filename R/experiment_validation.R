# Analysis for experimental validation
# Fan Zhang
# 2018-10-23

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
source("../amp_phase1_ra_viewer/R/install-packages.R")
source("R/meta_colors.R")

# -------------------------
# T cell experimental data from Helena
dat <- read.xls("../../HMS/amp/resource_paper/amp_paper_revision/experiment_validation/AMP Phase 1 revision data for Fan.xlsx")
dat_1 <- dat[which(dat$Experiment == "TNF production"), c(1:3)]
dat_1$group <- c(seq(1,7), seq(1,7))
dat_1$value <- dat_1$value / 100
dat_1$cell_type <- as.character(dat_1$cell_type)
dat_1$cell_type[which(dat_1$cell_type == "CD4")] <- "CD4 T cells"
dat_1$cell_type[which(dat_1$cell_type == "CD8")] <- "CD8 T cells"

t.test(dat_1$value[which(dat_1$cell_type == "CD8 T cells")], dat_1$value[which(dat_1$cell_type == "CD4 T cells")], alternative = "greater")

ggplot(
  data = dat_1,
  aes(x = cell_type, y = value, group = group)
  ) +
  geom_line() +
  geom_point(
    aes(x = cell_type, y = value),
    size = 3, stroke = 0.2, shape = 19) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = pretty_breaks()) +
  labs(
    x = "",
    y = expression(paste("Percent of IFN", gamma^"+", " cells")),
    # title = bquote("IFN"~gamma~"production by \n synovial T cells")
    title = expression(paste("IFN", gamma, " production")),
    subtitle = "by synovial T cells"
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(color="black", size=20)
  ) 
ggsave(file = paste("IFN", ".pdf", sep = ""),
       width = 4, height = 6, dpi = 300)
dev.off()


dat_2 <- dat[which(dat$Experiment == "IFN production"), c(1:3)]
dat_2$group <- c(seq(1,7), seq(1,7))
dat_2$value <- dat_2$value / 100
dat_2$cell_type <- as.character(dat_2$cell_type)
dat_2$cell_type[which(dat_2$cell_type == "CD4")] <- "CD4 T cells"
dat_2$cell_type[which(dat_2$cell_type == "CD8")] <- "CD8 T cells"

ggplot(
  data = dat_2,
  aes(x = cell_type, y = value, group = group)
  ) +
  geom_line() +
  geom_point(
    aes(x = cell_type, y = value),
    size = 3, stroke = 0.2, shape = 19) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = pretty_breaks()) +
  labs(
    x = "",
    y = expression(paste("Percent of ", "TNF"^"+", " cells")),
    title = "TNF production",
    subtitle = expression(paste("by synovial T cells"))
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=20), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(color="black", size=20)
  ) 
ggsave(file = paste("TNF", ".pdf", sep = ""),
       width = 4, height = 6, dpi = 300)
dev.off()


dat_3 <- dat[c(1:3), c(4:ncol(dat))]
colnames(dat_3) <- c("GzmB+", "GzmB+GzmK+", "GzmK+", "Gzm-", "Experiment_Gzm")
dat_3 <- melt(dat_3, na.rm = TRUE, id="Experiment_Gzm")
dat_3$value <- dat_3$value / 100
dat_median <- dat_3 %>% group_by(variable) %>% summarise(median = median(value))

ggplot(
  data = dat_3,
  aes(x = variable, y = value)
  ) +
  # geom_quasirandom(
  #   shape = 19, size = 3, stroke = 0.2
  # ) +
  stat_summary(
    fun.y = median, fun.ymin = median, fun.ymax = median,
    geom = "crossbar", width = 0.5
  ) +
  # geom_boxplot() +
  geom_point(
    aes(x = variable, y = value),
    size = 3, stroke = 0.2, shape = 19) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                   breaks = pretty_breaks()) +
  labs(
    x = "",
    y = expression(paste("Percent among HLA-", "DR"^"+", " CD8 T cells")),
    title = "Granzyme expression by ",
    subtitle = expression(paste("HLA-", "DR"^"+", " CD8 T cells"))
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black", size=18), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(color="black", size=18)
  ) 
ggsave(file = paste("Gzm", ".pdf", sep = ""),
       width = 4, height = 6, dpi = 300)
dev.off()
