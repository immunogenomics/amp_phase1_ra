#' ---
#' title: "Plot overlapped samples for all 4 data types"
#' author: "Fan Zhang"
#' date: "2018-03-20"
#' Figure 1C
#' ---


setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/code_visualization")

library(ggplot2)
library(reshape2)
library(dplyr)
require(gdata)

source("../R/pure_functions.R")
source("../R/meta_colors.R")

all_postQC <- readRDS("../data/all_postQC_for_barplot.rds")

inflam_label <- read.xls("../data-raw/postQC_all_samples.xlsx")
colnames(inflam_label)[1] <- "Sample"

dat_merge <- merge(all_postQC, inflam_label, by = "Sample")

ggplot() +
  geom_bar(
    data = dat_merge,
    mapping = aes(x = reorder(Sample, value), weight = value, fill = Type),
    colour="black"
  ) +
  # facet_grid(variable ~ DiseaseAssay, scales = "free", space = "free_x") +
  facet_grid(variable ~ lym_25, scales = "free", space = "free_x") +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = meta_colors$type, name = "Cell Type") +
  theme(
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x=element_text(angle=45,hjust=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
       # x = "Donors",
       # y = "Number of Cells (Log scale)"
       x = NULL, 
       y = NULL
       )
ggsave(file = paste("barplot_all_postQC", ".pdf", sep = ""), width = 18, height = 8, dpi = 200)
dev.off()



