#' ---
#' title: "Plot cells from all single-cell RNA-seq clusters for each patient"
#' author: "Fan Zhang"
#' date: "2018-03-25"
#' Figure 2
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/")

source("R/meta_colors.R")

# Read post-QC single-cell RNA-seq meta data
meta <- readRDS("data/celseq_synovium_meta.rds")
meta <- meta[order(meta$cell_name),]
meta$cell_name <- as.character(meta$cell_name)
dim(meta)

# Read cluster labels 
x <- readRDS("/Users/fanzhang/Documents/GitHub/ampviewer/data/all_cells_fine_cluster_label.rds")
x <- x[order(rownames(x)),]

meta <- meta[meta$cell_name %in% rownames(x),]
dim(meta)

all(rownames(x) == meta$cell_name)
meta$fine_cluster <- x$fine_cluster
mat <- meta[, c("sample", "disease", "plate", "fine_cluster")]
non_samples <- names(which(table(mat$sample) == 0))

dat <- table(mat$fine_cluster, mat$sample)
dat <- as.data.frame(dat)
colnames(dat) <- c("fine_cluster", "sample", "freq")
dat <- dat[-which(dat$sample %in% non_samples),]
dat$sample <- as.character(dat$sample)

table(dat$sample)
dim(dat)

# Read Case.Control labels for OA, non-inflamed RA, and RA
inflam_label <- read.xls("data-raw/postQC_all_samples.xlsx")
inflam_label <- inflam_label[c(1:51),]
table(inflam_label$Case.Control)
inflam_label <- inflam_label[which(inflam_label$Patient %in% dat$sample), ]
inflam_label <- inflam_label[, c("Patient", "Disease", "Tissue.type", "Case.Control")]
colnames(inflam_label)[1] <- "sample"
plot_final <- merge(dat, inflam_label, by = "sample")
dim(plot_final)

plot_final$Case_Control_order = factor(plot_final$Case.Control, levels=c('OA','non-inflamed RA','inflamed RA'))

# saveRDS(plot_final, "barplot_sc_cluster_per_patient.rds")
plot_final <- readRDS("data/barplot_sc_cluster_per_patient.rds")
  
# Plot cells from all single-cell RNA-seq clusters for each patient
ggplot(
  data=plot_final,
  aes(x=reorder(sample, freq), y= freq, fill = fine_cluster)
  ) +
  geom_bar(stat="identity",
           # position = "fill"
           position = "stack"
           ) +
  facet_grid(~ Case_Control_order, scales = "free", space = "free_x") +
  scale_fill_manual("", values = meta_colors$fine_cluster) +
  labs(
    x = "Donors",
    # y = "Proportion of all clusters"
    y = "Number of cells"
  ) +
  theme_bw(base_size = 25) +
  theme(    
        # axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 25),
        axis.text.x=element_text(angle=35, hjust=1),
        axis.text.y = element_text(size = 25))
ggsave(file = paste("barplot_sc_cluster_per_patient", ".pdf", sep = ""),
       width = 15, height = 12, dpi = 300)
dev.off()


