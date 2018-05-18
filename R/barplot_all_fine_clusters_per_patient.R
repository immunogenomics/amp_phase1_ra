#' ---
#' title: "Plot cells from all single-cell RNA-seq clusters for each patient"
#' author: "Fan Zhang"
#' date: "2018-03-25"
#' Figure 2
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

source("meta_colors.R")
library(gdata) 
library(ggplot2)


# Read post-QC single-cell RNA-seq meta data
meta <- readRDS("../data/celseq_synovium_meta_5452cells_paper.rds")
meta <- meta[order(meta$cell_name),]
meta$cell_name <- as.character(meta$cell_name)
dim(meta)

mat <- meta[, c("sample", "disease", "plate", "fine_cluster", "type")]
non_samples <- names(which(table(mat$sample) == 0))

dat <- table(mat$fine_cluster, mat$sample) 
dat <- as.data.frame(dat)
colnames(dat) <- c("fine_cluster", "sample", "freq")
dat <- dat[-which(dat$sample %in% non_samples),]
dat$sample <- as.character(dat$sample)
table(dat$sample)
dim(dat)

# Read Case.Control labels for OA, non-inflamed RA, and RA
inflam_label <- read.xls("../data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "inflamed RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "non-inflamed RA"

table(inflam_label$Mahalanobis_20)
inflam_label <- inflam_label[which(inflam_label$Patient %in% dat$sample), ]
inflam_label <- inflam_label[, c("Patient", "Disease", "Tissue.type", "lym_25")]
colnames(inflam_label)[1] <- "sample"
plot_final <- merge(dat, inflam_label, by = "sample")
dim(plot_final)
plot_final[1:4,]

# Add cell type column
type <- rep("Fibroblast", nrow(plot_final))
type[grep("B-", plot_final$fine_cluster)] <- "B cell"
type[grep("T-", plot_final$fine_cluster)] <- "T cell"
type[grep("M-", plot_final$fine_cluster)] <- "Monocyte"
plot_final$type <- type

plot_final$lym_25_order = factor(plot_final$lym_25, levels=c('OA','non-inflamed RA','inflamed RA'))
plot_final$type_order = factor(plot_final$type, levels=c('Fibroblast','Monocyte', 'T cell','B cell'))

# saveRDS(plot_final, "barplot_sc_cluster_per_patient.rds")
# plot_final <- readRDS("../data/barplot_sc_cluster_per_patient.rds")
  
# Plot cells from all single-cell RNA-seq clusters for each patient
ggplot(
  data=plot_final,
  aes(x=reorder(sample, freq), y= freq, fill = fine_cluster)
  ) +
  geom_bar(stat="identity",
           # position = "fill"
           position = "stack",
           width = 0.85
           ) +
  facet_grid(type_order ~ lym_25_order, scales = "free", space = "free_x") +
  scale_fill_manual("", values = meta_colors$fine_cluster) +
  labs(
    x = "Donors",
    # y = "Proportion of all clusters"
    y = "Number of cells"
  ) +
  theme_bw(base_size = 27) +
  theme(    
        # axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 25),
        axis.text.x=element_text(angle=30, hjust=1),
        axis.text.y = element_text(size = 30))
ggsave(file = paste("barplot_sc_cluster_per_patient_goodcells", ".pdf", sep = ""),
       width = 14, height = 11, dpi = 300)
dev.off()


# ------
# Plot heatmap of cell count for each cluster and each donor
dat <- table(mat$fine_cluster, mat$sample) 
meta_col <- meta[-which(duplicated(meta$sample)),]
sc_inflam_label <- inflam_label[which(inflam_label$Patient %in% meta_col$sample),]
sc_inflam_label$Patient <- as.character(sc_inflam_label$Patient)
colnames(sc_inflam_label)[1] <- "sample"
# all(meta_col$sample == sc_inflam_label$sample)
meta_col <- merge(meta_col, sc_inflam_label, by = "sample")

meta_row <- meta[-which(duplicated(meta$fine_cluster)),]
annotation_col <- meta_col[, c("Mahalanobis_20", "sample")]
colnames(annotation_col)[1] <- "Case.Control"
annotation_row <- meta_row[, c("type", "fine_cluster")]
rownames(annotation_col) <- as.character(meta_col$sample)
rownames(annotation_row) <- as.character(meta_row$fine_cluster)

dat <- dat[, -which(colnames(dat) %in% names(which(colSums(dat) == 0)))]
# dat <- dat[,order(annotation_col$Case.Control)]
# dat <- dat[c("F-2", "T-5", "B-4", "M-1", "B-1", "B-2", "B-3", "T-2","F-1", "T-6", "T-7",  "F-4", "M-2", "F-3", "M-3", "M-4", "T-1", "T-3", "T-4"),]

pdf("sc_clusters_heatmap_per_donor.pdf", width=10, height=9)
pheatmap(
  mat = dat,
  border_color = NA,
  color = cet_pal(100, name = "kgy"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = meta_colors,
  cellwidth = 12,
  cellheight = 14,
  fontsize = 10,
  fontsize_row = 10,
  scale = 'none',
  legend = TRUE
)
dev.off()



# -------
load("../data/synData.Fibro.downsample.SNE.RData")
load("../data/synData.Bcell.downsample.SNE.RData")
load("../data/synData.Mono.downsample.SNE.RData")
load("../data/synData.Tcell.downsample.SNE.RData")

# Plot # cells from clusters per patient for mass cytometry
fibro_plot <- table(synData.Fibro.downsample$markers, synData.Fibro.downsample$sampleID)
fibro_plot <- as.data.frame(fibro_plot)
colnames(fibro_plot) <- c("cluster","donor", "cells")
fibro_plot$cell_type <- rep("Fibroblast", nrow(fibro_plot))
fibro_plot$cluster <- as.character(fibro_plot$cluster)
# fibro_plot$cluster[which(fibro_plot$cluster == "CD90– CD34– HLA-DR+")] <- "CD90- CD34- HLA-DR+"
# fibro_plot$cluster[which(fibro_plot$cluster == "CD90+ CD34– HLA-DR-")] <- "CD90+ CD34- HLA-DR-"

mono_plot <- table(synData.Mono.downsample$markers, synData.Mono.downsample$sampleID)
mono_plot <- as.data.frame(mono_plot)
colnames(mono_plot) <- c("cluster","donor", "cells")
mono_plot$cell_type <- rep("Monocyte", nrow(mono_plot))
mono_plot$cluster <- as.character(mono_plot$cluster)
mono_plot$cluster[which(mono_plot$cluster == "CD11c+ CD38+ CCR2-")] <- "CD11c+ CD38+"
mono_plot$cluster[which(mono_plot$cluster == "CD11c+ CD38- CCR2- CD64+")] <- "CD11c+ CD38- CD64+"
mono_plot$cluster[which(mono_plot$cluster == "CD11c+ CD38- CCR2-")] <- "CD11c+ CD38-"


tcell_plot <- table(synData.Tcell.downsample$markers, synData.Tcell.downsample$sampleID)
tcell_plot <- as.data.frame(tcell_plot)
colnames(tcell_plot) <- c("cluster","donor", "cells")
tcell_plot$cell_type <- rep("T cell", nrow(tcell_plot))

# synData.Bcell.downsample$merge <- synData.Bcell.downsample$SNE.cluster 
# synData.Bcell.downsample$merge <- as.character(synData.Bcell.downsample$merge)
# synData.Bcell.downsample$merge[synData.Bcell.downsample$SNE.cluster %in% c("6", "8")] <- "mean_6_8"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$SNE.cluster %in% c("1", "3")] <- "mean_1_3"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$SNE.cluster %in% c("10", "12")] <- "mean_10_12"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$SNE.cluster %in% c("16", "17")] <- "mean_16_17"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$SNE.cluster %in% c("14", "15")] <- "mean_14_15"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$SNE.cluster %in% c("4", "9", "7", "13")] <- "mean_4_9_7_13"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$merge == "2"] <- "CM-B2"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$merge == "5"] <- "CM-B5"
# synData.Bcell.downsample$merge[synData.Bcell.downsample$merge == "11"] <- "CM-B11"
# synData.Bcell.downsample <- merge(synData.Bcell.downsample, bcell_merg, by = "merge")
bcell_plot <- table(synData.Bcell.downsample$markers, synData.Bcell.downsample$sampleID)
bcell_plot <- as.data.frame(bcell_plot)
colnames(bcell_plot) <- c("cluster","donor", "cells")
bcell_plot$cell_type <- rep("B cell", nrow(bcell_plot))

cytof_plot <- rbind.data.frame(fibro_plot, mono_plot, tcell_plot, bcell_plot)
cytof_plot$donor <- as.character(cytof_plot$donor)
cytof_plot$donor[which(cytof_plot$donor == "300-0486C")] <- "300-0486"
cytof_plot$donor[which(cytof_plot$donor == "300-0511A")] <- "300-0511"

# Add disease type
inflam_label <- read.xls("../data/postQC_all_samples.xlsx")
inflam_label$Mahalanobis_20 <- rep("OA", nrow(inflam_label))
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis > 20)] <- "inflamed RA"
inflam_label$Mahalanobis_20[which(inflam_label$Mahalanobis < 20 & inflam_label$Disease != "OA")] <- "non-inflamed RA"
table(inflam_label$Mahalanobis_20)
inflam_label <- inflam_label[which(inflam_label$Patient %in% cytof_plot$donor), ]
inflam_label <- inflam_label[, c("Patient", "Disease", "Tissue.type", "Mahalanobis_20")]
colnames(inflam_label)[1] <- "donor"
cytof_plot <- merge(cytof_plot, inflam_label, by = "donor")

cytof_plot$Mahalanobis_20 = factor(cytof_plot$Mahalanobis_20, levels=c('OA','non-inflamed RA','inflamed RA'))
cytof_plot$cell_type = factor(cytof_plot$cell_type, levels=c('Fibroblast', 'Monocyte', 'T cell', 'B cell'))

dm <- cytof_plot[,c("donor", "cluster", "cells")]
dm <- acast(dm, cluster~donor, value.var="cells")
dm[is.na(dm)] <- 0
hc <- hclust(dist(t(dm)))
cytof_plot$donor <- factor(cytof_plot$donor, levels = reorder(hc$labels, cytof_plot$donor))

ggplot(
  data=cytof_plot,
  aes(x=donor, y= cells, fill = cluster)
  ) +
  geom_bar(
           stat="identity",
           # position = "fill"
           position = "stack",
           width = 0.85
  ) +
  facet_grid(cell_type ~ Mahalanobis_20, scales = "free", space = "free_x") +
  scale_fill_manual(values = meta_colors$cytof_cluster) +
  labs(
    x = "Donors",
    y = "Number of cells"
  ) +
  theme_bw(base_size = 15) +
  theme(    
    # axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 15),
    axis.text.x=element_text(angle=35, hjust=1),
    axis.text.y = element_text(size = 15))
ggsave(file = paste("barplot_cytof_cluster_per_patient", ".pdf", sep = ""),
       width = 15, height = 7, dpi = 300)
dev.off()


