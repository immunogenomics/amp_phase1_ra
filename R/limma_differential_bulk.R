#' ---
#' title: "Differential analysis (inflamed RA vs OA) of each cell type by bulk RNA-seq"
#' author: "Fan Zhang"
#' date: "2018-03-28"
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R")

library(ggplot2)
library(liger)
library(ggrepel)
library(gdata) 
library(limma)
library(ggbeeswarm)

source("pure_functions.R")
source("meta_colors.R")

log2tpm <- readRDS("../data//filtered_log2tpm_lowinput_phase_1.rds")
clin_meta <- readRDS("../data//filtered_meta_lowinput_phase_1.rds")

# Read lym_25 labels for OA, non-inflamed RA, and RA
inflam_label <- read.xls("../data/postQC_all_samples.xlsx")

# -----------------

inflam_label <- inflam_label[, c("Patient", "Mahalanobis_20")]
table(inflam_label$Mahalanobis_20)
colnames(inflam_label)[1] <- "Donor.ID"

clin_merge <- merge(clin_meta, inflam_label, by = "Donor.ID")
table(clin_merge$Mahalanobis_20)

clin_merge <- clin_merge[order(match(clin_merge$Sample.ID, colnames(log2tpm))), ]
all(clin_merge$Sample.ID == colnames(log2tpm))


file_mean_sd <- "celseq_synovium_log2tpm_mean_sd.rds"
if (!file.exists(file_mean_sd)) {
  # dat_mean_sd <- data.frame(
  #   mean  = Matrix::rowMeans(dat),
  #   sd    = apply(dat, 1, sd),
  #   count = apply(dat, 1, function(x) sum(x > 0))
  # )
  dat_mean_sd <- data.frame(
    mean  = Matrix::rowMeans(log2tpm),
    sd    = apply(log2tpm, 1, sd),
    count = apply(log2tpm, 1, function(x) sum(x > 0))
  )
  xs <- dat_mean_sd$mean
  xs[xs == 0] <- 1
  dat_mean_sd$cv <- dat_mean_sd$sd / xs
  rm(xs)
  dat_mean_sd$density_cv <- with(dat_mean_sd, get_density(mean, cv))
  dat_mean_sd$density_sd <- with(dat_mean_sd, get_density(mean, sd))
  saveRDS(dat_mean_sd, file_mean_sd)
} else {
  dat_mean_sd <- readRDS(file_mean_sd)
}

top_percent <- 0.45
mat_idx <- with(
  dat_mean_sd,
  mean > quantile(mean, top_percent) &
    sd > quantile(sd, top_percent)
)

all_samples <- log2tpm[mat_idx,]
dim(all_samples)
all(clin_merge$Sample.ID == colnames(all_samples))

# ---------------------
# Fibroblast samples
# ---------------------
cell_type <- "Fibroblast"
cell_index <- clin_merge$Cell.type==cell_type
bulk_all <- log2tpm[cell_index]
bulk_m <- clin_merge[cell_index,]
bulk_samples <- bulk_all[mat_idx,]
dim(bulk_m)
dim(bulk_samples)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# -------------------------------------------------
# inflamed RA vs OA differential analysis 
bulk_samples <- bulk_samples[, -which(bulk_m$Mahalanobis_20 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$Mahalanobis_20 == "non-inflamed RA"), ]
bulk_m$Mahalanobis_20 <- as.character(bulk_m$Mahalanobis_20)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# # inflamed RA vs non-inflamed differential analysis 
# bulk_m$Mahalanobis_20 <- as.character(bulk_m$Mahalanobis_20)
# bulk_m$Mahalanobis_20[which(bulk_m$Mahalanobis_20 == "non-inflamed RA")] <- "non-inflamed"
# bulk_m$Mahalanobis_20[which(bulk_m$Mahalanobis_20 == "OA")] <- "non-inflamed"
# all(bulk_m$Sample.ID == colnames(bulk_samples))

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ lym_25 + Data.Access.Group + cDNA.plate)
  model.matrix(~ Mahalanobis_20)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
toptable_fdr <-topTable(fit, 
                        coef = "Mahalanobis_20OA",
                        # coef = "Mahalanobis_20non-inflamed",
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")

dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-2),])

outFile = "fibro_diff_inflamedRA_OA_mah.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
gene_order <- c("PTGFR", "FOS", "C3", "CD34", # SC-F1
                "HLA.DRA", "HLA.DPA1", "HLA.DRB1", "IFI30", "IL6", # SC-F2
                "DKK3", "CADM1", "CAPG", "AKR1C2", "COL8A2", # SC-F3
                "HBEGF", "CLIC5", "HTRA4", "ITGA6", "DNASE1L3" #SC-F4
                )
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                         # c("HBEGF", "CLIC5", "HTRA1", "PRG4", "CD55", "DNASE1L3",
                         #  "FOS", "F3", "HLA.DRA", # "C3", "PTGFR",
                         #   "IL6", "HLA.DPA1", "HLA.DRB1", # "IFI30", 
                         #   "CADM1", # "DKK3", # "CAPG", "AKR1C2", "COL8A2",
                         #   "ACTA2", "CD34", # "MCAM", "MYH11"
                         #   "IRF1", "CXCL12" # "PDGFRB")
                         gene_order
                         ),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- dat_plot$logFC * (-1)
# dat_plot$FC <- 2^abs(dat_plot$logFC) * sign(dat_plot$logFC)
dat_plot$CI.L <- dat_plot$CI.L * (-1)
dat_plot$CI.R <- dat_plot$CI.R * (-1)
# dat_plot$CI.L_FC <- (dat_plot$FC + 2^(dat_plot$CI.L - dat_plot$logFC))
# dat_plot$CI.R_FC <- (dat_plot$FC - 2^(dat_plot$logFC - dat_plot$CI.R)) 

# Use the same order with the single-cell genes order
dat_plot <- dat_plot[match(gene_order, dat_plot$gene),]
dat_plot$cluster <- c(rep("SC-F1", 4), rep("SC-F2", 5), rep("SC-F3", 5), rep("SC-F4", 5))

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DRB1"),]$gene <- "HLA-DRB1"
dat_plot[which(dat_plot$gene == "HLA.DPA1"),]$gene <- "HLA-DPA1"

# Label the direction
dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"
# Remove stars
# dat_plot$gene[which(dat_plot$P.Value < 1e-4)] <- paste(dat_plot$gene[which(dat_plot$P.Value < 1e-4)], "****", sep = "")
# dat_plot$gene[which(dat_plot$P.Value < 1e-3 & dat_plot$P.Value > 1e-4)] <- paste(dat_plot$gene[which(dat_plot$P.Value < 1e-3 & dat_plot$P.Value > 1e-4)], "***", sep = "")
# dat_plot$gene[which(dat_plot$P.Value < 1e-2 & dat_plot$P.Value > 1e-3)] <- paste(dat_plot$gene[which(dat_plot$P.Value < 1e-2 & dat_plot$P.Value > 1e-3)], "**", sep = "")

# ---
# Plot sc cluster markers in bulk
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
                # mapping=aes(x=gene, ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             mapping=aes(x=reorder(gene, logFC), y= logFC),
             # mapping=aes(x=gene, y= logFC, color = up_down),
             size = 3
  ) +
  facet_grid(cluster ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  coord_flip() +
  labs(x = NULL, 
       y = NULL
       # y = "Fold change"
  ) +
  theme_bw(base_size = 20) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 20, color = "black"),
    axis.text.y = element_text(hjust=1, size=20, face = "italic"),
    # axis.text.x=element_text(angle=45, hjust=1, size = 25)
    strip.background = element_rect(fill=meta_colors$fine_cluster)
  )
# dev.print("fibro_DE_genes_FC.pdf", width = 9.5, height = 5, dev = pdf)
dev.print("fibro_sc_markers_in_bulk.pdf", width = 5, height = 6, dev = pdf)
dev.off()



# ------------------------------------------------------------
# Only plot one gene
bulk_m$gene <- as.numeric(bulk_samples[which(rownames(bulk_samples) == "HLA.DRA"),])
dat_median <- bulk_m %>% group_by(Mahalanobis_20) %>% summarise(median = median(gene))
bulk_m$Mahalanobis_20 <- as.character(bulk_m$Mahalanobis_20)
bulk_m$Mahalanobis_20[which(bulk_m$Mahalanobis_20 == "inflamed RA")] <- "leukocyte-rich RA"
bulk_m$Mahalanobis_20[which(bulk_m$Mahalanobis_20 == "non-inflamed RA")] <- "leukocyte-poor RA"
bulk_m$Mahalanobis_20 <- factor(bulk_m$Mahalanobis_20,
                                levels = c('OA','leukocyte-poor RA', "leukocyte-rich RA"),ordered = TRUE)

bulk_m <- bulk_m[-which(bulk_m$gene > 10 & bulk_m$Mahalanobis_20 == "leukocyte-poor RA"),]

ggplot(data=bulk_m, 
       mapping = aes(x=Mahalanobis_20, y=gene, fill=Mahalanobis_20)
  ) +
  geom_quasirandom(
    shape = 21, size = 4.5, stroke = 0.35
  ) +
  stat_summary(
    fun.y = median, fun.ymin = median, fun.ymax = median,
    geom = "crossbar", width = 0.8
  ) +
  scale_fill_manual(values = meta_colors$Case.Control) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  xlab('')+ylab('Gene expression')+
  theme_bw(base_size = 25) +
  labs(title = "HLA-DRA") +
  theme(
    legend.position = "none",
    axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 25, color = "black"),
    axis.text.x = element_text(angle=35, hjust=1, size=25),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(color="black", size=25)) +
  theme(legend.text=element_text(size=25)) +
  coord_cartesian(ylim = c(2.2, 12))
  ggsave(
  file = "bulk_HLA-DRA_fibro_3disease.pdf",
  width = 4.5, height = 7
)

t.test(bulk_m$gene[which(bulk_m$Mahalanobis_20 == "leukocyte-rich RA")],
       bulk_m$gene[which(bulk_m$Mahalanobis_20 == "leukocyte-poor RA")],
       alternative ="greater")

# ---







# Take the top markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-F1", "C-F2", "C-F3", "C-F4")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=100) %>% filter(auc > 0.7)
x$cluster <- factor(x$cluster)
table(x$cluster)
sc_markers <- as.character(unique(x$gene))
sc_markers <- gsub("[-]", ".", sc_markers)

bulk_sum <- bulk_all[,which(colnames(bulk_all) %in% bulk_m$Sample.ID)]
bulk_sum <- bulk_sum[which(rownames(bulk_sum) %in% sc_markers),]
dim(bulk_sum)

mean_F1 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-F1"),]$gene),])))
mean_F2 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-F2"),]$gene),])))
mean_F3 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-F3"),]$gene),])))
mean_F4 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-F4"),]$gene),])))
mean_sc <- rbind.data.frame(mean_F1, mean_F2, mean_F3, mean_F4)
rownames(mean_sc) <- c("C-F1", "C-F2", "C-F3", "C-F4")
all(colnames(mean_sc) == bulk_m$Sample.ID)

des <- with(
  bulk_m,
  model.matrix(~ Mahalanobis_20)
)
fit_2 <- lmFit(object = mean_sc, design = des)
fit_2 <- eBayes(fit_2)
toptable_fdr <- topTable(fit_2, coef = "Mahalanobis_20OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
# single-cell RNA-seq clusters log2(FC) and p-value
toptable_fdr[, c("logFC", "P.Value")]

# ---
# Limma result is the same as t.test
bulk_m$F4 <- as.numeric(mean_F4[1,])
t.test(F4 ~ Mahalanobis_20, data = bulk_m) 
#y = bulk_m$Mahalanobis_20 == "inflamed RA", x = mean_F1[1,])

# ---
# single-cell RNA-seq clusters log2(FC) and p-value
d_fibro <- toptable_fdr[, c("logFC", "CI.L", "CI.R", "t", "P.Value")]
d_fibro$cluster <- rownames(d_fibro)
d_fibro$cluster <- paste("S", d_fibro$cluster, sep="")
d_fibro$logFC <- d_fibro$logFC * (-1)
# d_fibro$FC <- 2^abs(d_fibro$logFC) * sign(d_fibro$logFC)
d_fibro$CI.L <- d_fibro$CI.L * (-1)
d_fibro$CI.R <- d_fibro$CI.R * (-1)
# d_fibro$CI.L_FC <- (d_fibro$FC + 2^(d_fibro$CI.L - d_fibro$logFC))
# d_fibro$CI.R_FC <- (d_fibro$FC - 2^(d_fibro$logFC - d_fibro$CI.R)) 


ggplot(
  d_fibro, 
  aes(logFC, -log10(P.Value), label = cluster, fill = cluster)
  ) +
  geom_point(size = 4.5,  shape = 21) +
  geom_errorbarh(aes(xmin=CI.L, xmax=CI.R, y = -log10(P.Value)),
                  size=0.25, color = "black"
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(labels = function(x) round(2^abs(x), 1)) +
  # scale_x_continuous(breaks=c(-4,-2,0,2,4)) +
  # geom_rect(aes(xmin = 0.2, xmax = 5, ymin = 1.4, ymax = 7),
  #           fill = "orange", alpha = 0.1) +
  # geom_rect(aes(xmin = -5, xmax = -0.2, ymin = 4.5, ymax = 7),
  #           fill = "purple", alpha = 0.1) +
  geom_text_repel(
    data = d_fibro,
    aes(x = logFC, y = -log10(P.Value), label = cluster),
    size = 5, color = "black",
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.4, "lines")
  ) +
  theme_classic(base_size = 16) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x=element_text(size = 16, color = "black")
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme(legend.position = "none") +
  labs(x = ("Fold change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("fibro_fc_pvalue", ".pdf", sep = ""), width = 3.2, height = 3)
dev.off()




# ---------------------
# Mono samples
# ---------------------
cell_type <- "Monocyte"
cell_index <- clin_merge$Cell.type==cell_type
bulk_all <- log2tpm[cell_index]
bulk_m <- clin_merge[cell_index,]
bulk_samples <- bulk_all[mat_idx,]
dim(bulk_m)
dim(bulk_samples)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# -------------------------------------------------
# inflamed RA vs OA differential analysis 
bulk_samples <- bulk_samples[, -which(bulk_m$Mahalanobis_20 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$Mahalanobis_20 == "non-inflamed RA"), ]
bulk_m$Mahalanobis_20 <- as.character(bulk_m$Mahalanobis_20)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ lym_25 + Data.Access.Group + cDNA.plate)
  model.matrix(~ Mahalanobis_20)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
toptable_fdr <-topTable(fit, coef = "Mahalanobis_20OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[order(toptable_fdr$logFC, decreasing = T),][1:30,]

# 
dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-2),])


# outFile = "mono_diff_inflamedRA_OA_mah.txt"
# write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
gene_order <- c("NR4A2", "ATF3", "PLAUR", "HBEGF", "IL1B", "RGS2", # SC-M1
                "VSIG4", "NUPR1", "GPNMB","HTRA1", "MERTK", "CTSK", # SC-M2
                "C1QA", "MARCO", "CD14", # SC-M3
                "IFITM3", "SPP1", "IFI6", "LY6E" # SC-M4
)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% gene_order),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- dat_plot$logFC * (-1)
# dat_plot$FC <- 2^abs(dat_plot$logFC) * sign(dat_plot$logFC)
dat_plot$CI.L <- dat_plot$CI.L * (-1)
dat_plot$CI.R <- dat_plot$CI.R * (-1)
# dat_plot$CI.L_FC <- (dat_plot$FC + 2^(dat_plot$CI.L - dat_plot$logFC))
# dat_plot$CI.R_FC <- (dat_plot$FC - 2^(dat_plot$logFC - dat_plot$CI.R)) 

# Use the same order with the single-cell genes order
dat_plot <- dat_plot[match(gene_order, dat_plot$gene),]
dat_plot$cluster <- c(rep("SC-M1", 6), rep("SC-M2", 6), rep("SC-M3", 3), rep("SC-M4", 4))

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DPA1"),]$gene <- "HLA-DPA1"
dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"
# dat_plot$gene[which(dat_plot$P.Value < 1e-4)] <- paste(dat_plot$gene[which(dat_plot$P.Value < 1e-4)], "****", sep = "")
# dat_plot$gene[which(dat_plot$P.Value < 1e-3 & dat_plot$P.Value > 1e-4)] <- paste(dat_plot$gene[which(dat_plot$P.Value < 1e-3 & dat_plot$P.Value > 1e-4)], "***", sep = "")
# dat_plot$gene[which(dat_plot$P.Value < 1e-2 & dat_plot$P.Value > 1e-3)] <- paste(dat_plot$gene[which(dat_plot$P.Value < 1e-2 & dat_plot$P.Value > 1e-3)], "**", sep = "")

# Plot genes based on logFC
# Plot sc cluster markers in bulk
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
                # mapping=aes(x=gene, ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             # mapping=aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             mapping=aes(x=reorder(gene, logFC), y= logFC),
             # mapping=aes(x=gene, y= logFC, color = up_down),
             size = 3
  ) +
  facet_grid(cluster ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  coord_flip() +
  labs(x = NULL, 
       y = NULL
       # y = "Fold change"
  ) +
  theme_bw(base_size = 21) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 21, color = "black"),
    axis.text.y = element_text(hjust=1, size=21, face = "italic")
    # axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # strip.background = element_rect(fill=meta_colors$fine_cluster)
  )
# dev.print("fibro_DE_genes_FC.pdf", width = 9.5, height = 5, dev = pdf)
dev.print("mono_sc_markers_in_bulk.pdf", width = 5, height = 7.5, dev = pdf)
dev.off()




# Take the top markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-M1", "C-M2", "C-M3", "C-M4")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=100) %>% filter(auc > 0.7)
x$cluster <- factor(x$cluster)
table(x$cluster)
sc_markers <- as.character(unique(x$gene))
sc_markers <- gsub("[-]", ".", sc_markers)

bulk_sum <- bulk_all[,which(colnames(bulk_all) %in% bulk_m$Sample.ID)]
bulk_sum <- bulk_sum[which(rownames(bulk_sum) %in% sc_markers),]
dim(bulk_sum)

mean_M1 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-M1"),]$gene),])))
mean_M2 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-M2"),]$gene),])))
mean_M3 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-M3"),]$gene),])))
mean_M4 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-M4"),]$gene),])))
mean_sc <- rbind.data.frame(mean_M1, mean_M2, mean_M3, mean_M4)
rownames(mean_sc) <- c("C-M1", "C-M2", "C-M3", "C-M4")
all(colnames(mean_sc) == bulk_m$Sample.ID)

des <- with(
  bulk_m,
  model.matrix(~ Mahalanobis_20)
)
fit_2 <- lmFit(object = mean_sc, design = des)
fit_2 <- eBayes(fit_2)
toptable_fdr <- topTable(fit_2, coef = "Mahalanobis_20OA", 
                         number=nrow(bulk_samples), adjust.method="BH", 
                         confint=TRUE,
                         sort.by="p")
# single-cell RNA-seq clusters log2(FC) and p-value
toptable_fdr[, c("logFC", "P.Value")]

# ---
# Limma result is the same as t.test
bulk_m$M2 <- as.numeric(mean_M2[1,])
t.test(M2 ~ Mahalanobis_20, data = bulk_m, alternative = "two.sided") 


# ---
# single-cell RNA-seq clusters log2(FC) and p-value
d_mono <- toptable_fdr[, c("logFC", "CI.L", "CI.R", "t", "P.Value")]
d_mono$cluster <- rownames(d_mono)
d_mono$cluster <- paste("S", d_mono$cluster, sep="")
d_mono$logFC <- d_mono$logFC * (-1)
# d_mono$FC <- 2^abs(d_mono$logFC) * sign(d_mono$logFC)
d_mono$CI.L <- d_mono$CI.L * (-1)
d_mono$CI.R <- d_mono$CI.R * (-1)
# d_mono$CI.L_FC <- (d_mono$FC + 2^(d_mono$CI.L - d_mono$logFC))
# d_mono$CI.R_FC <- (d_mono$FC - 2^(d_mono$logFC - d_mono$CI.R)) 


ggplot(
  d_mono, 
  aes(logFC, -log10(P.Value), label = cluster, fill = cluster)
  ) +
  # geom_rect(aes(xmin = 0.1, xmax = 4, ymin = 1.4, ymax = 5.5),
  #           fill = "orange", alpha = 0.1) +
  # geom_rect(aes(xmin = -4, xmax = -0.1, ymin = 3.5, ymax = 5.5),
  #           fill = "purple", alpha = 0.1) +
  geom_point(size = 4.5,  shape = 21) +
  geom_errorbarh(aes(xmin=CI.L, xmax=CI.R, y = -log10(P.Value)),
                 size=0.25, color = "black"
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(labels = function(x) round(2^abs(x), 1)) +
  # scale_x_continuous(breaks=c(-4,-2,0,2,4)) +
  geom_text_repel(
    data = d_mono,
    aes(x = logFC, y = -log10(P.Value), label = cluster),
    size = 5, color = "black",
    box.padding = unit(0.7, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_classic(base_size = 16) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x=element_text(size = 16, color = "black")
  ) +
  # coord_cartesian(xlim = c(-2,2), ylim = c(0,5)) +
  scale_y_continuous(breaks=c(0,1,3,5)) +
  theme(legend.position = "none") +
  labs(x = ("Fold change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("mono_fc_pvalue", ".pdf", sep = ""), width = 3.3, height = 3.1)
dev.off()



# ---------------------
# B cell samples
# ---------------------
cell_type <- "B cell"
cell_index <- clin_merge$Cell.type==cell_type
bulk_all <- log2tpm[cell_index]
bulk_m <- clin_merge[cell_index,]
bulk_samples <- bulk_all[mat_idx,]
dim(bulk_m)
dim(bulk_samples)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# -------------------------------------------------
# inflamed RA vs OA differential analysis 
bulk_samples <- bulk_samples[, -which(bulk_m$Mahalanobis_20 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$Mahalanobis_20 == "non-inflamed RA"), ]
bulk_m$Mahalanobis_20 <- as.character(bulk_m$Mahalanobis_20)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# # inflamed RA vs (OA + non-inflamed RA) differential analysis 
# bulk_m[which(bulk_m$lym_25 == "non-inflamed RA"), ]$lym_25 <- "OA"
# bulk_m$lym_25 <- as.character(bulk_m$lym_25)
# all(bulk_m$Sample.ID == colnames(bulk_samples))
# table(bulk_m$lym_25)

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ lym_25 + Data.Access.Group + cDNA.plate)
  model.matrix(~ Mahalanobis_20)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
colnames(fit$coefficients)
toptable_fdr <-topTable(fit, coef = "Mahalanobis_20OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[order(toptable_fdr$logFC, decreasing = F),][1:50,]

# 
dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-1),])

# outFile = "bcell_diff_inflamedRA_OA_mah.txt"
# write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
gene_order <- c("CD83", "CD38", "IL6", "IGHD", "IGHM", "BACH2", # SC-B1
                "HLA.DPB1", "HLA.DRA", "MS4A1", # SC-B2
                "ITGAX", "TBX21", "ZEB2", "ACTB", # SC-B3
                "XBP1", "MZB1", "FKBP11", "SSR4", "DERL3" # SC-B4
)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% gene_order),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- dat_plot$logFC * (-1)
# dat_plot$FC <- 2^abs(dat_plot$logFC) * sign(dat_plot$logFC)
dat_plot$CI.L <- dat_plot$CI.L * (-1)
dat_plot$CI.R <- dat_plot$CI.R * (-1)

# Use the same order with the single-cell genes order
dat_plot <- dat_plot[match(gene_order, dat_plot$gene),]
dat_plot$cluster <- c(rep("SC-B1", 6), rep("SC-B2", 3), rep("SC-B3", 4), rep("SC-B4", 5))

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DPB1"),]$gene <- "HLA-DPB1"

dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"


# Plot genes based on logFC
# Plot sc cluster markers in bulk
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
                # mapping=aes(x=gene, ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             # mapping=aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             mapping=aes(x=reorder(gene, logFC), y= logFC),
             # mapping=aes(x=gene, y= logFC, color = up_down),
             size = 3
  ) +
  facet_grid(cluster ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  coord_flip() +
  labs(x = NULL, 
       y = NULL
       # y = "Fold change"
  ) +
  theme_bw(base_size = 21) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 21, color = "black"),
    axis.text.y = element_text(hjust=1, size=21, face = "italic")
    # axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # strip.background = element_rect(fill=meta_colors$fine_cluster)
  )
# dev.print("fibro_DE_genes_FC.pdf", width = 9.5, height = 5, dev = pdf)
dev.print("bcell_sc_markers_in_bulk.pdf", width = 5, height = 7.5, dev = pdf)
dev.off()




# Take the top markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-B1", "C-B2", "C-B3", "C-B4")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=100) %>% filter(auc > 0.7)
x$cluster <- factor(x$cluster)
table(x$cluster)
sc_markers <- as.character(unique(x$gene))
sc_markers <- gsub("[-]", ".", sc_markers)

bulk_sum <- bulk_all[,which(colnames(bulk_all) %in% bulk_m$Sample.ID)]
bulk_sum <- bulk_sum[which(rownames(bulk_sum) %in% sc_markers),]
dim(bulk_sum)

mean_B1 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-B1"),]$gene),])))
mean_B2 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-B2"),]$gene),])))
mean_B3 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-B3"),]$gene),])))
mean_B4 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-B4"),]$gene),])))
mean_sc <- rbind.data.frame(mean_B1, mean_B2, mean_B3, mean_B4)
rownames(mean_sc) <- c("C-B1", "C-B2", "C-B3", "C-B4")
all(colnames(mean_sc) == bulk_m$Sample.ID)

des <- with(
  bulk_m,
  model.matrix(~ Mahalanobis_20)
)
fit_2 <- lmFit(object = mean_sc, design = des)
fit_2 <- eBayes(fit_2)
toptable_fdr <- topTable(fit_2, coef = "Mahalanobis_20OA", 
                         number=nrow(bulk_samples), adjust.method="BH", 
                         confint=TRUE,
                         sort.by="p")
# single-cell RNA-seq clusters log2(FC) and p-value
toptable_fdr[, c("logFC", "P.Value")]


# single-cell RNA-seq clusters log2(FC) and p-value
d_bcell <- toptable_fdr[, c("logFC", "P.Value")]
d_bcell$cluster <- rownames(d_bcell)
d_bcell$logFC <- (-1) * d_bcell$logFC

ggplot(d_bcell, aes(logFC, -log10(P.Value), label = cluster, fill = cluster)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_text_repel(size = 5) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = bquote("Log"[2]~"fold-change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("bcell_fc_pvalue", ".pdf", sep = ""), width = 3.5, height = 3)
dev.off()



# ---------------------
# T cell samples
# ---------------------
cell_type <- "T cell"
cell_index <- clin_merge$Cell.type==cell_type
bulk_all <- log2tpm[cell_index]
bulk_m <- clin_merge[cell_index,]
bulk_samples <- bulk_all[mat_idx,]
dim(bulk_m)
dim(bulk_samples)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# -------------------------------------------------
# inflamed RA vs OA differential analysis 
bulk_samples <- bulk_samples[, -which(bulk_m$Mahalanobis_20 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$Mahalanobis_20 == "non-inflamed RA"), ]
bulk_m$Mahalanobis_20 <- as.character(bulk_m$Mahalanobis_20)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# # inflamed RA vs (OA + non-inflamed RA) differential analysis 
# bulk_m[which(bulk_m$lym_25 == "non-inflamed RA"), ]$lym_25 <- "OA"
# bulk_m$lym_25 <- as.character(bulk_m$lym_25)
# all(bulk_m$Sample.ID == colnames(bulk_samples))
# table(bulk_m$lym_25)

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ lym_25 + Data.Access.Group + cDNA.plate)
  model.matrix(~ Mahalanobis_20)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
colnames(fit$coefficients)
toptable_fdr <-topTable(fit, coef = "Mahalanobis_20OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")

dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-2),])

# outFile = "tcell_diff_inflamedRA_OA_mah.txt"
# write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
gene_order <- c("IL7R", "SELL", "CCR7", "NFKBIZ", "LEF1", # SC-T1
                "FOXP3", "CTLA4", "TIGIT", "DUSP4",  # SC-T2
                "CXCL13", "PDCD1", "CD200", # SC-T3
                "GZMK", "CD8A", "NKG7", "GZMA", # SC-T4
                "GNLY", "GZMB", "ZNF683", "PRF1", # SC-T5
                "HLA.DRB1", "GZMB", "IFNG" # SC-T6
)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% gene_order),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- dat_plot$logFC * (-1)
# dat_plot$FC <- 2^abs(dat_plot$logFC) * sign(dat_plot$logFC)
dat_plot$CI.L <- dat_plot$CI.L * (-1)
dat_plot$CI.R <- dat_plot$CI.R * (-1)

# Use the same order with the single-cell genes order
dat_plot <- dat_plot[match(gene_order, dat_plot$gene),]
dat_plot$cluster <- c(rep("SC-T1", 5), rep("SC-T2", 4), rep("SC-T3", 3), rep("SC-T4", 4),
                      rep("SC-T5", 4), rep("SC-T6", 3))

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRB5"),]$gene <- "HLA-DRB5"
dat_plot[which(dat_plot$gene == "HLA.DQA2"),]$gene <- "HLA-DQA2"
dat_plot[which(dat_plot$gene == "HLA.DRB1"),]$gene <- "HLA-DRB1"
dat_plot[which(dat_plot$gene == "HLA.DQA1"),]$gene <- "HLA-DQA1"

dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"


# Plot genes based on logFC
# Plot sc cluster markers in bulk
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
                # mapping=aes(x=gene, ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             # mapping=aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             mapping=aes(x=reorder(gene, logFC), y= logFC),
             # mapping=aes(x=gene, y= logFC, color = up_down),
             size = 3
  ) +
  facet_grid(cluster ~ ., scales = "free", space = "free_x") +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  coord_flip() +
  labs(x = NULL, 
       y = NULL
       # y = "Fold change"
  ) +
  theme_bw(base_size = 21) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 21, color = "black"),
    axis.text.y = element_text(hjust=1, size=21, face = "italic")
    # axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # strip.background = element_rect(fill=meta_colors$fine_cluster)
  )
# dev.print("fibro_DE_genes_FC.pdf", width = 9.5, height = 5, dev = pdf)
dev.print("Tcell_sc_markers_in_bulk.pdf", width = 5, height = 8.5, dev = pdf)
dev.off()



# Take the top markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=100) %>% filter(auc > 0.7)
x$cluster <- factor(x$cluster)
table(x$cluster)
sc_markers <- as.character(unique(x$gene))
sc_markers <- gsub("[-]", ".", sc_markers)

bulk_sum <- bulk_all[,which(colnames(bulk_all) %in% bulk_m$Sample.ID)]
bulk_sum <- bulk_sum[which(rownames(bulk_sum) %in% sc_markers),]
dim(bulk_sum)

mean_T1 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T1"),]$gene),])))
mean_T2 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T2"),]$gene),])))
mean_T3 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T3"),]$gene),])))
mean_T4 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T4"),]$gene),])))
mean_T5 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T5"),]$gene),])))
mean_T6 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T6"),]$gene),])))
mean_T7 <- t(as.data.frame(colMeans(bulk_sum[which(rownames(bulk_sum) %in% x[which(x$cluster == "C-T7"),]$gene),])))
mean_sc <- rbind.data.frame(mean_T1, mean_T2, mean_T3, mean_T4, mean_T5, mean_T6, mean_T7)
rownames(mean_sc) <- c("C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7")
all(colnames(mean_sc) == bulk_m$Sample.ID)

des <- with(
  bulk_m,
  model.matrix(~ Mahalanobis_20)
)
fit_2 <- lmFit(object = mean_sc, design = des)
fit_2 <- eBayes(fit_2)
toptable_fdr <- topTable(fit_2, coef = "Mahalanobis_20OA", 
                         number=nrow(bulk_samples), adjust.method="BH", 
                         confint=TRUE,
                         sort.by="p")

# single-cell RNA-seq clusters log2(FC) and p-value
d_tcell <- toptable_fdr[, c("logFC", "P.Value")]

# ---
# Limma result is the same as t.test
bulk_m$T4 <- as.numeric(mean_T4[1,])
t.test(T4 ~ Mahalanobis_20, data = bulk_m, alternative = "two.sided") 


# ---
# single-cell RNA-seq clusters log2(FC) and p-value
d_tcell <- toptable_fdr[, c("logFC", "CI.L", "CI.R", "t", "P.Value")]
d_tcell$cluster <- rownames(d_tcell)
d_tcell$cluster <- paste("S", d_tcell$cluster, sep="")
d_tcell$logFC <- d_tcell$logFC * (-1)
# d_tcell$FC <- 2^abs(d_tcell$logFC) * sign(d_tcell$logFC)
d_tcell$CI.L <- d_tcell$CI.L * (-1)
d_tcell$CI.R <- d_tcell$CI.R * (-1)
# d_tcell$CI.L_FC <- (d_tcell$FC + 2^(d_tcell$CI.L - d_tcell$logFC))
# d_tcell$CI.R_FC <- (d_tcell$FC - 2^(d_tcell$logFC - d_tcell$CI.R)) 


ggplot(
  d_tcell, 
  aes(logFC, -log10(P.Value), label = cluster, fill = cluster)
  ) +
  # geom_rect(aes(xmin = 0.1, xmax = 4, ymin = 1.5, ymax = 2.5),
            # fill = "orange", alpha = 0.1) +
  geom_point(size = 4.5,  shape = 21) +
  geom_errorbarh(aes(xmin=CI.L, xmax=CI.R, y = -log10(P.Value)),
                 size=0.2, color = "black"
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(labels = function(x) round(2^abs(x), 1)) +
  geom_text_repel(
    data = d_tcell,
    aes(x = logFC, y = -log10(P.Value), label = cluster),
    size = 5, color = "black",
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  theme_classic(base_size = 16) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x=element_text(size = 16, color = "black")
  ) +
  coord_cartesian(ylim = c(0,2)) +
  scale_y_continuous(breaks=c(0,1,2)) +
  theme(legend.position = "none") +
  labs(x = ("Fold change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("tcell_fc_pvalue", ".pdf", sep = ""), width = 3.5, height = 3.5)
dev.off()



