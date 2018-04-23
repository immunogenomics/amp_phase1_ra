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

source("pure_functions.R")
source("meta_colors.R")

log2tpm <- readRDS("../data//filtered_log2tpm_lowinput_phase_1.rds")
clin_meta <- readRDS("../data//filtered_meta_lowinput_phase_1.rds")

# Read lym_25 labels for OA, non-inflamed RA, and RA
inflam_label <- read.xls("../data/postQC_all_samples.xlsx")

# # -----------------
# # Plot Figure S1

ggplot() +
  geom_point(
    data = inflam_label,
    mapping = aes(Lymphocytes.Live * 100, CD45pos.ov.Live, fill = lym_25),
    # mapping = aes(Lymphocytes.Live * 100, Monocytes * 100, fill = disease_tissue),
    shape = 21, size = 4, stroke = 0.1
  ) +
  scale_fill_manual(values = meta_colors$tissue_type) +
  geom_vline(xintercept = 25, linetype = "dashed") +
  geom_hline(yintercept = 47, linetype = "dashed") +
  labs(x = "Lymphocytes/Live by flow (%)",
       y = "CD45+/Live by flow (%)"
  ) +
  theme_bw(base_size = 20) +
  theme(
    # axis.ticks = element_blank(),
    # legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(size = 20, color = "black"),
    axis.text.x=element_text(size = 20)
    # axis.text.y = element_text()
  )
ggsave(file = paste("postQC_samplse_mono_lym", ".pdf", sep = ""), width = 7.5, height = 5, dpi = 200)
dev.off()

# inflam_label$lym_25 <- inflam_label$Disease
# inflam_label$lym_25 <- as.character(inflam_label$lym_25)
# inflam_label$lym_25[which(inflam_label$Lymphocytes.Live > 0.25 &
#                            inflam_label$Disease != "OA")] <- "inflamed RA"
# inflam_label$lym_25[which(inflam_label$lym_25 == "RA")] <- "non-inflamed RA"
# table(inflam_label$lym_25)
# 
# inflam_label <- inflam_label[order(inflam_label$disease_tissue),]
# 
# write.table(inflam_label, file = "postQC_all_samples.txt",
#              row.names=T,col.names=T, quote = F, sep = "\t")
# 


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

top_percent <- 0.65
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
cell_type <- "Fibro"
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

# 58 markers
dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-3),])

outFile = "fibro_diff_inflamedRA_OA_lym_25.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("fibro_diff_inflamedRA_OA_lym_25.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                         c("HBEGF", "CLIC5", "HTRA1", "PRG4", "CD55", "DNASE1L3",
                           "PTGFR", "FOS", "F3", "HLA.DRA", # "C3", 
                           "IL6", "HLA.DPA1", "HLA.DRB1", # "IFI30", 
                           "CADM1", "DKK3", # "CAPG", "AKR1C2", "COL8A2",
                           "ACTA2", "CD34", # "MCAM", "MYH11"
                           "IRF1", "CXCL12" # "PDGFRB"
                               )),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- (-1) * dat_plot$logFC
dat_plot$CI.L <- (-1) * dat_plot$CI.L
dat_plot$CI.R <- (-1) * dat_plot$CI.R

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DRB1"),]$gene <- "HLA-DRB1"
dat_plot[which(dat_plot$gene == "HLA.DPA1"),]$gene <- "HLA-DPA1"

# Plot genes based on logFC
# ggplot() +
#   geom_errorbar(data=dat_plot, 
#                 mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
#                 width=0.5, size=0.9, color = "black"
#                 ) +
#   geom_point(data=dat_plot, 
#             aes(x=reorder(gene, logFC),y= logFC),
#             size = 2
#   ) +
#   # scale_color_manual(values = c('grey70', 'black')) +
#   geom_hline(yintercept = 0, size = 0.3) +
#   scale_x_discrete(position = "top") +
#   #scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
#   coord_flip() +
#   labs(x = NULL, 
#        y = "Log2 FC",
#        subtitle = "Inflamed RA versus OA"
#        # title = "abs(Log2 FC) > 0.5 and 5% FDR"
#        ) +
#   theme_bw(base_size = 28) +
#   theme(    
#     # axis.ticks = element_blank(), 
#     panel.grid = element_blank(),
#     axis.text = element_text(size = 30, color = "black"),
#     axis.text.y = element_text(size = 28, face="italic"))
# dev.print("fibro_DE_genes_log2FC.pdf", width = 6, height = 9, dev = pdf)
# dev.off()

dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"


ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, -logFC), ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  # coord_flip() +
  labs(x = NULL, 
       y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 24) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 24, color = "black"),
    axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # axis.text.y = element_text()
  )
dev.print("fibro_DE_genes_log2FC.pdf", width = 11, height = 4, dev = pdf)
dev.off()


# Take the top 20 markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-F1", "C-F2", "C-F3", "C-F4")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=20)# %>% filter(auc > 0.75)
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

# Limma result is the same as t.test
#bulk_m$F1 <- as.numeric(mean_F1[1,])
#t.test(F1 ~ Mahalanobis_20, data = bulk_m) #y = bulk_m$Mahalanobis_20 == "inflamed RA", x = mean_F1[1,])

# single-cell RNA-seq clusters log2(FC) and p-value
d_fibro <- toptable_fdr[, c("logFC", "P.Value")]
d_fibro$cluster <- rownames(d_fibro)
d_fibro$logFC <- -1 * d_fibro$logFC

ggplot(d_fibro, aes(logFC, -log10(P.Value), label = cluster, fill = cluster)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_text_repel(size = 5) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  # coord_cartesian(xlim = c(0, 2.5)) +
  theme(legend.position = "none") +
  labs(x = bquote("Log"[2]~"fold-change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("fibro_fc_pvalue", ".pdf", sep = ""), width = 3.4, height = 3)
dev.off()




# ---------------------
# Mono samples
# ---------------------
cell_type <- "Mono"
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
# 42 markers
dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-3),])


outFile = "mono_diff_inflamedRA_OA_lym_25.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("mono_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                                 c("NR4A2", "ATF3", "PLAUR", "HBEGF", "IFITM3", "CD14", "HLA.DRA", "HLA.DPA1",
                                   "TIMP2", "NUPR1", "ITGB5","HTRA1")),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- (-1) * dat_plot$logFC
dat_plot$CI.L <- (-1) * dat_plot$CI.L
dat_plot$CI.R <- (-1) * dat_plot$CI.R

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DPA1"),]$gene <- "HLA-DPA1"

dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"
dat_plot <- dat_plot[match(c("NR4A2", "ATF3", "PLAUR", "HBEGF", "IFITM3", "CD14", "HLA.DRA", "HLA.DPA1",
                               "TIMP2", "NUPR1", "ITGB5","HTRA1"), rownames(dat_plot)),]
rownames(dat_plot) <- dat_plot$gene

# Plot genes based on logFC
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, -logFC), ymin=CI.L, ymax=CI.R), 
                width=0.5, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  # coord_flip() +
  labs(x = NULL, 
       y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 28) +
  theme(    
    legend.position = "none",
    # axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 26, color = "black"),
    axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # axis.text.y = element_text()
  )
dev.print("mono_DE_genes_log2FC_lym_25.pdf", width = 7, height = 4, dev = pdf)
dev.off()


# Take the top 20 markers for each single-cell cluster and summarize 
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


# single-cell RNA-seq clusters log2(FC) and p-value
d_mono <- toptable_fdr[, c("logFC", "P.Value")]
d_mono$cluster <- rownames(d_mono)
d_mono$logFC <- (-1) * d_mono$logFC

ggplot(d_mono, aes(logFC, -log10(P.Value), label = cluster, fill = cluster)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_text_repel(size = 5) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  labs(x = bquote("Log"[2]~"fold-change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("mono_fc_pvalue", ".pdf", sep = ""), width = 3.2, height = 3)
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

# 0 markers
dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-3),])

outFile = "bcell_diff_inflamedRA_OA_lym_25.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("bcell_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                                 c("FCRL4", "XBP1","IGHG1", "CD38", "SDC1",
                                   "ITGAX","CD19", "IGHM", "MS4A1",
                                   "HLA.DRA", "CXCR5", "HLA.DPB1", "ADIRF",
                                   "RNASE1", "MZB1", "FKBP11" # "DERL3"
                                 )),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- (-1) * dat_plot$logFC
dat_plot$CI.L <- (-1) * dat_plot$CI.L
dat_plot$CI.R <- (-1) * dat_plot$CI.R

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DPB1"),]$gene <- "HLA-DPB1"

dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"

# Plot genes based on logFC
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, -logFC), ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  # scale_color_manual(values = c('grey70', 'black')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  # coord_flip() +
  labs(x = NULL, 
       y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 26) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(size = 24, color = "black"),
    axis.text.x=element_text(angle=45, hjust=1, size = 24)
    # axis.text.y = element_text()
  )
dev.print("bcell_DE_genes_log2FC.pdf", width = 8.5, height = 4, dev = pdf)
dev.off()


# Take the top 20 markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-B1", "C-B2", "C-B3", "C-B4")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=200) %>% filter(auc > 0.75)
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

# 3 markers
dim(toptable_fdr[which(toptable_fdr$logFC < -1 & toptable_fdr$adj.P.Val < 1e-3),])

outFile = "tcell_diff_inflamedRA_OA_lym_25.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("tcell_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                                 c("HLA.DQA2", "CXCL13", "IFNG", # "CCL4L2",
                                   "PDCD1", "NFKBIZ", "CTLA4", # "IL1B",
                                   "TIGIT", "ICOS", "CCR7", # "SELL", # "CD200",
                                   "NFKBID", "NR4A2", "CLECL1","PTGR1" # "SPINT1", "STX3", "SPINT" "GEM", 
                                 )),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- (-1) * dat_plot$logFC
dat_plot$CI.L <- (-1) * dat_plot$CI.L
dat_plot$CI.R <- (-1) * dat_plot$CI.R

# Change HLA.DRA to HLA-DRA
# dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DQA2"),]$gene <- "HLA-DQA2"

dat_plot$up_down <- dat_plot$logFC
dat_plot$up_down[which(dat_plot$up_down > 0)] <- "up"
dat_plot$up_down[which(dat_plot$up_down < 0)] <- "down"

# Plot genes based on logFC
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             aes(x=reorder(gene, logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  # coord_flip() +
  labs(x = NULL, 
       y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 25) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 26, color = "black"),
    axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # axis.text.y = element_text()
  )
dev.print("tcell_DE_genes_log2FC.pdf", width = 9, height = 4, dev = pdf)
dev.off()



# Take the top 20 markers for each single-cell cluster and summarize 
dat_table <- readRDS("../data/cluster_marker_table.rds")
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dat_table <- dat_table[which(dat_table$cluster %in% c("C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7")),]
# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
genes_rp <- grep(pattern = "^RP[SL]", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2, genes_rp)), ]

x <- dat_table %>%
  group_by(cluster) %>% top_n(auc, n=100) %>% filter(auc > 0.75)
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
d_tcell$cluster <- rownames(d_tcell)
d_tcell$logFC <- -1 * d_tcell$logFC

ggplot(d_tcell, aes(logFC, -log10(P.Value), label = cluster, fill = cluster)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values = meta_colors$fine_cluster) +
  geom_text_repel(size = 3) +
  theme_classic(base_size = 17) +
  theme(legend.position = "none") +
  labs(x = bquote("Log"[2]~"fold-change"), y = bquote("-Log"[10]~"P"))
ggsave(file = paste("tcell_fc_pvalue", ".pdf", sep = ""), width = 3.5, height = 3)
dev.off()


