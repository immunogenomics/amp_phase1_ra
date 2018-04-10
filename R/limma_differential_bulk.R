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


# Read CD45_49 labels for OA, non-inflamed RA, and RA
inflam_label <- read.xls("../data-raw/postQC_all_samples.xlsx")

# # -----------------
# # Plot Figure S1

ggplot() +
  geom_point(
    data = inflam_label,
    mapping = aes(Lymphocytes.Live * 100, CD45pos.ov.Live, fill = CD45_49),
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

inflam_label <- inflam_label[, c("Patient", "CD45_49")]
table(inflam_label$CD45_49)
colnames(inflam_label)[1] <- "Donor.ID"

clin_merge <- merge(clin_meta, inflam_label, by = "Donor.ID")
table(clin_merge$CD45_49)

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
bulk_samples <- bulk_samples[, -which(bulk_m$CD45_49 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$CD45_49 == "non-inflamed RA"), ]
bulk_m$CD45_49 <- as.character(bulk_m$CD45_49)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ CD45_49 + Data.Access.Group + cDNA.plate)
  model.matrix(~ CD45_49)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
toptable_fdr <-topTable(fit, coef = "CD45_49OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[1:4,]

outFile = "fibro_diff_inflamedRA_OA_CD45_49.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("fibro_diff_inflamedRA_OA_CD45_49.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                         c("HBEGF", "CLIC5", "HTRA1", "PRG4", "CD55", "DNASE1L3",
                           "PTGFR", "FOS", "F3", "HLA.DRA", # "C3", 
                           "IL6", "HLA.DPA1", "HLA.DRB1", # "IFI30", 
                           "CADM1", # "DKK3", "CAPG", "AKR1C2", "COL8A2",
                           "ACTA2", "CD34", # "MCAM", "MYH11"
                           "IRF1", "CXCL12", "PDGFRB"
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
                width=0.5, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             aes(x=reorder(gene, -logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
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
dev.print("fibro_DE_genes_log2FC_CD45_49.pdf", width = 11.5, height = 4.5, dev = pdf)
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
bulk_samples <- bulk_samples[, -which(bulk_m$CD45_49 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$CD45_49 == "non-inflamed RA"), ]
bulk_m$CD45_49 <- as.character(bulk_m$CD45_49)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ CD45_49 + Data.Access.Group + cDNA.plate)
  model.matrix(~ CD45_49)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
toptable_fdr <-topTable(fit, coef = "CD45_49OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[1:4,]
toptable_fdr[order(toptable_fdr$logFC, decreasing = T),][1:30,]

outFile = "mono_diff_inflamedRA_OA_CD45_49.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("mono_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                                 c("PLAUR", "CD38", "IFITM3", "HBEGF", "HLA.DRA",
                                   "HLA.DPA1", "ATF3",
                                   "TIMP2", "ITGB5", "NUPR1", "HTRA1"
                                 )),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- (-1) * dat_plot$logFC
dat_plot$CI.L <- (-1) * dat_plot$CI.L
dat_plot$CI.R <- (-1) * dat_plot$CI.R

# Change HLA.DRA to HLA-DRA
dat_plot[which(dat_plot$gene == "HLA.DRA"),]$gene <- "HLA-DRA"
dat_plot[which(dat_plot$gene == "HLA.DQA1"),]$gene <- "HLA-DQA1"
dat_plot[which(dat_plot$gene == "HLA.DPA1"),]$gene <- "HLA-DPA1"

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
             aes(x=reorder(gene, -logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
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
dev.print("mono_DE_genes_log2FC_CD45_49.pdf", width = 9, height = 4, dev = pdf)
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
bulk_samples <- bulk_samples[, -which(bulk_m$CD45_49 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$CD45_49 == "non-inflamed RA"), ]
bulk_m$CD45_49 <- as.character(bulk_m$CD45_49)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# # inflamed RA vs (OA + non-inflamed RA) differential analysis 
# bulk_m[which(bulk_m$CD45_49 == "non-inflamed RA"), ]$CD45_49 <- "OA"
# bulk_m$CD45_49 <- as.character(bulk_m$CD45_49)
# all(bulk_m$Sample.ID == colnames(bulk_samples))
# table(bulk_m$CD45_49)

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ CD45_49 + Data.Access.Group + cDNA.plate)
  model.matrix(~ CD45_49)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
colnames(fit$coefficients)
toptable_fdr <-topTable(fit, coef = "CD45_49OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[1:4,]
toptable_fdr[order(toptable_fdr$logFC, decreasing = F),][1:50,]

outFile = "bcell_diff_inflamedRA_OA_CD45_49.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("bcell_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                                 c("FCRL4", "XBP1","IGHG1", "CD38", "SDC1",
                                   "ITGAX","CD19", "IGHA1", "IGHM", "MS4A1",
                                   "HLA.DRA", "CXCR5", "HLA.DPB1", "ADIRF",
                                   "RNASE1"
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
             aes(x=reorder(gene, -logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  # scale_color_manual(values = c('grey70', 'black')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  labs(x = NULL, 
       y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 28) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_text(size = 26, color = "black"),
    axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # axis.text.y = element_text()
  )
dev.print("bcell_DE_genes_log2FC_CD45_49.pdf", width = 10, height = 3.7, dev = pdf)
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
bulk_samples <- bulk_samples[, -which(bulk_m$CD45_49 == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$CD45_49 == "non-inflamed RA"), ]
bulk_m$CD45_49 <- as.character(bulk_m$CD45_49)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# # inflamed RA vs (OA + non-inflamed RA) differential analysis 
# bulk_m[which(bulk_m$CD45_49 == "non-inflamed RA"), ]$CD45_49 <- "OA"
# bulk_m$CD45_49 <- as.character(bulk_m$CD45_49)
# all(bulk_m$Sample.ID == colnames(bulk_samples))
# table(bulk_m$CD45_49)

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ CD45_49 + Data.Access.Group + cDNA.plate)
  model.matrix(~ CD45_49)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
colnames(fit$coefficients)
toptable_fdr <-topTable(fit, coef = "CD45_49OA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[order(toptable_fdr$logFC, decreasing = T),][1:50,]

outFile = "tcell_diff_inflamedRA_OA.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("tcell_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                                 c("HLA.DQA2", "CXCL13", "CCL4L2", "IFNG",
                                   "PDCD1", "NFKBIZ", "CTLA4", # "IL1B",
                                   "TIGIT", "ICOS", "CCR7", "SELL", # "CD200",
                                   "GEM", "NFKBID", "NR4A2", "CLECL1","PTGR1" # "SPINT1", "STX3", "SPINT"
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
                mapping=aes(x=reorder(gene, -logFC), ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.9, color = "black"
  ) +
  geom_point(data=dat_plot, 
             aes(x=reorder(gene, -logFC), y= logFC, color = up_down),
             size = 3
  ) +
  scale_color_manual(values = c('#6A3D9A', '#FF7F00')) +
  geom_hline(yintercept = 0, size = 0.3) +
  # scale_x_discrete(position = "top") +
  # scale_y_reverse() +
  # scale_x_discrete(position = "top") +
  # scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  labs(x = NULL, 
       y = "Log2 (FC)"
  ) +
  theme_bw(base_size = 28) +
  theme(    
    # axis.ticks = element_blank(), 
    legend.position = "nonoe",
    panel.grid = element_blank(),
    axis.text = element_text(size = 26, color = "black"),
    axis.text.x=element_text(angle=45, hjust=1, size = 25)
    # axis.text.y = element_text()
  )
dev.print("tcell_DE_genes_log2FC.pdf", width = 10, height = 4, dev = pdf)
dev.off()

