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


# Read Case.Control labels for OA, non-inflamed RA, and RA
inflam_label <- read.xls("../data-raw/postQC_all_samples.xlsx")
inflam_label <- inflam_label[c(1:51),]
inflam_label <- inflam_label[, c("Patient", "Case.Control")]
table(inflam_label$Case.Control)
colnames(inflam_label)[1] <- "Donor.ID"

clin_merge <- merge(clin_meta, inflam_label, by = "Donor.ID")
table(clin_merge$Case.Control)

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

top_percent <- 0.7
mat_idx <- with(
  dat_mean_sd,
  mean > quantile(mean, top_percent) &
    sd > quantile(sd, top_percent)
)

all_samples <- log2tpm[mat_idx,]
dim(all_samples)
all(clin_merge$Sample.ID == colnames(all_samples))


# ----------------------------------------------------------------------------------------
# fibroblast samples
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
bulk_samples <- bulk_samples[, -which(bulk_m$Case.Control == "non-inflamed RA")]
bulk_m <- bulk_m[-which(bulk_m$Case.Control == "non-inflamed RA"), ]
bulk_m$Case.Control <- as.character(bulk_m$Case.Control)
all(bulk_m$Sample.ID == colnames(bulk_samples))

# Use Limma R function by fitting data access group and plates as covariates 
des <- with(
  bulk_m,
  # model.matrix(~ Case.Control + Data.Access.Group + cDNA.plate)
  model.matrix(~ Case.Control)
)

# Fit the model to each gene
fit <- lmFit(object = bulk_samples, design = des)

# Share variance across genes
fit <- eBayes(fit)

# BH FDR control is used
toptable_fdr <-topTable(fit, coef = "Case.ControlOA", 
                        number=nrow(bulk_samples), adjust.method="BH", 
                        confint=TRUE,
                        sort.by="p")
toptable_fdr[1:4,]

outFile = "fibro_diff_inflamedRA_OA.txt"
write.table(toptable_fdr,outFile,row.names=T,col.names=T,quote=F, sep = "\t")
# toptable_fdr <- read.table("fibro_diff_inflamedRA_OA.txt")

# Show CCA genes (also the single-cell RNA-seq cluster markers)
dat_plot <- toptable_fdr[which(rownames(toptable_fdr) %in% 
                    c("HBEGF", "CLIC5", "HTRA1", "PRG4", "CD55", "DNASE1L3",
                      "PTGFR", "FOS", "F3", "HLA.DRA", "HLA.DPA1", "HLA.DRB1", "IL6", # "IFI30", 
                      "DKK3", "CADM1", "CAPG", # "AKR1C2", "COL8A2",
                      "ACTA2", "CD34", # "MCAM", "MYH11"
                      "IRF1", "CXCL12",
                      "PDGFRB"
                   )),]

dat_plot$gene <- rownames(dat_plot)
dat_plot$logFC <- (-1) * dat_plot$logFC
dat_plot$CI.L <- (-1) * dat_plot$CI.L
dat_plot$CI.R <- (-1) * dat_plot$CI.R


# Plot genes based on logFC
ggplot() +
  geom_errorbar(data=dat_plot, 
                mapping=aes(x=reorder(gene, logFC), ymin=CI.L, ymax=CI.R), 
                width=0.3, size=0.6, color = "black"
                ) +
  geom_point(data=dat_plot, 
            aes(x=reorder(gene, logFC),y= logFC),
            size = 1
  ) +
  # scale_color_manual(values = c('grey70', 'black')) +
  geom_hline(yintercept = 0, size = 0.3) +
  scale_x_discrete(position = "top") +
  #scale_y_continuous(labels = function(x) round(2^abs(x), 1)) +
  coord_flip() +
  labs(x = NULL, 
       y = "Log2 FC",
       subtitle = "Inflamed RA versus OA"
       # title = "abs(Log2 FC) > 0.5 and 5% FDR"
       ) +
  theme_bw(base_size = 26) +
  theme(    
    # axis.ticks = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_text(size = 26, color = "black"),
    axis.text.y = element_text(size = 26, face="italic"))
dev.print("fibro_DE_genes_log2FC.pdf", width = 6, height = 10, dev = pdf)
dev.off()



