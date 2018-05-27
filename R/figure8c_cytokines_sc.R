#' ---
#' title: "Figure 8c"
#' author: "Fan Zhang"
#' date: "2018-05-24"
#' ---

setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

source("meta_colors.R")
pacman::p_load(
  ggplot2,
  patchwork,
  dplyr,
  magrittr,
  cetcolor,
  seriation,
  data.table,
  pheatmap
)

# Read post-QC single-cell RNA-seq meta data
log2cpm <- readRDS("../data/celseq_synovium_log2_5452cells_paper.rds")
meta <- readRDS("../data/celseq_synovium_meta_5452cells_paper.rds")
meta <- meta[order(meta$cell_name),]
meta$cell_name <- as.character(meta$cell_name)
all(colnames(log2cpm) == meta$cell_name)

genes <- read.xls("../../../HMS/amp/results/2017_01_12_Low_input_RA_RNA_seq_data_analysis/Figure3B_gene_list_order_remove_IFNB1.xlsx")
log2cpm_plot <- log2cpm[which(rownames(log2cpm) %in% c(as.character(genes$gene), "TLR10", "TLR8")), ]
meta_plot <- meta[which(meta$cell_name %in% colnames(log2cpm_plot)), ]
all(colnames(log2cpm_plot) == meta_plot$cell_name)
all <- cbind.data.frame(meta_plot, t(log2cpm_plot))
dim(all)
all <- all[, c(35, 37:78)]


gene_names <- colnames(all)[-1]
all_nonzero_perc <- list()
for (x in c(1:length(gene_names))){
  print(x)
  temp <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(cat(gene_names[x], "\n") != 0))
  all_nonzero_perc <- append(all_nonzero_perc, list(temp))
}
# all %>% group_by(fine_cluster) %>% summarise(sum(cat(gene_names[x], "\n") != 0))
# test <- apply(all, 2, function(x) all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(x != 0)) ) 


gene_IL1B <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IL1B != 0))
gene_IL1B$gene <- rep("IL1B", length(table(all$fine_cluster)))
gene_IL1B$nonzero_perc <- (gene_IL1B$nonzero_count/sum(gene_IL1B$nonzero_count)) * 100
gene_IL1B <- as.data.frame(gene_IL1B)

gene_IL15 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IL15 != 0))
gene_IL15$gene <- rep("IL15", length(table(all$fine_cluster)))
gene_IL15$nonzero_perc <- (gene_IL15$nonzero_count/sum(gene_IL15$nonzero_count)) * 100
gene_IL15 <- as.data.frame(gene_IL15)

gene_CXCL11 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(CXCL11 != 0))
gene_CXCL11$gene <- rep("CXCL11", length(table(all$fine_cluster)))
gene_CXCL11$nonzero_perc <- (gene_CXCL11$nonzero_count/sum(gene_CXCL11$nonzero_count)) * 100
gene_CXCL11 <- as.data.frame(gene_CXCL11)

gene_NFKB1 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(NFKB1 != 0))
gene_NFKB1$gene <- rep("NFKB1", length(table(all$fine_cluster)))
gene_NFKB1$nonzero_perc <- (gene_NFKB1$nonzero_count/sum(gene_NFKB1$nonzero_count)) * 100
gene_NFKB1 <- as.data.frame(gene_NFKB1)

gene_PTPN2 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(PTPN2 != 0))
gene_PTPN2$gene <- rep("PTPN2", length(table(all$fine_cluster)))
gene_PTPN2$nonzero_perc <- (gene_PTPN2$nonzero_count/sum(gene_PTPN2$nonzero_count)) * 100
gene_PTPN2 <- as.data.frame(gene_PTPN2)

gene_PTGS2 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(PTGS2 != 0))
gene_PTGS2$gene <- rep("PTGS2", length(table(all$fine_cluster)))
gene_PTGS2$nonzero_perc <- (gene_PTGS2$nonzero_count/sum(gene_PTGS2$nonzero_count)) * 100
gene_PTGS2 <- as.data.frame(gene_PTGS2)

gene_ICAM1 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(ICAM1 != 0))
gene_ICAM1$gene <- rep("ICAM1", length(table(all$fine_cluster)))
gene_ICAM1$nonzero_perc <- (gene_ICAM1$nonzero_count/sum(gene_ICAM1$nonzero_count)) * 100
gene_ICAM1 <- as.data.frame(gene_ICAM1)

gene_IFIT2 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IFIT2 != 0))
gene_IFIT2$gene <- rep("IFIT2", length(table(all$fine_cluster)))
gene_IFIT2$nonzero_perc <- (gene_IFIT2$nonzero_count/sum(gene_IFIT2$nonzero_count)) * 100
gene_IFIT2 <- as.data.frame(gene_IFIT2)

gene_RSAD2 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(RSAD2 != 0))
gene_RSAD2$gene <- rep("RSAD2", length(table(all$fine_cluster)))
gene_RSAD2$nonzero_perc <- (gene_RSAD2$nonzero_count/sum(gene_RSAD2$nonzero_count)) * 100
gene_RSAD2 <- as.data.frame(gene_RSAD2)

gene_STAT1 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(STAT1 != 0))
gene_STAT1$gene <- rep("STAT1", length(table(all$fine_cluster)))
gene_STAT1$nonzero_perc <- (gene_STAT1$nonzero_count/sum(gene_STAT1$nonzero_count)) * 100
gene_STAT1 <- as.data.frame(gene_STAT1)

gene_CXCL9 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(CXCL9 != 0))
gene_CXCL9$gene <- rep("CXCL9", length(table(all$fine_cluster)))
gene_CXCL9$nonzero_perc <- (gene_CXCL9$nonzero_count/sum(gene_CXCL9$nonzero_count)) * 100
gene_CXCL9 <- as.data.frame(gene_CXCL9)

gene_CXCL12 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(CXCL12 != 0))
gene_CXCL12$gene <- rep("CXCL12", length(table(all$fine_cluster)))
gene_CXCL12$nonzero_perc <- (gene_CXCL12$nonzero_count/sum(gene_CXCL12$nonzero_count)) * 100
gene_CXCL12 <- as.data.frame(gene_CXCL12)

gene_CX3CL1 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(CX3CL1 != 0))
gene_CX3CL1$gene <- rep("CX3CL1", length(table(all$fine_cluster)))
gene_CX3CL1$nonzero_perc <- (gene_CX3CL1$nonzero_count/sum(gene_CX3CL1$nonzero_count)) * 100
gene_CX3CL1 <- as.data.frame(gene_CX3CL1)

gene_IL6 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IL6 != 0))
gene_IL6$gene <- rep("IL6", length(table(all$fine_cluster)))
gene_IL6$nonzero_perc <- (gene_IL6$nonzero_count/sum(gene_IL6$nonzero_count)) * 100
gene_IL6 <- as.data.frame(gene_IL6)

gene_PDGFRB <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(PDGFRB != 0))
gene_PDGFRB$gene <- rep("PDGFRB", length(table(all$fine_cluster)))
gene_PDGFRB$nonzero_perc <- (gene_PDGFRB$nonzero_count/sum(gene_PDGFRB$nonzero_count)) * 100
gene_PDGFRB <- as.data.frame(gene_PDGFRB)

gene_SLAMF7 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(SLAMF7 != 0))
gene_SLAMF7$gene <- rep("SLAMF7", length(table(all$fine_cluster)))
gene_SLAMF7$nonzero_perc <- (gene_SLAMF7$nonzero_count/sum(gene_SLAMF7$nonzero_count)) * 100
gene_SLAMF7 <- as.data.frame(gene_SLAMF7)

gene_IRF8 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IRF8 != 0))
gene_IRF8$gene <- rep("IRF8", length(table(all$fine_cluster)))
gene_IRF8$nonzero_perc <- (gene_IRF8$nonzero_count/sum(gene_IRF8$nonzero_count)) * 100
gene_IRF8 <- as.data.frame(gene_IRF8)

gene_IRF7 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IRF7 != 0))
gene_IRF7$gene <- rep("IRF7", length(table(all$fine_cluster)))
gene_IRF7$nonzero_perc <- (gene_IRF7$nonzero_count/sum(gene_IRF7$nonzero_count)) * 100
gene_IRF7 <- as.data.frame(gene_IRF7)

gene_IRF1 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IRF1 != 0))
gene_IRF1$gene <- rep("IRF1", length(table(all$fine_cluster)))
gene_IRF1$nonzero_perc <- (gene_IRF1$nonzero_count/sum(gene_IRF1$nonzero_count)) * 100
gene_IRF1 <- as.data.frame(gene_IRF1)

gene_IRF9 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IRF9 != 0))
gene_IRF9$gene <- rep("IRF9", length(table(all$fine_cluster)))
gene_IRF9$nonzero_perc <- (gene_IRF9$nonzero_count/sum(gene_IRF9$nonzero_count)) * 100
gene_IRF9 <- as.data.frame(gene_IRF9)

gene_IRF3 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(IRF3 != 0))
gene_IRF3$gene <- rep("IRF3", length(table(all$fine_cluster)))
gene_IRF3$nonzero_perc <- (gene_IRF3$nonzero_count/sum(gene_IRF3$nonzero_count)) * 100
gene_IRF3 <- as.data.frame(gene_IRF3)

gene_HSPD1 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(HSPD1 != 0))
gene_HSPD1$gene <- rep("HSPD1", length(table(all$fine_cluster)))
gene_HSPD1$nonzero_perc <- (gene_HSPD1$nonzero_count/sum(gene_HSPD1$nonzero_count)) * 100
gene_HSPD1 <- as.data.frame(gene_HSPD1)

gene_CXCL13 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(CXCL13 != 0))
gene_CXCL13$gene <- rep("CXCL13", length(table(all$fine_cluster)))
gene_CXCL13$nonzero_perc <- (gene_CXCL13$nonzero_count/sum(gene_CXCL13$nonzero_count)) * 100
gene_CXCL13 <- as.data.frame(gene_CXCL13)

gene_TNF <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(TNF != 0))
gene_TNF$gene <- rep("TNF", length(table(all$fine_cluster)))
gene_TNF$nonzero_perc <- (gene_TNF$nonzero_count/sum(gene_TNF$nonzero_count)) * 100
gene_TNF <- as.data.frame(gene_TNF)

gene_TLR10 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(TLR10 != 0))
gene_TLR10$gene <- rep("TLR10", length(table(all$fine_cluster)))
gene_TLR10$nonzero_perc <- (gene_TLR10$nonzero_count/sum(gene_TLR10$nonzero_count)) * 100
gene_TLR10 <- as.data.frame(gene_TLR10)

gene_TLR8 <- all %>% group_by(fine_cluster) %>% summarise(nonzero_count = sum(TLR8 != 0))
gene_TLR8$gene <- rep("TLR8", length(table(all$fine_cluster)))
gene_TLR8$nonzero_perc <- (gene_TLR8$nonzero_count/sum(gene_TLR8$nonzero_count)) * 100
gene_TLR8 <- as.data.frame(gene_TLR8)


mat <- rbind.data.frame(gene_IL15, gene_NFKB1, gene_PTPN2, gene_PTGS2, gene_ICAM1,
                        gene_IFIT2, gene_RSAD2, gene_STAT1, gene_CXCL9, gene_CXCL12, gene_CX3CL1,
                        gene_IL6, gene_PDGFRB,  gene_SLAMF7, gene_IRF8, gene_IRF7, gene_IRF1, gene_IRF9,
                        gene_IRF3, gene_HSPD1, gene_TNF, gene_TLR10, gene_TLR8
                        )

seriate_cols <- function(d, col1, col2, value.var = "percent") {
  mat <- dcast(
    data = d,
    formula = as.formula(sprintf("%s ~ %s", col1, col2)),
    value.var = value.var
  )
  rownames(mat) <- mat[[1]]
  mat[[1]] <- NULL
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat_order <- seriate(mat, method = "BEA_TSP")
  d[[col1]] <- factor(as.character(d[[col1]]), rownames(mat)[mat_order[[1]]])
  d[[col2]] <- factor(as.character(d[[col2]]), colnames(mat)[mat_order[[2]]])
  return(d)
}

# mat$fine_cluster = factor(mat$fine_cluster, levels=c('F-1', 'F-2', 'F-3', 'F-4',
#                                                      'M-1', 'M-2', 'M-3', 'M-4',
#                                                      'B-1', 'B-2', 'B-3', 'B-4',
#                                                      'T-1', 'T-2', 'T-3', 'T-4', 'T-5', 'T-6', 'T-7'))
this_d <- seriate_cols(mat, "gene", "fine_cluster","nonzero_count") # "nonzero_count", "nonzero_perc" 
this_d %<>% arrange(fine_cluster)

ggplot(
  data=this_d,
  mapping = aes(x=gene, y= nonzero_count, fill = fine_cluster)
  ) +
  geom_col() +
  geom_bar(stat="identity",
           position = "fill"
           # position = "stack",
           # width = 0.15
  ) +
  coord_flip() +
  scale_fill_manual("", values = meta_colors$fine_cluster) 


