# Read top n markers for each identified clusters
# These markers are used for heritability analysis
# Fan Zhang

dat_table <- readRDS("../data/cluster_marker_table.rds")
table(dat_table$cluster)
test <- dat_table$cell_type
test[which(test != "All cells")] <- "One cell type cells"
dat_table$test <- test

# "All cells" are used for comparing one cluster versus all (18) the other clusters
dat_table <- dat_table[-which(dat_table$cell_type == "All cells"),]
dim(dat_table)

# Remove the mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = dat_table$gene, value = TRUE)
mito.genes_2 <- grep(pattern = "^MTRNR", x = dat_table$gene, value = TRUE)
dat_table <- dat_table[-which(dat_table$gene %in% c(mito.genes_1, mito.genes_2)), ] 
dim(dat_table)

# Take top 20 markers for each cluster
x <- dat_table %>%
  group_by(cluster) %>%
  top_n(20, wt = auc)

x$cluster = factor(x$cluster, 
                   levels=c("C-F1", "C-F2", "C-F3", "C-F4",
                            "C-M1", "C-M2", "C-M3", "C-M4",
                            "C-B1", "C-B2", "C-B3", "C-B4",
                            "C-T1", "C-T2", "C-T3", "C-T4", "C-T5", "C-T6", "C-T7"
                   ))

# Sort the genes based on auc score
x <- x[order(x$auc, decreasing = TRUE),]
x <- x[order(x$cluster, decreasing = FALSE),]
x <- as.data.frame(x)
dim(x)

