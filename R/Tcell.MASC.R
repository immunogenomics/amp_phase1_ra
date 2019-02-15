library(lme4)
library(useful)
library(readxl)
library(plyr)

load("synData.Tcell.downsample.SNE.RData")

source("graphFn.R")


### ADD QC-ED METADATA
metaDF <- read_excel("postQC_all_samples.xlsx", sheet = "R_import")
metaDF$sampleID <- as.character(metaDF$Patient)
# Check that all samples in dataset have associated metadata
allSampleID <- unique(as.character(synData.Tcell.downsample$sampleID))
stopifnot(all(allSampleID %in% unique(as.character(metaDF$Patient))))
# Join dataframes
synData.Tcell.downsample <- join(synData.Tcell.downsample, metaDF, by = "sampleID")
# Relabel status.tissue to match Case/Control column (samples do not change between groups)
synData.Tcell.downsample$case <- factor(synData.Tcell.downsample$status.tissue,
                                        levels = c("OA.Arthoplasty", "RA.Arthoplasty", "RA.Biopsy"),
                                        labels = c("OA", "non-inflamed RA", "inflamed RA"))
synData.Tcell.downsample$case.2way <- factor(ifelse(synData.Tcell.downsample$case == "inflamed RA", "leuko-high", "leuko-low"),
                                             levels = c("leuko-low", "leuko-high"),
                                             labels = c("Leukocyte-poor RA and OA", "Leukocyte-rich RA"))
tcell_mrg <- readRDS("tcell_merg.rds")
tcell_mrg$assign_lbl <- LETTERS[1:nrow(tcell_mrg)]
# Give clusters assignment labels based on merging strategy
synData.Tcell.downsample$assign_lbl <- recode_factor(synData.Tcell.downsample$SNE.cluster,
                                                     "1" = "I", "2" = "G", "3" = "I",
                                                     "4" = "I", "5" = "A", "6" = "H",
                                                     "7" = "G", "8" = "F", "9" = "H",
                                                     "10" = "F", "11" = "F", "12" = "H",
                                                     "13" = "G", "14" = "B", "15" = "D",
                                                     "16" = "D", "17" = "C", "18" = "E",
                                                     "19" = "D", "20" = "D", "21" = "E")
synData.Tcell.downsample <- join(synData.Tcell.downsample, tcell_mrg)

###############################################################################
# 1) Plot SNE
###############################################################################
cluster.centers <- t(sapply(split(synData.Tcell.downsample, synData.Tcell.downsample$markers), function(x) c(mean(x[["SNE1"]]), mean(x[["SNE2"]]))))
cluster.centers <- data.frame(cluster = rownames(cluster.centers),
                              SNE1 = cluster.centers[,1], 
                              SNE2 = cluster.centers[,2])
cluster.centers[2, c("SNE1", "SNE2")] <- c(-25, -4)
cluster.centers[5, c("SNE1", "SNE2")] <- c(-10, -26)
cluster.centers[6, c("SNE1", "SNE2")] <- c(24, 27)
cluster.centers[7, c("SNE1", "SNE2")] <- c(11, 19)
cluster.centers[9, c("SNE1", "SNE2")] <- c(22, -19)

cluster.colors <- as.vector(ggthemes_data$tableau$colors$tableau10medium)
p <- ggplot(synData.Tcell.downsample, aes(x = SNE1, y = SNE2)) + geom_point(aes(color = markers), size = 0.25)
p <- p + coord_fixed(xlim = c(-50, 50), ylim = c(-50, 50)) + GenericTheme + sneTheme
p <- p + geom_label(data = cluster.centers, aes(label = cluster), fill = "white", color = "black", size = 4, fontface = "bold")
p <- p + scale_color_manual(values = cluster.colors, guide = F) + scale_fill_manual(values = cluster.colors, guide = F)
p 
ggsave(p, filename = "output/Tcell.SNE.clusters.merged.pdf", width = 5, height = 5, useDingbats = FALSE)


###############################################################################
# 2) Test for  leuko-rich sample differential abundance with MASC
###############################################################################
nCluster <- nlevels(droplevels(synData.Tcell.downsample$assign_lbl))
tcell.clusters.df <- data.frame(cluster = 1:nCluster,
                                markers = levels(droplevels(synData.Tcell.downsample$markers)),
                                nCells = sapply(split(synData.Tcell.downsample, synData.Tcell.downsample$markers, drop = T), nrow))
tcell.clusters.df <- join(tcell.clusters.df, tcell_mrg)
cluster.groupPct <- t(sapply(split(synData.Tcell.downsample, synData.Tcell.downsample$markers), function(x) 100 * prop.table(table(x$case))))
colnames(cluster.groupPct) <- c("pct_OA", "pct_non.inflamed.RA", "pct_inflamed.RA")
tcell.clusters.df <- cbind(tcell.clusters.df, as.data.frame(cluster.groupPct))
# Logistic mixed-effect model association testing
modelData <- synData.Tcell.downsample[c("sampleID", "runDate", "case.2way", "Age", "Sex", "assign_lbl")]
# Relevel case variable
# modelData$case2 <- factor(ifelse(modelData$case == "inflamed RA", "inflamed", "non-inflamed"),
#                           levels = c("non-inflamed", "inflamed"))
# Rescale age variable
# modelData$Age.sc <- scale(modelData$Age)

for (i in 1:nCluster) {
  name <- paste("cluster", i, sep = "")
  modelData[[name]] <- factor(ifelse(modelData$assign_lbl == tcell.clusters.df[i, "assign_lbl"], 1, 0))
}

lmmData <- vector(length = nCluster, mode = "list")
names(lmmData) <- paste("cluster", 1:nCluster, sep = "")


for (i in 1:nCluster) {
  message(paste("Creating logistic mixed models for cluster", i))
  null.model <- as.formula(paste("cluster", i, " ~ 1 + (1|sampleID) + (1|runDate)", sep = ""))
  full.model <- as.formula(paste("cluster", i, " ~ case.2way + (1|sampleID)+ (1|runDate)", sep = ""))
  lmmData[[i]]$null <- glmer(formula = null.model, data = modelData, family = binomial,
                             control = glmerControl(optimizer = "bobyqa") ,nAGQ = 1,  verbose = 1)
  lmmData[[i]]$full <- glmer(formula = full.model, data = modelData, family = binomial,
                             control = glmerControl(optimizer = "bobyqa") ,nAGQ = 1,  verbose = 1)
  lmmData[[i]]$anova <- anova(lmmData[[i]]$null, lmmData[[i]]$full)
  # calculate confidence intervals for the case-control beta
  lmmData[[i]]$confint <- confint.merMod(lmmData[[i]]$full, parm = "case.2wayleuko-high", method = "profile")
}

# Add logistic mixed-effect model p-values and odds ratios for each cluster
tcell.clusters.df$masc.pval <- sapply(lmmData, function(x) x$anova[["Pr(>Chisq)"]][2])
tcell.clusters.df$fdr.pval <- p.adjust(tcell.clusters.df$masc.pval, method = "fdr")
tcell.clusters.df$leukohigh.pval <- sapply(lmmData, function(x) coef(summary(x$full))["case.2wayleuko-high", "Pr(>|z|)"])
tcell.clusters.df$leukohigh.or <- sapply(lmmData, function(x) exp(fixef(x$full)["case.2wayleuko-high"]))
tcell.clusters.df$leukohigh.or.lower <- sapply(lmmData, function(x) exp(x$confint["case.2wayleuko-high", "2.5 %"]))
tcell.clusters.df$leukohigh.or.upper <- sapply(lmmData, function(x) exp(x$confint["case.2wayleuko-high", "97.5 %"]))

# Signifcant clusters
tcell.clusters.df[tcell.clusters.df$fdr.pval < 0.05,]

write.table(tcell.clusters.df, file = pipe("pbcopy"), col.names = T, row.names = F, quote = F, sep = "\t")
