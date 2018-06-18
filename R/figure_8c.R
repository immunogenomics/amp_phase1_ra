
setwd("/Users/fanzhang/Documents/GitHub/amp_phase1_ra/R/")

source("meta_colors.R")
library(gdata) 
library(ggplot2)
library(plyr)
library(dplyr)

# Read post-QC single-cell RNA-seq meta data
log2cpm <- readRDS("../data/celseq_synovium_log2_5452cells_paper.rds")
meta <- readRDS("../data/celseq_synovium_meta_5452cells_paper.rds")
meta$cell_name <- as.character(meta$cell_name)

log2cpm <- log2cpm[, -which(meta$fine_cluster == "T-1")]
meta <- meta[-which(meta$fine_cluster == "T-1"),]
all(meta$cell_name == colnames(log2cpm))
log2cpm <- as.data.frame(log2cpm)

clusters <- c(
  "F-1", "F-2", "F-3", "F-4",
  "M-1",  "M-2", "M-3", "M-4",
  "T-2", "T-3", "T-4", "T-5", "T-6", "T-7",
  "B-1", "B-2", "B-3", "B-4"
)
meta$fine_cluster <- factor(meta$fine_cluster, levels = clusters)
fontsize <- 30

# -----
marker <- "IL6"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "IL6"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p1 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


# -----
marker <- "IL1B"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "IL1B"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p2 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "CXCL9"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "CXCL9"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p3 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


# -----
marker <- "CXCL13"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "CXCL13"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p4 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


# -----
marker <- "TNF"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "TNF"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p5 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


# -----
marker <- "TNFRSF1A"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "TNFRSF1A"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p6 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


# -----
marker <- "IFNG"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "IFNG"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p7 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "IRF8"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "IRF8"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p8 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


# -----
marker <- "SLAMF7"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "SLAMF7"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p9 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "TLR10"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "TLR10"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p10 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "TLR8"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "TLR8"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p11 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "PTGS2"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "PTGS2"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p12 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "RSAD2"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "RSAD2"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p13 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )

# -----
marker <- "IFITM3"
meta$marker <- as.numeric(log2cpm[which(rownames(log2cpm) == "IFITM3"),])
dat_percent <- meta %>%
  group_by(fine_cluster) %>%
  summarise(percent = sum(marker > 0) / length(marker) * 100)

p14 <- ggplot() +
  geom_col(
    data = dat_percent,
    mapping = aes(x = fine_cluster, y = percent, fill = fine_cluster),
    color = "grey50", size = 0.1
  ) +
  scale_fill_manual(values = meta_colors$fine_cluster, name = "Cluster") +
  scale_x_discrete(limits = rev(levels(meta$fine_cluster))) +
  # scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  geom_vline(xintercept = 14.5, linetype="dashed") +
  geom_vline(xintercept = 10.5, linetype="dashed") +
  geom_vline(xintercept = 4.5, linetype="dashed") +
  coord_flip() +
  labs(x = NULL, y = NULL, title = marker) +
  theme_bw(base_size = fontsize) + theme(
    legend.position = "none",
    axis.text.y       = element_blank(),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(size = 0.5),
    plot.title = element_text(color="black", size=fontsize, face = "italic")
  )


g <- egg::ggarrange(
  p1, p2, p3, p4, p5, p6,
  p7, p8, p9, p10, p11, p12,
  ncol = 6
)
ggsave("figure8c.pdf", g, height = 8, width = 22, dpi = 300)
dev.off()

g <- egg::ggarrange(
  p13, p14,
  ncol = 2
)
ggsave("INF.pdf", g, height = 4.5, width = 8, dpi = 300)
dev.off()

