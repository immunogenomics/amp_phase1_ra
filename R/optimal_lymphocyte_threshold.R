#' ---
#' title: "Choose an optimal lymphocyte percent threshold"
#' author: "Kamil Slowikowski"
#' date: "2018-03-20"
#' ---

#+ setup, include=FALSE
library(knitr)
library(ggplot2)
library(mixtools)

opts_chunk$set(
  echo = TRUE
)
source("R/meta_colors.R")

#' Read the data tables:
#+ read-data
d1 <- readxl::read_excel("data-raw/Phase1/metadata/postQC_15OA_18RA_clin_flow.xlsx")
d2 <- readxl::read_excel("data-raw/Phase1/metadata/postQC_18bx_clin_flow.xlsx")

keep_cols <- c(
  "Patient", "Disease", "Sex", "Age", "Race", "Lymphocytes"
)

# d <- rbind(d1[,keep_cols], d2[,keep_cols])
d <- data.table::rbindlist(list(d1, d2))

d$Type <- ifelse(startsWith(d$Patient, "300"), "biopsy", "arthro")
d$disease_assay <- paste(d$Disease, d$Type)

#+ histogram
ggplot(d) +
  geom_histogram(aes(x = Lymphocytes.Live, fill = disease_tissue), bins = 30) +
  scale_fill_manual(values = meta_colors$tissue_type)

#' Fit two normal distributions:
#+ fit

fit <- normalmixEM(x = d$Lymphocytes)

sdnorm <- function(x, mean = 0, sd = 1, lambda = 1){
  lambda * dnorm(x, mean = mean, sd = sd)
}

dnorm2 <- function(fit) {
  return(function(x) {
    dnorm(x, mean = fit$mu[1], sd = fit$sigma[1]) * fit$lambda[1] -
      dnorm(x, mean = fit$mu[2], sd = fit$sigma[2]) * fit$lambda[2]
  })
}

root <- uniroot(dnorm2(fit), interval = c(0, 1))

d$Active <- d$Lymphocytes.Live > root$root
table(d$disease_tissue, d$Active)

ggplot(d) +
  geom_histogram(
    aes(x = Lymphocytes.Live, fill = disease_tissue),
    # bins = 30
    binwidth = 0.02
  ) +
  scale_fill_manual(
    values = meta_colors$tissue_type,
    name = NULL
  ) +
  # stat_function(
  #   fun = sdnorm,
  #   args = list(
  #     mean = fit$mu[1],
  #     sd = fit$sigma[1],
  #     lambda = fit$lambda[1]
  #   ),
  #   color = "black", fill = NA, geom = "polygon"
  # ) +
  # stat_function(
  #   fun = sdnorm,
  #   args = list(
  #     mean = fit$mu[2],
  #     sd = fit$sigma[2],
  #     lambda = fit$lambda[2]
  #   ),
  #   color = "black", fill = NA, geom = "polygon"
  # ) +
  # geom_vline(xintercept = root$root) +
  geom_vline(xintercept = 0.25) +
  annotate(
    geom = "text",
    # x = root$root,
    x = 0.25,
    y = 8,
    # label = signif(100 * root$root, 3),
    label = signif(100 * 0.25, 3),
    hjust = 0, size = 6
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_x_continuous(
    label = function(x) 100 * x,
    limits = c(-0.1, 0.8)
  ) +
  labs(x = "Lymphocytes\n(% of synovial cells)", y = "Count") +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.85, 0.8),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  )
ggsave(file = paste("postQC_samplse_lym_25", ".pdf", sep = ""), width = 7.5, height = 5, dpi = 200)
dev.off()

ggsave(
  filename = "figures/Phase1/optimal_lymphocyte_threshold.pdf",
  width = 6, height = 4
)

ggplot(d) +
  geom_histogram(aes(x = CD45pos.ov.Live, fill = disease_tissue), bins = 30) +
  scale_fill_manual(values = meta_colors$tissue_type)

d$CD45pos.ov.Live[is.na(d$CD45pos.ov.Live)] <- 0
d$CD45pos.ov.Live <- d$CD45pos.ov.Live / 100
fit <- normalmixEM(x = d$CD45pos.ov.Live, mu = c(0.3, 0.8))

root <- uniroot(dnorm2(fit), interval = c(0, 1))

d$Active2 <- d$CD45pos.ov.Live > root$root
table(d$disease_tissue, d$Active2)

ggplot(d) +
  geom_histogram(
    aes(x = CD45pos.ov.Live, fill = disease_tissue),
    # bins = 60
    binwidth = 0.018
  ) +
  scale_fill_manual(
    values = meta_colors$tissue_type,
    name = NULL
  ) +
  stat_function(
    fun = sdnorm,
    args = list(
      mean = fit$mu[1],
      sd = fit$sigma[1],
      lambda = fit$lambda[1]
    ),
    color = "black", fill = NA, geom = "polygon"
  ) +
  stat_function(
    fun = sdnorm,
    args = list(
      mean = fit$mu[2],
      sd = fit$sigma[2],
      lambda = fit$lambda[2]
    ),
    color = "black", fill = NA, geom = "polygon"
  ) +
  # geom_vline(xintercept = root$root) +
  geom_vline(xintercept = 0.50) +
  annotate(
    geom = "text",
    # x = root$root,
    x = 0.50,
    y = 8,
    # label = signif(100 * root$root, 3),
    label = signif(100 * 0.50, 3),
    hjust = -0.1, size = 6
  ) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_x_continuous(
    label = function(x) 100 * x,
    limits = c(-0.1, 1.1)
  ) +
  labs(x = "CD45+\n(% of synovial cells)", y = "Count") +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.85, 0.8),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  )
ggsave(file = paste("postQC_samplse_CD45_49", ".pdf", sep = ""), width = 7.5, height = 5, dpi = 200)
dev.off()

ggsave(
  filename = "figures/Phase1/optimal_cd45_threshold.pdf",
  width = 6, height = 4
)

ggplot(d) +
  geom_point(
    aes(x = Lymphocytes, y = CD45pos.ov.Live, fill = disease_assay),
    size = 4, shape = 21, stroke = 0.2
  ) +
  scale_fill_manual(
    values = meta_colors$disease_assay,
    name = NULL
  ) +
  labs(y = "CD45+", x = "Lymphocytes") +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.85, 0.3),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  )

ggplot(d) +
  geom_point(
    aes(x = B.cells, y = T.cells, fill = disease_assay),
    size = 4, shape = 21, stroke = 0.2
  ) +
  scale_fill_manual(
    values = meta_colors$disease_assay,
    name = NULL
  ) +
  labs(x = "B cells\n(% of synovial cells)", y = "T cells (%)") +
  scale_x_continuous(label = function(x) 100 * x) +
  scale_y_continuous(label = function(x) 100 * x) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = c(0.85, 0.8),
    legend.box.background = element_rect(size = 0.5),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5)
  )
ggsave(
  filename = "figures/Phase1/bcells_vs_tcells.pdf",
  width = 6, height = 4
)
