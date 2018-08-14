# plot_data_overview.R
# Kamil Slowikowski
# 2017-08-30
#
# Make a figure that summarizes which types of data we have for each patient.
#
# This figure only includes samples that pass quality-control filters.
#

pacman::p_load(
  dplyr,
  reshape2,
  readr,
  readxl,
  cowplot,
  ggplot2,
  gridExtra,
  scales,
  sitools,
  stringr,
  dendsort
)

# 2017-08-29
# Kevin Wei says:
# 300-0134 should be deleted, did not get enough cells from this person

# 2017-08-30
# Fan Zhang says:
#
# "300-0133" "300-0211" "300-0212" "300-0452" "301-0241" these are the failed QC patients
#
# Before QC:
#   /data/srlab/fzhang/amp/results/2017_02_09_flow_cytometry_clinical_data/AMP_lowInput_masterMetadata_020617.txt
#
# After QC:
#   /data/srlab/fzhang/amp/results/2017_06_26_CCA_bulk_singlecell/filtered_meta_lowinput_phase_1.rds
#
# This is the after QC log2tpm matrix:
#   /srlab/fzhang/amp/results/2017_06_26_CCA_bulk_singlecell/filtered_log2tpm_lowinput_phase_1.rds

source("R/meta_colors.R")

# Low-input RNA-seq and single-cell RNA-seq
d1 <- read_excel(
  path = "data-raw/Phase1/metadata/2017-07-19_AMP-Phase1-Sample-Accounting.xlsx",
  sheet = 2
)
d1 <- janitor::clean_names(d1)

# CyTOF
d2 <- read_excel(
  path = "data-raw/Phase1/metadata/2017-07-19_AMP-Phase1-Sample-Accounting.xlsx",
  sheet = 3
)
d2 <- janitor::clean_names(d2)
colnames(d2)[3:6] <- sprintf("%s_cytof", colnames(d2)[3:6])

# Flow cytometry info
d3 <- read_csv(
  # file = "data-raw/Phase1/metadata/2017-07-28_RA-flow-cytometry-counts.csv",
  file = "data-raw/Phase1/metadata/2017-08-29_AMP-RA-flow-cell-count.csv",
  col_types = cols(
    Sample = col_character(),
    Tcell = col_integer(),
    Bcell = col_double(),
    Monocyte = col_double(),
    Fibroblast = col_double(),
    `Other cells` = col_double(),
    `All viable cells` = col_integer()
  )
)
d3 <- janitor::clean_names(d3)
colnames(d3)[2:7] <- sprintf("%s_flow", colnames(d3)[2:7])

# Joint of origin, drugs
d4 <- read_excel(
  path = "data-raw/Phase1/metadata/2017-04-21_AMP-RA-Arthroplasty.xlsx",
  sheet = 1
)
d4 <- janitor::clean_names(d4)
d5 <- read_excel(
  path = "data-raw/Phase1/metadata/2017-04-21_AMP-RA-Biopsy.xlsx",
  sheet = 1
)
d5 <- janitor::clean_names(d5)
d5 <- subset(d5, visit_idx_y_z_a == "WK 0")

joints <- c(
  setNames(d4$operative_joint, d4$clinical_barcode),
  setNames(d5$biopsy_site, d5$clinical_barcode)
)
x <- str_split_fixed(joints, " ", 2)
x[,2][x[,2] == ""] <- x[,1][x[,2] == ""]
x[,2] <- tolower(x[,2])
joints <- setNames(x[,2], names(joints))

# Single-cell RNA-seq
d6 <- read_csv(
  file = "data-raw/Phase1/metadata/2017-06-12_AMP-RA-singlecells.csv",
  col_types = cols(
    sample = col_character(),
    disease_assay = col_character(),
    type = col_character(),
    all_cells = col_integer(),
    good_cells = col_integer()
  )
)
d6 <- janitor::clean_names(d6)
d6 <- dcast(
  data = d6, formula = sample + disease_assay ~ type, value.var = "good_cells",
  fill = 0
)
d6[["Empty"]] <- NULL
d6 <- subset(d6, sample != "none")

d <- left_join(d1, d2, by = c("subject_id" = "sample"))
d <- left_join(d, d3, by = c("subject_id" = "sample"))

d <- d[!is.na(d$subject_id),]
d <- d[d$sample_type != "Blood",]
d$joint <- joints[d$subject_id]
d$joint[d$joint == ""] <- "unknown"

d$sex <- c(
  setNames(d4$sex, d4$clinical_barcode),
  setNames(d5$sex, d5$clinical_barcode)
)[d$subject_id]

d$disease_assay <- paste(d$diagnosis, d$sample_type)
d$disease_assay[d$disease_assay == "OA Arthroplasty"] <- "OA Arthro."

dd <- d[,c(
  "subject_id", "disease_assay", "joint", "sex",
  "t_cells_for_single_cell",
  "b_cells_for_single_cell",
  "monoytes_for_single_cell",
  "fibroblasts_for_single_cell",
  "low_input_fibroblasts",
  "low_input_monoytes",
  "low_input_t_cells",
  "low_input_b_cells",
  "tcell_cytof",
  "bcell_cytof",
  "monocyte_cytof",
  "fibroblast_cytof",
  "tcell_flow",
  "bcell_flow",
  "monocyte_flow",
  "fibroblast_flow"
)]
# Collapse duplicate subjects.
dd <- dd %>% group_by(subject_id, disease_assay, joint, sex) %>% summarise_all(.funs = sum)
dd <- as.data.frame(dd)
rownames(dd) <- dd$subject_id

sorted_hclust <- function(...) as.hclust(dendsort(as.dendrogram(hclust(...))))

xs <- c("OA Arthro.", "RA Arthroplasty", "RA Biopsy")
ordered_subjects <- unlist(lapply(xs, function(x) {
  dd_temp <- subset(dd, disease_assay == x)
  hc <- sorted_hclust(dist(dd_temp[,5:ncol(dd)]))
  rownames(dd_temp)[hc$order]
}))

dd$subject_id <- factor(dd$subject_id, ordered_subjects)

x <- table(dd$disease_assay)
dd$disease_assay <- sprintf(
  "%s\n(n=%s)", dd$disease_assay, x[dd$disease_assay]
)

# Why is there no flow data for 301-0161 301-0241 300-0134
# dd[dd$subject_id == "301-0161",]
# dd[dd$subject_id == "301-0247",]
# dd[dd$subject_id == "301-0154",]
# dd[dd$subject_id == "300-0134",]

# Flow cytometry
# -----------------------------------------------------------------------------
md_flow <- melt(
  data = dd,
  id.vars = c("subject_id", "disease_assay"),
  measure.vars = c("tcell_flow", "bcell_flow", "monocyte_flow", "fibroblast_flow")
)
md_flow$variable <- c(
  "tcell_flow" = "T cell",
  "bcell_flow" = "B cell",
  "monocyte_flow" = "Monocyte",
  "fibroblast_flow" = "Fibroblast"
)[md_flow$variable]
md_flow$variable <- factor(
  md_flow$variable, c("B cell", "T cell", "Fibroblast", "Monocyte")
)
p_flow <- ggplot(
  md_flow, aes(x = subject_id, weight = value, fill = variable)
) +
  geom_bar() +
  # facet_grid(Type ~ DiseaseAssay, scale = "free", space = "free_x") +
  facet_grid(variable ~ disease_assay, scale = "free", space = "free_x") +
  scale_y_log10(labels = f2si) +
  scale_fill_manual(values = meta_colors$type) +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Flow Cytometry")

# CyTOF
# -----------------------------------------------------------------------------
md_cytof <- melt(
  data = dd,
  id.vars = c("subject_id", "disease_assay"),
  measure.vars = c("tcell_cytof", "bcell_cytof", "monocyte_cytof", "fibroblast_cytof")
)
# md_cytof$value[is.na(md_cytof$value)] <- 0.1
md_cytof$variable <- c(
  "tcell_cytof" = "T cell",
  "bcell_cytof" = "B cell",
  "monocyte_cytof" = "Monocyte",
  "fibroblast_cytof" = "Fibroblast"
)[md_cytof$variable]
md_cytof$variable <- factor(
  md_cytof$variable, c("B cell", "T cell", "Fibroblast", "Monocyte")
)
p_cytof <- ggplot(
  md_cytof, aes(x = subject_id, weight = value, fill = variable)
) +
  geom_bar() +
  # facet_grid(Type ~ DiseaseAssay, scale = "free", space = "free_x") +
  facet_grid(variable ~ disease_assay, scale = "free", space = "free_x") +
  scale_y_log10(labels = f2si) +
  scale_fill_manual(values = meta_colors$type) +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Mass Cytometry")

# Low-input RNA-seq
# -----------------------------------------------------------------------------
md_lirna <- melt(
  data = dd,
  id.vars = c("subject_id", "disease_assay"),
  measure.vars = c(
    "low_input_fibroblasts", "low_input_monoytes",
    "low_input_t_cells", "low_input_b_cells"
  )
)
md_lirna$variable <- c(
  "low_input_fibroblasts" = "Fibroblast",
  "low_input_monoytes" = "Monocyte",
  "low_input_t_cells" = "T cell",
  "low_input_b_cells" = "B cell"
)[md_lirna$variable]
md_lirna$variable <- factor(
  md_lirna$variable, c("B cell", "T cell", "Fibroblast", "Monocyte")
)
p_lirna <- ggplot(
  md_lirna, aes(x = subject_id, weight = value, fill = variable)
) +
  geom_bar() +
  facet_grid(variable ~ disease_assay, scale = "free", space = "free_x") +
  scale_y_log10(labels = f2si) +
  scale_fill_manual(values = meta_colors$type) +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Low-input RNA-seq")

# Single-cell RNA-seq
# -----------------------------------------------------------------------------
md_scrna <- melt(
  data = dd,
  id.vars = c("subject_id", "disease_assay"),
  measure.vars = c(
    "t_cells_for_single_cell",
    "b_cells_for_single_cell",
    "monoytes_for_single_cell",
    "fibroblasts_for_single_cell"
  )
)
md_scrna$variable <- c(
  "t_cells_for_single_cell" = "T cell",
  "b_cells_for_single_cell" = "B cell",
  "monoytes_for_single_cell" = "Monocyte",
  "fibroblasts_for_single_cell" = "Fibroblast"
)[md_scrna$variable]
md_scrna$variable <- factor(
  md_scrna$variable, c("B cell", "T cell", "Fibroblast", "Monocyte")
)
p_scrna <- ggplot(
  md_scrna, aes(x = subject_id, weight = value, fill = variable)
) +
  geom_bar() +
  facet_grid(variable ~ disease_assay, scale = "free", space = "free_x") +
  # scale_y_log10(labels = f2si) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  scale_fill_manual(values = meta_colors$type, name = NULL) +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = "Single-cell RNA-seq")
# p_scrna

# Joint, Sex
# -----------------------------------------------------------------------------
md_meta <- melt(
  data = dd,
  id.vars = c("subject_id", "disease_assay"),
  measure.vars = c(
    "joint", "sex"
  )
)
p_meta <- ggplot() +
  geom_tile(
    data = md_meta,
    mapping = aes(x = subject_id, y = variable, fill = value)
  ) +
  scale_fill_manual(values = c(meta_colors$joint, meta_colors$sex)) +
  theme_bw(base_size = 20) +
  facet_grid(~ disease_assay, scale = "free", space = "free_x") +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    # strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Metadata")

# Combined plot
# -----------------------------------------------------------------------------
# p_all <- rbind(
#   ggplotGrob(p_flow),
#   ggplotGrob(p_lirna),
#   ggplotGrob(p_cytof),
#   ggplotGrob(p_scrna)
# )
# pdf("data-overview.pdf", width = 7, height = 12)
# grid::grid.newpage()
# grid::grid.draw(p_all)
# dev.off()

pdf("figures/Phase1/data-overview-detailed.pdf", width = 18, height = 15)
# png("data-overview.png", width = 14, height = 13)
plot_grid(
  p_flow,
  p_lirna,
  p_cytof,
  p_scrna,
  p_meta,
  align = "v",
  axis = "l",
  nrow = 5,
  rel_heights = c(2.5, 2, 2, 3.5, 1)
)
dev.off()

# Simple version
# -----------------------------------------------------------------------------

md_flow$fill <- md_flow$variable
md_flow$fill[md_flow$value == 0 | is.na(md_flow$value)] <- NA
p_flow <- ggplot() +
  geom_tile(
    data = md_flow,
    mapping = aes(x = subject_id, y = variable, fill = fill)
  ) +
  scale_fill_manual(values = meta_colors$type) +
  facet_grid(~ disease_assay, scale = "free", space = "free_x") +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    # strip.background = element_blank(),
    # strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Flow Cytometry")

md_lirna$fill <- md_lirna$variable
md_lirna$fill[md_lirna$value == 0 | is.na(md_lirna$value)] <- NA
p_lirna <- ggplot() +
  geom_tile(
    data = md_lirna,
    mapping = aes(x = subject_id, y = variable, fill = fill)
  ) +
  scale_fill_manual(values = meta_colors$type) +
  facet_grid(~ disease_assay, scale = "free", space = "free_x") +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    # strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Low-input RNA-seq")

md_cytof$fill <- md_cytof$variable
md_cytof$fill[md_cytof$value == 0 | is.na(md_cytof$value)] <- NA
p_cytof <- ggplot() +
  geom_tile(
    data = md_cytof,
    mapping = aes(x = subject_id, y = variable, fill = fill)
  ) +
  scale_fill_manual(values = meta_colors$type) +
  facet_grid(~ disease_assay, scale = "free", space = "free_x") +
  theme_bw(base_size = 20) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    # strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Mass Cytometry")

md_scrna$fill <- md_scrna$variable
md_scrna$fill[md_scrna$value == 0 | is.na(md_scrna$value)] <- NA
p_scrna <- ggplot() +
  geom_tile(
    data = md_scrna,
    mapping = aes(x = subject_id, y = variable, fill = fill)
  ) +
  scale_fill_manual(values = meta_colors$type) +
  theme_bw(base_size = 20) +
  facet_grid(~ disease_assay, scale = "free", space = "free_x") +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 90),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    # strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Single-cell RNA-seq")


pdf("figures/Phase1/data-overview-simple.pdf", width = 8, height = 8)
plot_grid(
  p_flow,
  p_lirna,
  p_cytof,
  p_scrna,
  p_meta,
  align = "v",
  axis = "l",
  nrow = 5,
  rel_heights = c(3, 2, 2, 2, 1.5)
)
dev.off()

ggplot2::ggsave(
  plot = p_meta + theme(legend.position = "right"),
  device = "pdf",
  filename = "figures/Phase1/data-overview-meta.pdf",
  height = 5
)
