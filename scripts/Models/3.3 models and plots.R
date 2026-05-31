############################################################
## 3.3.2 Cross-colony comparison of NN1-NN8 spacing
## Absolute spacing, clustering, and normalized scaling
############################################################

## ---------------------------------------------------------
## 0. Packages and paths
## ---------------------------------------------------------

library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)

setwd("/Users/amalia/Documents/GitHub/pattern_formation")


input_dir  <- "inputs"
output_dir <- "outputs"
plot_dir   <- "plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------
## 1. Read data
## ---------------------------------------------------------

nn_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN1to8_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

nn_summary_raw <- read_csv(
  file.path(input_dir, "Colony_NN1to8_pooled_mean_sd_se_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

cat("\n--- Absolute NN dataset columns ---\n")
print(names(nn_raw))

cat("\n--- Summary dataset columns ---\n")
print(names(nn_summary_raw))

cat("\n--- Ratio dataset columns ---\n")
print(names(nn_ratios_raw))

## ---------------------------------------------------------
## 2. Prepare absolute NN1-NN8 long dataset
## ---------------------------------------------------------

geno_levels <- c("1", "2", "3", "4", "6", "7", "8", "9", "10")

nn_long <- nn_raw %>%
  mutate(
    genotype = case_when(
      "colony_number" %in% names(.) ~ as.character(colony_number),
      "colony" %in% names(.) ~ str_remove(as.character(colony), "^SC"),
      TRUE ~ NA_character_
    ),
    genotype = factor(genotype, levels = geno_levels),
    slide = case_when(
      "slide" %in% names(.) ~ as.character(slide),
      "sample" %in% names(.) ~ as.character(sample),
      TRUE ~ "unknown_slide"
    ),
    polyp_id = case_when(
      "polyp_id" %in% names(.) ~ as.character(polyp_id),
      TRUE ~ as.character(row_number())
    ),
    nubbin_id = interaction(genotype, slide, drop = TRUE)
  ) %>%
  pivot_longer(
    cols = matches("^nn[1-8]_distance$"),
    names_to = "nn_rank",
    values_to = "distance"
  ) %>%
  mutate(
    nn_rank_num = as.numeric(str_extract(nn_rank, "[1-8]")),
    nn_rank = factor(
      paste0("NN", nn_rank_num),
      levels = paste0("NN", 1:8)
    )
  ) %>%
  filter(
    !is.na(genotype),
    !is.na(distance),
    is.finite(distance),
    distance > 0
  )

cat("\n--- Absolute NN long data ---\n")
print(glimpse(nn_long))

cat("\n--- Counts by genotype and NN rank ---\n")
print(table(nn_long$genotype, nn_long$nn_rank))

## ---------------------------------------------------------
## 3. Model 1: colony differences without assuming groups
## ---------------------------------------------------------
## This tests whether colonies differ in absolute spacing and
## whether distance increases with NN rank.

m_abs_genotype <- lmer(
  distance ~ genotype * nn_rank_num + (1 | nubbin_id),
  data = nn_long
)

cat("\n==============================\n")
cat("ABSOLUTE SPACING MODEL: genotype * NN rank\n")
cat("==============================\n")
print(anova(m_abs_genotype))
print(summary(m_abs_genotype))

abs_genotype_emm <- emmeans(m_abs_genotype, ~ genotype)
abs_genotype_slopes <- emtrends(m_abs_genotype, ~ genotype, var = "nn_rank_num")

cat("\n--- Estimated mean spacing by genotype ---\n")
print(abs_genotype_emm)

cat("\n--- Estimated NN-rank slopes by genotype ---\n")
print(abs_genotype_slopes)

write_csv(
  as.data.frame(anova(m_abs_genotype)) %>% rownames_to_column("term"),
  file.path(output_dir, "33_abs_spacing_genotype_anova.csv")
)

write_csv(
  as.data.frame(abs_genotype_emm),
  file.path(output_dir, "33_abs_spacing_genotype_emmeans.csv")
)

write_csv(
  as.data.frame(abs_genotype_slopes),
  file.path(output_dir, "33_abs_spacing_genotype_slopes.csv")
)

## ---------------------------------------------------------
## 4. Cluster colonies by NN1-NN8 absolute spacing profile
## ---------------------------------------------------------
## This tests whether the apparent two groups emerge from the data.

profile_df <- nn_long %>%
  group_by(genotype, nn_rank_num) %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  )

profile_mat <- profile_df %>%
  pivot_wider(
    names_from = nn_rank_num,
    values_from = mean_distance,
    names_prefix = "NN"
  ) %>%
  column_to_rownames("genotype") %>%
  as.matrix()

cat("\n==============================\n")
cat("COLONY MEAN PROFILE MATRIX\n")
cat("==============================\n")
print(profile_mat)

hc_abs <- hclust(dist(profile_mat), method = "ward.D2")
cluster_abs <- cutree(hc_abs, k = 2)

cluster_df <- tibble(
  genotype = names(cluster_abs),
  spacing_cluster_raw = paste0("Cluster_", cluster_abs)
) %>%
  mutate(
    proposed_group = case_when(
      genotype %in% c("2", "3", "4", "6") ~ "High-spacing",
      genotype %in% c("1", "7", "8", "9", "10") ~ "Low-spacing",
      TRUE ~ NA_character_
    )
  )

## Give clusters interpretable labels based on mean distance
cluster_means <- nn_long %>%
  left_join(cluster_df, by = "genotype") %>%
  group_by(spacing_cluster_raw) %>%
  summarise(
    cluster_mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(cluster_mean_distance) %>%
  mutate(spacing_cluster = c("Low-spacing", "High-spacing"))

cluster_df <- cluster_df %>%
  left_join(cluster_means, by = "spacing_cluster_raw") %>%
  mutate(
    spacing_cluster = factor(spacing_cluster, levels = c("Low-spacing", "High-spacing")),
    proposed_group = factor(proposed_group, levels = c("Low-spacing", "High-spacing"))
  )

cat("\n==============================\n")
cat("CLUSTER ASSIGNMENTS\n")
cat("==============================\n")
print(cluster_df)

cat("\n--- Cluster vs David proposed grouping ---\n")
print(table(cluster_df$spacing_cluster, cluster_df$proposed_group))

write_csv(
  cluster_df,
  file.path(output_dir, "33_spacing_cluster_assignments.csv")
)

## Optional dendrogram
# ---------------------------------------------------------
# Preview dendrogram in R plotting pane
# ---------------------------------------------------------

# Save current graphics settings
old_par <- par(no.readonly = TRUE)

# Set margins: bottom, left, top, right
par(mar = c(7, 5, 5, 2))

plot(
  hc_abs,
  main = "Clustering of colonies by NN1–NN8 absolute spacing profiles",
  xlab = "",
  sub = "",
  cex.main = 1.2,
  cex.axis = 1,
  hang = -1
)

rect.hclust(
  hc_abs,
  k = 2,
  border = "red"
)

mtext(
  "Colony",
  side = 1,
  line = 5,
  cex = 1.1
)
# =========================================================
# Dendrogram with colored colony labels instead of rectangles
# =========================================================

library(dendextend)
library(dplyr)

# Make dendrogram object
dend <- as.dendrogram(hc_abs)

# Get cluster assignment
cluster_lookup <- cluster_df %>%
  mutate(
    genotype = as.character(genotype),
    spacing_cluster = as.character(spacing_cluster)
  ) %>%
  dplyr::select(genotype, spacing_cluster)

# Color labels by cluster
label_cols <- ifelse(
  labels(dend) %in% cluster_lookup$genotype[cluster_lookup$spacing_cluster == "High-spacing"],
  "#CC79A7",   # high-spacing
  "#56B4E9"    # low-spacing
)

labels_colors(dend) <- label_cols
labels_cex(dend) <- 1.1

# Color branches into two clusters
dend <- color_branches(
  dend,
  k = 2,
  col = c("#56B4E9", "#CC79A7")
)

# Make branches thicker
dend <- dendextend::set(dend, "branches_lwd", 2.5)
# ---------------------------------------------------------
# Preview
# ---------------------------------------------------------
old_par <- par(no.readonly = TRUE)

par(mar = c(6, 5, 5, 2))

plot(
  dend,
  main = "Clustering of colonies by NN1–NN8 absolute spacing profiles",
  ylab = "Height",
  cex.main = 1.2
)

mtext("Colony", side = 1, line = 4.5, cex = 1.1)

legend(
  "topleft",
  legend = c("Low-spacing", "High-spacing"),
  col = c("#56B4E9", "#CC79A7"),
  lwd = 2.5,
  bty = "n",
  cex = 0.9
)

par(old_par)

# ---------------------------------------------------------
# Save PNG
# ---------------------------------------------------------
png(
  file.path(plot_dir, "33_absolute_spacing_cluster_dendrogram_colored.png"),
  width = 2000,
  height = 1600,
  res = 250
)

par(mar = c(6, 5, 5, 2))

plot(
  dend,
  main = "Clustering of colonies by NN1–NN8 absolute spacing profiles",
  ylab = "Height",
  cex.main = 1.2
)

mtext("Colony", side = 1, line = 4.5, cex = 1.1)

legend(
  "topright",
  legend = c("Low-spacing", "High-spacing"),
  col = c("#56B4E9", "#CC79A7"),
  lwd = 2.5,
  bty = "n",
  cex = 0.9
)

dev.off()

# ---------------------------------------------------------
# Save PDF
# ---------------------------------------------------------
pdf(
  file.path(plot_dir, "33_absolute_spacing_cluster_dendrogram_colored.pdf"),
  width = 8,
  height = 6
)

par(mar = c(6, 5, 5, 2))

plot(
  dend,
  main = "Clustering of colonies by NN1–NN8 absolute spacing profiles",
  ylab = "Height",
  cex.main = 1.2
)

mtext("Colony", side = 1, line = 4.5, cex = 1.1)

legend(
  "topright",
  legend = c("Low-spacing", "High-spacing"),
  col = c("#56B4E9", "#CC79A7"),
  lwd = 2.5,
  bty = "n",
  cex = 0.9
)

dev.off()

## ---------------------------------------------------------
## 5. Model 2: test data-derived high/low spacing clusters
## ---------------------------------------------------------

nn_long_cluster <- nn_long %>%
  left_join(cluster_df %>% select(genotype, spacing_cluster), by = "genotype") %>%
  mutate(
    spacing_cluster = factor(spacing_cluster, levels = c("Low-spacing", "High-spacing"))
  )

m_abs_cluster <- lmer(
  distance ~ spacing_cluster * nn_rank_num + (1 | nubbin_id),
  data = nn_long_cluster
)

cat("\n==============================\n")
cat("ABSOLUTE SPACING MODEL: cluster * NN rank\n")
cat("==============================\n")
print(anova(m_abs_cluster))
print(summary(m_abs_cluster))

abs_cluster_emm <- emmeans(m_abs_cluster, ~ spacing_cluster)
abs_cluster_slopes <- emtrends(m_abs_cluster, ~ spacing_cluster, var = "nn_rank_num")

cat("\n--- Estimated mean spacing by cluster ---\n")
print(abs_cluster_emm)
print(pairs(abs_cluster_emm))

cat("\n--- Estimated NN-rank slopes by cluster ---\n")
print(abs_cluster_slopes)
print(pairs(abs_cluster_slopes))

write_csv(
  as.data.frame(anova(m_abs_cluster)) %>% rownames_to_column("term"),
  file.path(output_dir, "33_abs_spacing_cluster_anova.csv")
)

write_csv(
  as.data.frame(abs_cluster_emm),
  file.path(output_dir, "33_abs_spacing_cluster_emmeans.csv")
)

write_csv(
  as.data.frame(abs_cluster_slopes),
  file.path(output_dir, "33_abs_spacing_cluster_slopes.csv")
)

## ---------------------------------------------------------
## 6. Prepare normalized NNk / NN1 dataset
## ---------------------------------------------------------

ratio_long <- nn_ratios_raw %>%
  mutate(
    genotype = case_when(
      "colony_number" %in% names(.) ~ as.character(colony_number),
      "colony" %in% names(.) ~ str_remove(as.character(colony), "^SC"),
      TRUE ~ NA_character_
    ),
    genotype = factor(genotype, levels = geno_levels),
    slide = case_when(
      "slide" %in% names(.) ~ as.character(slide),
      "sample" %in% names(.) ~ as.character(sample),
      TRUE ~ "unknown_slide"
    ),
    polyp_id = case_when(
      "polyp_id" %in% names(.) ~ as.character(polyp_id),
      TRUE ~ as.character(row_number())
    ),
    nubbin_id = interaction(genotype, slide, drop = TRUE)
  )

## Flexible handling depending on whether ratios are already long or wide
if ("ratio_value" %in% names(ratio_long)) {
  
  ratio_long <- ratio_long %>%
    mutate(
      ratio_name = case_when(
        "ratio_name" %in% names(.) ~ as.character(ratio_name),
        TRUE ~ NA_character_
      ),
      nn_rank_num = as.numeric(str_extract(ratio_name, "(?<=NN)\\d+|\\d+(?=/NN1)")),
      nn_rank = factor(paste0("NN", nn_rank_num), levels = paste0("NN", 2:8)),
      ratio_to_NN1 = ratio_value
    ) %>%
    filter(
      !is.na(genotype),
      !is.na(nn_rank_num),
      nn_rank_num >= 2,
      nn_rank_num <= 8,
      !is.na(ratio_to_NN1),
      is.finite(ratio_to_NN1),
      ratio_to_NN1 > 0
    )
  
} else {
  
  ratio_long <- ratio_long %>%
    pivot_longer(
      cols = matches("^nn[2-8].*nn1|^nn[2-8]_ratio$|^nn[2-8]_to_nn1$"),
      names_to = "ratio_name",
      values_to = "ratio_to_NN1"
    ) %>%
    mutate(
      nn_rank_num = as.numeric(str_extract(ratio_name, "[2-8]")),
      nn_rank = factor(paste0("NN", nn_rank_num), levels = paste0("NN", 2:8))
    ) %>%
    filter(
      !is.na(genotype),
      !is.na(nn_rank_num),
      !is.na(ratio_to_NN1),
      is.finite(ratio_to_NN1),
      ratio_to_NN1 > 0
    )
}

ratio_long_cluster <- ratio_long %>%
  left_join(cluster_df %>% select(genotype, spacing_cluster), by = "genotype") %>%
  mutate(
    spacing_cluster = factor(spacing_cluster, levels = c("Low-spacing", "High-spacing"))
  )

cat("\n--- Normalized ratio long data ---\n")
print(glimpse(ratio_long_cluster))

cat("\n--- Counts by genotype and NN rank for ratios ---\n")
print(table(ratio_long_cluster$genotype, ratio_long_cluster$nn_rank))

## ---------------------------------------------------------
## 7. Model 3: normalized ratios by genotype
## ---------------------------------------------------------

m_norm_genotype <- lmer(
  ratio_to_NN1 ~ genotype * nn_rank_num + (1 | nubbin_id),
  data = ratio_long_cluster
)

cat("\n==============================\n")
cat("NORMALIZED SPACING MODEL: genotype * NN rank\n")
cat("==============================\n")
print(anova(m_norm_genotype))
print(summary(m_norm_genotype))

norm_genotype_emm <- emmeans(m_norm_genotype, ~ genotype)
norm_genotype_slopes <- emtrends(m_norm_genotype, ~ genotype, var = "nn_rank_num")

cat("\n--- Normalized ratio means by genotype ---\n")
print(norm_genotype_emm)

cat("\n--- Normalized ratio slopes by genotype ---\n")
print(norm_genotype_slopes)

write_csv(
  as.data.frame(anova(m_norm_genotype)) %>% rownames_to_column("term"),
  file.path(output_dir, "33_norm_spacing_genotype_anova.csv")
)

write_csv(
  as.data.frame(norm_genotype_emm),
  file.path(output_dir, "33_norm_spacing_genotype_emmeans.csv")
)

write_csv(
  as.data.frame(norm_genotype_slopes),
  file.path(output_dir, "33_norm_spacing_genotype_slopes.csv")
)

## ---------------------------------------------------------
## 8. Model 4: do high/low clusters remain after normalization?
## ---------------------------------------------------------

m_norm_cluster <- lmer(
  ratio_to_NN1 ~ spacing_cluster * nn_rank_num + (1 | nubbin_id),
  data = ratio_long_cluster
)

cat("\n==============================\n")
cat("NORMALIZED SPACING MODEL: cluster * NN rank\n")
cat("==============================\n")
print(anova(m_norm_cluster))
print(summary(m_norm_cluster))

norm_cluster_emm <- emmeans(m_norm_cluster, ~ spacing_cluster)
norm_cluster_slopes <- emtrends(m_norm_cluster, ~ spacing_cluster, var = "nn_rank_num")

cat("\n--- Normalized ratio means by cluster ---\n")
print(norm_cluster_emm)
print(pairs(norm_cluster_emm))

cat("\n--- Normalized ratio slopes by cluster ---\n")
print(norm_cluster_slopes)
print(pairs(norm_cluster_slopes))

write_csv(
  as.data.frame(anova(m_norm_cluster)) %>% rownames_to_column("term"),
  file.path(output_dir, "33_norm_spacing_cluster_anova.csv")
)

write_csv(
  as.data.frame(norm_cluster_emm),
  file.path(output_dir, "33_norm_spacing_cluster_emmeans.csv")
)

write_csv(
  as.data.frame(norm_cluster_slopes),
  file.path(output_dir, "33_norm_spacing_cluster_slopes.csv")
)

## ---------------------------------------------------------
## 9. Descriptive summaries for text
## ---------------------------------------------------------

abs_summary <- nn_long_cluster %>%
  group_by(spacing_cluster, nn_rank) %>%
  summarise(
    n = n(),
    mean_distance = mean(distance, na.rm = TRUE),
    sd_distance = sd(distance, na.rm = TRUE),
    se_distance = sd_distance / sqrt(n),
    .groups = "drop"
  )

norm_summary <- ratio_long_cluster %>%
  group_by(spacing_cluster, nn_rank) %>%
  summarise(
    n = n(),
    mean_ratio = mean(ratio_to_NN1, na.rm = TRUE),
    sd_ratio = sd(ratio_to_NN1, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    .groups = "drop"
  )

write_csv(abs_summary, file.path(output_dir, "33_abs_spacing_summary_by_cluster_rank.csv"))
write_csv(norm_summary, file.path(output_dir, "33_norm_spacing_summary_by_cluster_rank.csv"))

cat("\n==============================\n")
cat("ABSOLUTE SUMMARY BY CLUSTER AND RANK\n")
cat("==============================\n")
print(abs_summary)

cat("\n==============================\n")
cat("NORMALIZED SUMMARY BY CLUSTER AND RANK\n")
cat("==============================\n")
print(norm_summary)

## ---------------------------------------------------------
## 10. Basic figures for checking
## ---------------------------------------------------------

p_abs_box <- ggplot(nn_long_cluster, aes(x = genotype, y = distance, fill = spacing_cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, color = "black") +
  geom_jitter(width = 0.12, alpha = 0.12, size = 0.6) +
  facet_wrap(~ nn_rank, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("Low-spacing" = "#56B4E9", "High-spacing" = "#CC79A7")) +
  labs(
    x = "Colony",
    y = "Nearest-neighbor distance (mm)",
    fill = "Spacing group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

p_abs_scaling <- ggplot(abs_summary,
                        aes(x = as.numeric(str_remove(nn_rank, "NN")),
                            y = mean_distance,
                            color = spacing_cluster,
                            group = spacing_cluster)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = mean_distance - se_distance,
        ymax = mean_distance + se_distance),
    width = 0.12,
    linewidth = 0.5
  ) +
  scale_color_manual(values = c("Low-spacing" = "#56B4E9", "High-spacing" = "#CC79A7")) +
  scale_x_continuous(breaks = 1:8) +
  labs(
    x = "Nearest-neighbor rank",
    y = "Mean distance (mm)",
    color = "Spacing group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

p_norm_scaling <- ggplot(norm_summary,
                         aes(x = as.numeric(str_remove(nn_rank, "NN")),
                             y = mean_ratio,
                             color = spacing_cluster,
                             group = spacing_cluster)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = mean_ratio - se_ratio,
        ymax = mean_ratio + se_ratio),
    width = 0.12,
    linewidth = 0.5
  ) +
  scale_color_manual(values = c("Low-spacing" = "#56B4E9", "High-spacing" = "#CC79A7")) +
  scale_x_continuous(breaks = 2:8) +
  labs(
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1]),
    color = "Spacing group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(p_abs_box)
print(p_abs_scaling)
print(p_norm_scaling)

ggsave(
  file.path(plot_dir, "33_absolute_spacing_boxplots_by_rank.png"),
  p_abs_box,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "33_absolute_spacing_scaling_by_cluster.png"),
  p_abs_scaling,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "33_normalized_spacing_scaling_by_cluster.png"),
  p_norm_scaling,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

pdf(file.path(plot_dir, "33_absolute_spacing_cluster_dendrogram.pdf"),
    width = 7, height = 5.5)

############################################################
## End 3.3.2 code
############################################################

# =========================================================
# Figure X. Absolute and proportional scaling of polyp spacing
# Panels labelled as (a), (b), (c), (d)
# =========================================================

library(tidyverse)
library(patchwork)
library(grid)

# ---------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------
# Assumes these already exist:
# output_dir <- "outputs"
# plot_dir   <- "plots"

# ---------------------------------------------------------
# 1. Read model summaries and descriptive summaries
# ---------------------------------------------------------

abs_summary <- read_csv(
  file.path(output_dir, "33_abs_spacing_summary_by_cluster_rank.csv"),
  show_col_types = FALSE
)

norm_summary <- read_csv(
  file.path(output_dir, "33_norm_spacing_summary_by_cluster_rank.csv"),
  show_col_types = FALSE
)

abs_cluster_emm <- read_csv(
  file.path(output_dir, "33_abs_spacing_cluster_emmeans.csv"),
  show_col_types = FALSE
)

abs_cluster_slopes <- read_csv(
  file.path(output_dir, "33_abs_spacing_cluster_slopes.csv"),
  show_col_types = FALSE
)

norm_cluster_emm <- read_csv(
  file.path(output_dir, "33_norm_spacing_cluster_emmeans.csv"),
  show_col_types = FALSE
)

norm_cluster_slopes <- read_csv(
  file.path(output_dir, "33_norm_spacing_cluster_slopes.csv"),
  show_col_types = FALSE
)

# ---------------------------------------------------------
# 2. Shared settings
# ---------------------------------------------------------

group_cols <- c(
  "Low-spacing"  = "#56B4E9",
  "High-spacing" = "#CC79A7"
)

theme_set(
  theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text  = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 18, face = "bold"),
      strip.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.key.size = unit(1.1, "lines"),
      panel.grid.minor = element_blank()
    )
)

# ---------------------------------------------------------
# 3. Prepare data
# ---------------------------------------------------------

abs_summary <- abs_summary %>%
  mutate(
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    ),
    nn_rank_num = as.numeric(str_remove(as.character(nn_rank), "NN"))
  )

norm_summary <- norm_summary %>%
  mutate(
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    ),
    nn_rank_num = as.numeric(str_remove(as.character(nn_rank), "NN"))
  )

# ---------------------------------------------------------
# 4. Panel (a): absolute spacing lines
# ---------------------------------------------------------

pA <- ggplot(
  abs_summary,
  aes(
    x = nn_rank_num,
    y = mean_distance,
    color = spacing_cluster,
    group = spacing_cluster
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_distance - se_distance,
      ymax = mean_distance + se_distance
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  scale_color_manual(values = group_cols) +
  scale_x_continuous(breaks = 1:8) +
  labs(
    x = "Nearest-neighbor rank",
    y = "Mean distance (mm)",
    color = "Spacing group",
    title = "(a) Absolute spacing"
  ) +
  theme(
    legend.position = "top"
  )

# ---------------------------------------------------------
# 5. Panel (b): absolute model estimates
# ---------------------------------------------------------

abs_estimates <- bind_rows(
  abs_cluster_emm %>%
    transmute(
      spacing_cluster,
      metric = "Mean distance",
      estimate = emmean,
      SE = SE,
      lower = asymp.LCL,
      upper = asymp.UCL
    ),
  abs_cluster_slopes %>%
    transmute(
      spacing_cluster,
      metric = "Scaling slope",
      estimate = nn_rank_num.trend,
      SE = SE,
      lower = asymp.LCL,
      upper = asymp.UCL
    )
) %>%
  mutate(
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    ),
    metric = factor(
      metric,
      levels = c("Mean distance", "Scaling slope")
    )
  )

pB <- ggplot(
  abs_estimates,
  aes(
    x = spacing_cluster,
    y = estimate,
    color = spacing_cluster
  )
) +
  geom_point(size = 3.5) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.7
  ) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = group_cols) +
  labs(
    x = NULL,
    y = "Model estimate",
    color = "Spacing group",
    title = "(b) Absolute model estimates"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# ---------------------------------------------------------
# 6. Panel (c): normalized ratio lines
# ---------------------------------------------------------

pC <- ggplot(
  norm_summary,
  aes(
    x = nn_rank_num,
    y = mean_ratio,
    color = spacing_cluster,
    group = spacing_cluster
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_ratio - se_ratio,
      ymax = mean_ratio + se_ratio
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  scale_color_manual(values = group_cols) +
  scale_x_continuous(breaks = 2:8) +
  labs(
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1]),
    color = "Spacing group",
    title = "(c) Proportional spacing"
  ) +
  theme(
    legend.position = "top"
  )

# ---------------------------------------------------------
# 7. Panel (d): normalized model estimates
# ---------------------------------------------------------

norm_estimates <- bind_rows(
  norm_cluster_emm %>%
    transmute(
      spacing_cluster,
      metric = "mean_ratio",
      estimate = emmean,
      SE = SE,
      lower = asymp.LCL,
      upper = asymp.UCL
    ),
  norm_cluster_slopes %>%
    transmute(
      spacing_cluster,
      metric = "slope",
      estimate = nn_rank_num.trend,
      SE = SE,
      lower = asymp.LCL,
      upper = asymp.UCL
    )
) %>%
  mutate(
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    ),
    metric = factor(
      metric,
      levels = c("mean_ratio", "slope")
    )
  )

facet_labels_d <- as_labeller(
  c(
    mean_ratio = "bold(Mean)~bold(NN)[k]/bold(NN)[1]",
    slope = "bold(Slope)"
  ),
  label_parsed
)

pD <- ggplot(
  norm_estimates,
  aes(
    x = spacing_cluster,
    y = estimate,
    color = spacing_cluster
  )
) +
  geom_point(size = 3.5) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.12,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1,
    labeller = facet_labels_d
  ) +
  scale_color_manual(values = group_cols) +
  labs(
    x = NULL,
    y = "Model estimate",
    color = "Spacing group",
    title = "(d) Normalized model estimates"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.text = element_text(size = 14, face = "bold")
  )
# ---------------------------------------------------------
# 8. Composite figure
# ---------------------------------------------------------

figure_spacing_scaling <- (pA | pB) / (pC | pD) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top"
  )

print(figure_spacing_scaling)

# ---------------------------------------------------------
# 9. Save composite
# ---------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "33_absolute_proportional_spacing_composite_panel_letters.png"),
  plot = figure_spacing_scaling,
  width = 12,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_absolute_proportional_spacing_composite_panel_letters.pdf"),
  plot = figure_spacing_scaling,
  width = 12,
  height = 9,
  bg = "white"
)

# ---------------------------------------------------------
# 10. Save individual panels
# ---------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "33_panel_a_absolute_spacing_lines.png"),
  plot = pA,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_a_absolute_spacing_lines.pdf"),
  plot = pA,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_b_absolute_model_estimates.png"),
  plot = pB,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_b_absolute_model_estimates.pdf"),
  plot = pB,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_c_normalized_spacing_lines.png"),
  plot = pC,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_c_normalized_spacing_lines.pdf"),
  plot = pC,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_d_normalized_model_estimates.png"),
  plot = pD,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "33_panel_d_normalized_model_estimates.pdf"),
  plot = pD,
  width = 7,
  height = 5,
  bg = "white"
)

# Average NNk/NN1 across all colonies, with SE
library(tidyverse)
library(janitor)

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_dir  <- "inputs"
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

avg_ratio_all_colonies <- nn_ratios_raw %>%
  mutate(
    ratio_name = as.character(ratio_name),
    nn_rank_num = as.numeric(str_extract(ratio_name, "(?<=NN)\\d+|\\d+(?=/NN1)")),
    nn_rank = paste0("NN", nn_rank_num),
    ratio_to_NN1 = ratio_value
  ) %>%
  filter(
    !is.na(nn_rank_num),
    nn_rank_num >= 2,
    nn_rank_num <= 8,
    !is.na(ratio_to_NN1),
    is.finite(ratio_to_NN1),
    ratio_to_NN1 > 0
  ) %>%
  group_by(nn_rank, nn_rank_num) %>%
  summarise(
    n = n(),
    mean_NNk_over_NN1 = mean(ratio_to_NN1, na.rm = TRUE),
    sd = sd(ratio_to_NN1, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(nn_rank_num)

print(avg_ratio_all_colonies)

write_csv(
  avg_ratio_all_colonies,
  file.path(output_dir, "average_NNk_over_NN1_all_colonies_with_SE.csv")
)

avg_ratio_all_colonies_pretty <- avg_ratio_all_colonies %>%
  mutate(
    mean_SE = sprintf("%.3f ± %.3f", mean_NNk_over_NN1, se)
  ) %>%
  select(nn_rank, n, mean_SE)

print(avg_ratio_all_colonies_pretty)

write_csv(
  avg_ratio_all_colonies_pretty,
  file.path(output_dir, "average_NNk_over_NN1_all_colonies_pretty.csv")
)



# Average NNk/NN1 across all colonies, with SE
library(tidyverse)
library(janitor)

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_dir  <- "inputs"
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

avg_ratio_all_colonies <- nn_ratios_raw %>%
  mutate(
    colony = as.character(colony),
    genotype = as.character(colony_number),
    numerator = as.character(numerator),
    denominator = as.character(denominator),
    nn_rank_num = as.numeric(str_remove(numerator, "NN")),
    nn_rank = numerator,
    ratio_to_NN1 = ratio_value
  ) %>%
  filter(
    denominator == "NN1",
    numerator %in% paste0("NN", 2:8),
    !is.na(ratio_to_NN1),
    is.finite(ratio_to_NN1),
    ratio_to_NN1 > 0
  ) %>%
  group_by(nn_rank, nn_rank_num) %>%
  summarise(
    n = n(),
    mean_NNk_over_NN1 = mean(ratio_to_NN1, na.rm = TRUE),
    sd = sd(ratio_to_NN1, na.rm = TRUE),
    se = sd / sqrt(n),
    mean_SE = sprintf("%.3f ± %.3f", mean_NNk_over_NN1, se),
    .groups = "drop"
  ) %>%
  arrange(nn_rank_num)

print(avg_ratio_all_colonies)

write_csv(
  avg_ratio_all_colonies,
  file.path(output_dir, "average_NNk_over_NN1_all_colonies_with_SE_CORRECTED.csv")
)

# Average NNk/NN1 by colony, with SE

avg_ratio_by_colony <- nn_ratios_raw %>%
  mutate(
    colony = as.character(colony),
    genotype = as.character(colony_number),
    numerator = as.character(numerator),
    denominator = as.character(denominator),
    nn_rank_num = as.numeric(str_remove(numerator, "NN")),
    nn_rank = numerator,
    ratio_to_NN1 = ratio_value
  ) %>%
  filter(
    denominator == "NN1",
    numerator %in% paste0("NN", 2:8),
    !is.na(ratio_to_NN1),
    is.finite(ratio_to_NN1),
    ratio_to_NN1 > 0
  ) %>%
  group_by(genotype, colony, nn_rank, nn_rank_num) %>%
  summarise(
    n = n(),
    mean_NNk_over_NN1 = mean(ratio_to_NN1, na.rm = TRUE),
    sd = sd(ratio_to_NN1, na.rm = TRUE),
    se = sd / sqrt(n),
    mean_SE = sprintf("%.3f ± %.3f", mean_NNk_over_NN1, se),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(genotype), nn_rank_num)

print(avg_ratio_by_colony)

write_csv(
  avg_ratio_by_colony,
  file.path(output_dir, "average_NNk_over_NN1_by_colony_with_SE_CORRECTED.csv")
)

avg_ratio_by_colony_wide <- avg_ratio_by_colony %>%
  select(colony, nn_rank, mean_SE) %>%
  pivot_wider(
    names_from = nn_rank,
    values_from = mean_SE
  ) %>%
  arrange(as.numeric(str_remove(colony, "^SC")))

print(avg_ratio_by_colony_wide)

write_csv(
  avg_ratio_by_colony_wide,
  file.path(output_dir, "average_NNk_over_NN1_by_colony_wide_CORRECTED.csv")
)

avg_ratio_all_colonies %>%
  arrange(nn_rank_num) %>%
  mutate(
    step_difference = mean_NNk_over_NN1 - lag(mean_NNk_over_NN1)
  )

lm_arithmetic <- lm(
  mean_NNk_over_NN1 ~ nn_rank_num,
  data = avg_ratio_all_colonies
)

summary(lm_arithmetic)

coef(lm_arithmetic)


ggplot(avg_ratio_all_colonies, aes(x = nn_rank_num, y = mean_NNk_over_NN1)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_NNk_over_NN1 - se,
      ymax = mean_NNk_over_NN1 + se
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = 2:8,
    labels = paste0("NN", 2:8)
  ) +
  labs(
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1])
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    panel.grid.minor = element_blank()
  )

library(tidyverse)
library(janitor)

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_dir  <- "inputs"
plot_dir   <- "plots"
output_dir <- "outputs"

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

# ---------------------------------------------------------
# 1. Keep ONLY NNk / NN1 ratios
# ---------------------------------------------------------

ratio_NNk_NN1 <- nn_ratios_raw %>%
  mutate(
    colony = as.character(colony),
    genotype = as.character(colony_number),
    numerator = as.character(numerator),
    denominator = as.character(denominator),
    nn_rank_num = as.numeric(str_remove(numerator, "NN")),
    nn_rank = factor(numerator, levels = paste0("NN", 2:8)),
    ratio_to_NN1 = ratio_value
  ) %>%
  filter(
    denominator == "NN1",
    numerator %in% paste0("NN", 2:8),
    !is.na(ratio_to_NN1),
    is.finite(ratio_to_NN1),
    ratio_to_NN1 > 0
  )

# ---------------------------------------------------------
# 2. Raw-data summary: mean ± SE across all colonies
# ---------------------------------------------------------

ratio_summary_all <- ratio_NNk_NN1 %>%
  group_by(nn_rank, nn_rank_num) %>%
  summarise(
    n = n(),
    mean_ratio = mean(ratio_to_NN1, na.rm = TRUE),
    sd_ratio = sd(ratio_to_NN1, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(nn_rank_num)

print(ratio_summary_all)

write_csv(
  ratio_summary_all,
  file.path(output_dir, "NNk_over_NN1_all_colonies_raw_mean_SE.csv")
)

# ---------------------------------------------------------
# 3. Plot raw-data means ± SE
# ---------------------------------------------------------

p_ratio_raw_mean_se <- ggplot(
  ratio_summary_all,
  aes(x = nn_rank_num, y = mean_ratio)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(
      ymin = mean_ratio - se_ratio,
      ymax = mean_ratio + se_ratio
    ),
    width = 0.12,
    linewidth = 0.6
  ) +
  scale_x_continuous(
    breaks = 2:8,
    labels = paste0("NN", 2:8)
  ) +
  labs(
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1]),
    title = "Mean normalized nearest-neighbor spacing"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_ratio_raw_mean_se)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_all_colonies_raw_mean_SE.png"),
  plot = p_ratio_raw_mean_se,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_all_colonies_raw_mean_SE.pdf"),
  plot = p_ratio_raw_mean_se,
  width = 7,
  height = 5,
  bg = "white"
)

p_ratio_raw <- ggplot(avg_ratio_all_colonies,
                      aes(x = nn_rank_num, y = mean_NNk_over_NN1)) +
  geom_errorbar(
    aes(
      ymin = mean_NNk_over_NN1 - se,
      ymax = mean_NNk_over_NN1 + se
    ),
    width = 0.18,
    linewidth = 1.1,
    color = "black"
  ) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_point(size = 3.5, color = "black") +
  scale_x_continuous(
    breaks = 2:8,
    labels = paste0("NN", 2:8)
  ) +
  labs(
    title = "Mean normalized nearest-neighbor spacing",
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1])
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_ratio_raw)


p_ratio_raw_sd <- ggplot(avg_ratio_all_colonies,
                         aes(x = nn_rank_num, y = mean_NNk_over_NN1)) +
  geom_errorbar(
    aes(
      ymin = mean_NNk_over_NN1 - sd,
      ymax = mean_NNk_over_NN1 + sd
    ),
    width = 0.18,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_point(size = 3.5, color = "black") +
  scale_x_continuous(
    breaks = 2:8,
    labels = paste0("NN", 2:8)
  ) +
  labs(
    title = "Mean normalized nearest-neighbor spacing",
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1])
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p_ratio_raw_sd)
