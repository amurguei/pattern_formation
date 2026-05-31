############################################################
## Figure S18. Absolute and proportional scaling of polyp spacing
## Corrected version: normalized ratios include ONLY NNk / NN1
############################################################

library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(grid)

# ---------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_dir <- "inputs"
plot_dir  <- "plots"

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 1. Read data
# ---------------------------------------------------------

nn_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN1to8_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

geno_levels <- c("1", "2", "3", "4", "6", "7", "8", "9", "10")

group_cols <- c(
  "Low-spacing"  = "#56B4E9",
  "High-spacing" = "#CC79A7"
)

# ---------------------------------------------------------
# 2. Prepare absolute NN1-NN8 data
# ---------------------------------------------------------

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
    nn_rank = factor(paste0("NN", nn_rank_num), levels = paste0("NN", 1:8))
  ) %>%
  filter(
    !is.na(genotype),
    !is.na(distance),
    is.finite(distance),
    distance > 0
  )

# ---------------------------------------------------------
# 3. Cluster colonies by mean absolute NN1-NN8 profile
# ---------------------------------------------------------

profile_mat <- nn_long %>%
  group_by(genotype, nn_rank_num) %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = nn_rank_num,
    values_from = mean_distance,
    names_prefix = "NN"
  ) %>%
  column_to_rownames("genotype") %>%
  as.matrix()

hc_abs <- hclust(dist(profile_mat), method = "ward.D2")
cluster_abs <- cutree(hc_abs, k = 2)

cluster_df <- tibble(
  genotype = names(cluster_abs),
  spacing_cluster_raw = paste0("Cluster_", cluster_abs)
)

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
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    )
  )

cat("\nCluster assignments:\n")
print(cluster_df)

# ---------------------------------------------------------
# 4. Absolute model estimates
# ---------------------------------------------------------

nn_long_cluster <- nn_long %>%
  left_join(
    cluster_df %>% dplyr::select(genotype, spacing_cluster),
    by = "genotype"
  ) %>%
  mutate(
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    )
  )

m_abs_cluster <- lmer(
  distance ~ spacing_cluster * nn_rank_num + (1 | nubbin_id),
  data = nn_long_cluster
)

cat("\nAbsolute model ANOVA:\n")
print(anova(m_abs_cluster))

abs_cluster_emm <- as.data.frame(
  emmeans(m_abs_cluster, ~ spacing_cluster)
)

abs_cluster_slopes <- as.data.frame(
  emtrends(m_abs_cluster, ~ spacing_cluster, var = "nn_rank_num")
)

cat("\nAbsolute model estimated means:\n")
print(abs_cluster_emm)

cat("\nAbsolute model estimated slopes:\n")
print(abs_cluster_slopes)

# ---------------------------------------------------------
# 5. Prepare corrected normalized NNk / NN1 data
# ---------------------------------------------------------
# IMPORTANT:
# This keeps only NN2/NN1, NN3/NN1, ..., NN8/NN1.
# It excludes all other pairwise ratios such as NN4/NN2.

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
    nubbin_id = interaction(genotype, slide, drop = TRUE),
    numerator = as.character(numerator),
    denominator = as.character(denominator),
    nn_rank_num = as.numeric(str_remove(numerator, "NN")),
    nn_rank = factor(numerator, levels = paste0("NN", 2:8)),
    ratio_to_NN1 = ratio_value
  ) %>%
  filter(
    denominator == "NN1",
    numerator %in% paste0("NN", 2:8),
    !is.na(genotype),
    !is.na(nn_rank_num),
    !is.na(ratio_to_NN1),
    is.finite(ratio_to_NN1),
    ratio_to_NN1 > 0
  ) %>%
  left_join(
    cluster_df %>% dplyr::select(genotype, spacing_cluster),
    by = "genotype"
  ) %>%
  mutate(
    spacing_cluster = factor(
      spacing_cluster,
      levels = c("Low-spacing", "High-spacing")
    )
  )

cat("\nCorrected ratio counts by spacing cluster and rank:\n")
print(table(ratio_long$spacing_cluster, ratio_long$nn_rank))

# ---------------------------------------------------------
# 6. Corrected normalized model estimates
# ---------------------------------------------------------

m_norm_cluster <- lmer(
  ratio_to_NN1 ~ spacing_cluster * nn_rank_num + (1 | nubbin_id),
  data = ratio_long
)

cat("\nCorrected normalized model ANOVA:\n")
print(anova(m_norm_cluster))

norm_cluster_emm <- as.data.frame(
  emmeans(m_norm_cluster, ~ spacing_cluster)
)

norm_cluster_slopes <- as.data.frame(
  emtrends(m_norm_cluster, ~ spacing_cluster, var = "nn_rank_num")
)

cat("\nCorrected normalized model estimated means:\n")
print(norm_cluster_emm)

cat("\nCorrected normalized model estimated slopes:\n")
print(norm_cluster_slopes)

# ---------------------------------------------------------
# 7. Descriptive summaries for line panels
# ---------------------------------------------------------

abs_summary <- nn_long_cluster %>%
  group_by(spacing_cluster, nn_rank) %>%
  summarise(
    n = n(),
    mean_distance = mean(distance, na.rm = TRUE),
    sd_distance = sd(distance, na.rm = TRUE),
    se_distance = sd_distance / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    nn_rank_num = as.numeric(str_remove(as.character(nn_rank), "NN"))
  )

norm_summary <- ratio_long %>%
  group_by(spacing_cluster, nn_rank) %>%
  summarise(
    n = n(),
    mean_ratio = mean(ratio_to_NN1, na.rm = TRUE),
    sd_ratio = sd(ratio_to_NN1, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    nn_rank_num = as.numeric(str_remove(as.character(nn_rank), "NN"))
  )

# ---------------------------------------------------------
# 8. Plot theme
# ---------------------------------------------------------

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
# 9. Panel (a): absolute spacing lines
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
# 10. Panel (b): absolute model estimates
# ---------------------------------------------------------

abs_estimates <- bind_rows(
  abs_cluster_emm %>%
    transmute(
      spacing_cluster,
      metric = "Mean distance",
      estimate = emmean,
      lower = asymp.LCL,
      upper = asymp.UCL
    ),
  abs_cluster_slopes %>%
    transmute(
      spacing_cluster,
      metric = "Scaling slope",
      estimate = nn_rank_num.trend,
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
# 11. Panel (c): corrected proportional spacing lines
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
# 12. Panel (d): corrected normalized model estimates
# ---------------------------------------------------------

norm_estimates <- bind_rows(
  norm_cluster_emm %>%
    transmute(
      spacing_cluster,
      metric = "mean_ratio",
      estimate = emmean,
      lower = asymp.LCL,
      upper = asymp.UCL
    ),
  norm_cluster_slopes %>%
    transmute(
      spacing_cluster,
      metric = "slope",
      estimate = nn_rank_num.trend,
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
# 13. Composite figure
# ---------------------------------------------------------

figure_spacing_scaling <- (pA | pB) / (pC | pD) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top"
  )

print(figure_spacing_scaling)

# ---------------------------------------------------------
# 14. Save composite
# ---------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "Fig_S18_absolute_proportional_scaling_CORRECTED.png"),
  plot = figure_spacing_scaling,
  width = 12,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "Fig_S18_absolute_proportional_scaling_CORRECTED.pdf"),
  plot = figure_spacing_scaling,
  width = 12,
  height = 9,
  bg = "white"
)

# ---------------------------------------------------------
# 15. Save individual panels
# ---------------------------------------------------------

ggsave(
  file.path(plot_dir, "Fig_S18a_absolute_spacing_CORRECTED.png"),
  pA,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18a_absolute_spacing_CORRECTED.pdf"),
  pA,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18b_absolute_model_estimates_CORRECTED.png"),
  pB,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18b_absolute_model_estimates_CORRECTED.pdf"),
  pB,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18c_proportional_spacing_CORRECTED.png"),
  pC,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18c_proportional_spacing_CORRECTED.pdf"),
  pC,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18d_normalized_model_estimates_CORRECTED.png"),
  pD,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "Fig_S18d_normalized_model_estimates_CORRECTED.pdf"),
  pD,
  width = 7,
  height = 5,
  bg = "white"
)

cat("\nDone. Corrected Figure S18 saved to plots/.\n")
