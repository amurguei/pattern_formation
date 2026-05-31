############################################################
## Corrected normalized NNk / NN1 analyses
## Use ONLY ratios where denominator == NN1
############################################################

library(tidyverse)
library(janitor)
library(lme4)
library(lmerTest)
library(emmeans)

# ---------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_dir  <- "inputs"
output_dir <- "outputs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 1. Read ratio data
# ---------------------------------------------------------

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

# ---------------------------------------------------------
# 2. Check available ratios
# ---------------------------------------------------------

cat("\n==============================\n")
cat("RATIO TABLE CHECK\n")
cat("==============================\n")

print(
  nn_ratios_raw %>%
    count(numerator, denominator) %>%
    arrange(numerator, denominator)
)

# ---------------------------------------------------------
# 3. Prepare corrected NNk / NN1 dataset
# ---------------------------------------------------------
# IMPORTANT:
# We keep only NN2/NN1, NN3/NN1, ..., NN8/NN1.
# This excludes ratios such as NN4/NN2 or NN8/NN6.

geno_levels <- c("1", "2", "3", "4", "6", "7", "8", "9", "10")

ratio_long_corrected <- nn_ratios_raw %>%
  mutate(
    colony = as.character(colony),
    genotype = case_when(
      "colony_number" %in% names(.) ~ as.character(colony_number),
      TRUE ~ str_remove(colony, "^SC")
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

cat("\n==============================\n")
cat("CORRECTED RATIO DATASET CHECK\n")
cat("==============================\n")

cat("\nCounts by genotype and NN rank:\n")
print(table(ratio_long_corrected$genotype, ratio_long_corrected$nn_rank))

cat("\nCounts by spacing cluster and NN rank:\n")
print(table(ratio_long_corrected$spacing_cluster, ratio_long_corrected$nn_rank))

cat("\nFirst rows:\n")
print(head(ratio_long_corrected))

# ---------------------------------------------------------
# 4. Corrected normalized genotype model
# ---------------------------------------------------------

m_norm_genotype_corrected <- lmer(
  ratio_to_NN1 ~ genotype * nn_rank_num + (1 | nubbin_id),
  data = ratio_long_corrected
)

cat("\n==============================\n")
cat("CORRECTED NORMALIZED MODEL: genotype * NN rank\n")
cat("==============================\n")

print(anova(m_norm_genotype_corrected))
print(summary(m_norm_genotype_corrected))

norm_genotype_emm_corrected <- emmeans(
  m_norm_genotype_corrected,
  ~ genotype
)

norm_genotype_slopes_corrected <- emtrends(
  m_norm_genotype_corrected,
  ~ genotype,
  var = "nn_rank_num"
)

cat("\n--- Corrected normalized ratio means by genotype ---\n")
print(norm_genotype_emm_corrected)

cat("\n--- Corrected normalized ratio slopes by genotype ---\n")
print(norm_genotype_slopes_corrected)

write_csv(
  as.data.frame(anova(m_norm_genotype_corrected)) %>%
    rownames_to_column("term"),
  file.path(output_dir, "33_norm_spacing_genotype_anova_CORRECTED_NNk_over_NN1.csv")
)

write_csv(
  as.data.frame(norm_genotype_emm_corrected),
  file.path(output_dir, "33_norm_spacing_genotype_emmeans_CORRECTED_NNk_over_NN1.csv")
)

write_csv(
  as.data.frame(norm_genotype_slopes_corrected),
  file.path(output_dir, "33_norm_spacing_genotype_slopes_CORRECTED_NNk_over_NN1.csv")
)

# ---------------------------------------------------------
# 5. Corrected normalized cluster model
# ---------------------------------------------------------

m_norm_cluster_corrected <- lmer(
  ratio_to_NN1 ~ spacing_cluster * nn_rank_num + (1 | nubbin_id),
  data = ratio_long_corrected
)

cat("\n==============================\n")
cat("CORRECTED NORMALIZED MODEL: cluster * NN rank\n")
cat("==============================\n")

print(anova(m_norm_cluster_corrected))
print(summary(m_norm_cluster_corrected))

norm_cluster_emm_corrected <- emmeans(
  m_norm_cluster_corrected,
  ~ spacing_cluster
)

norm_cluster_slopes_corrected <- emtrends(
  m_norm_cluster_corrected,
  ~ spacing_cluster,
  var = "nn_rank_num"
)

cat("\n--- Corrected normalized ratio means by cluster ---\n")
print(norm_cluster_emm_corrected)
print(pairs(norm_cluster_emm_corrected))

cat("\n--- Corrected normalized ratio slopes by cluster ---\n")
print(norm_cluster_slopes_corrected)
print(pairs(norm_cluster_slopes_corrected))

write_csv(
  as.data.frame(anova(m_norm_cluster_corrected)) %>%
    rownames_to_column("term"),
  file.path(output_dir, "33_norm_spacing_cluster_anova_CORRECTED_NNk_over_NN1.csv")
)

write_csv(
  as.data.frame(norm_cluster_emm_corrected),
  file.path(output_dir, "33_norm_spacing_cluster_emmeans_CORRECTED_NNk_over_NN1.csv")
)

write_csv(
  as.data.frame(norm_cluster_slopes_corrected),
  file.path(output_dir, "33_norm_spacing_cluster_slopes_CORRECTED_NNk_over_NN1.csv")
)

# ---------------------------------------------------------
# 6. Corrected summaries for plots and David
# ---------------------------------------------------------

norm_summary_corrected <- ratio_long_corrected %>%
  group_by(spacing_cluster, nn_rank) %>%
  summarise(
    n = n(),
    mean_ratio = mean(ratio_to_NN1, na.rm = TRUE),
    sd_ratio = sd(ratio_to_NN1, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    mean_SE = sprintf("%.3f ± %.3f", mean_ratio, se_ratio),
    .groups = "drop"
  ) %>%
  mutate(
    nn_rank_num = as.numeric(str_remove(as.character(nn_rank), "NN"))
  ) %>%
  arrange(spacing_cluster, nn_rank_num)

cat("\n==============================\n")
cat("CORRECTED NORMALIZED SUMMARY BY CLUSTER AND RANK\n")
cat("==============================\n")
print(norm_summary_corrected)

write_csv(
  norm_summary_corrected,
  file.path(output_dir, "33_norm_spacing_summary_by_cluster_rank_CORRECTED_NNk_over_NN1.csv")
)

# Across all colonies
avg_ratio_all_colonies_corrected <- ratio_long_corrected %>%
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

cat("\n==============================\n")
cat("CORRECTED AVERAGE NNk / NN1 ACROSS ALL COLONIES\n")
cat("==============================\n")
print(avg_ratio_all_colonies_corrected)

write_csv(
  avg_ratio_all_colonies_corrected,
  file.path(output_dir, "average_NNk_over_NN1_all_colonies_with_SE_CORRECTED.csv")
)

# By colony
avg_ratio_by_colony_corrected <- ratio_long_corrected %>%
  group_by(genotype, colony, nn_rank, nn_rank_num) %>%
  summarise(
    n = n(),
    mean_NNk_over_NN1 = mean(ratio_to_NN1, na.rm = TRUE),
    sd = sd(ratio_to_NN1, na.rm = TRUE),
    se = sd / sqrt(n),
    mean_SE = sprintf("%.3f ± %.3f", mean_NNk_over_NN1, se),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(as.character(genotype)), nn_rank_num)

cat("\n==============================\n")
cat("CORRECTED AVERAGE NNk / NN1 BY COLONY\n")
cat("==============================\n")
print(avg_ratio_by_colony_corrected)

write_csv(
  avg_ratio_by_colony_corrected,
  file.path(output_dir, "average_NNk_over_NN1_by_colony_with_SE_CORRECTED.csv")
)

avg_ratio_by_colony_wide_corrected <- avg_ratio_by_colony_corrected %>%
  dplyr::select(colony, nn_rank, mean_SE) %>%
  pivot_wider(
    names_from = nn_rank,
    values_from = mean_SE
  ) %>%
  arrange(as.numeric(str_remove(colony, "^SC")))

cat("\n==============================\n")
cat("CORRECTED AVERAGE NNk / NN1 BY COLONY, WIDE FORMAT\n")
cat("==============================\n")
print(avg_ratio_by_colony_wide_corrected)

write_csv(
  avg_ratio_by_colony_wide_corrected,
  file.path(output_dir, "average_NNk_over_NN1_by_colony_wide_CORRECTED.csv")
)

cat("\nDone. Corrected normalized-ratio analyses saved.\n")
