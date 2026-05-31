############################################################
## Arithmetic progression of mean NNk / NN1 ratios
## Across all colonies, with SE and SD plots
############################################################

library(tidyverse)
library(janitor)
library(broom)

# ---------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_dir  <- "inputs"
output_dir <- "outputs"
plot_dir   <- "plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 1. Read data
# ---------------------------------------------------------

nn_ratios_raw <- read_csv(
  file.path(input_dir, "All_polyps_combined_NN_ratios_day165.csv"),
  show_col_types = FALSE
) %>%
  clean_names()

# ---------------------------------------------------------
# 2. Keep only NNk / NN1 ratios
# ---------------------------------------------------------
# This keeps only:
# NN2/NN1, NN3/NN1, ..., NN8/NN1

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

cat("\nCounts by NN rank:\n")
print(table(ratio_NNk_NN1$nn_rank))

# ---------------------------------------------------------
# 3. Summary across all colonies
# ---------------------------------------------------------

ratio_summary_all <- ratio_NNk_NN1 %>%
  group_by(nn_rank, nn_rank_num) %>%
  summarise(
    n = n(),
    mean_ratio = mean(ratio_to_NN1, na.rm = TRUE),
    sd_ratio = sd(ratio_to_NN1, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    mean_SE = sprintf("%.3f ± %.3f", mean_ratio, se_ratio),
    mean_SD = sprintf("%.3f ± %.3f", mean_ratio, sd_ratio),
    .groups = "drop"
  ) %>%
  arrange(nn_rank_num) %>%
  mutate(
    step_difference = mean_ratio - lag(mean_ratio)
  )

cat("\nSummary across all colonies:\n")
print(ratio_summary_all)

write_csv(
  ratio_summary_all,
  file.path(output_dir, "NNk_over_NN1_all_colonies_mean_SD_SE_steps.csv")
)

# ---------------------------------------------------------
# 4. Linear model on rank means
# ---------------------------------------------------------
# This tests whether the mean NNk/NN1 ratios follow an
# approximately linear / arithmetic progression across rank.

lm_arithmetic <- lm(
  mean_ratio ~ nn_rank_num,
  data = ratio_summary_all
)

cat("\nLinear model for arithmetic progression:\n")
print(summary(lm_arithmetic))

cat("\nModel coefficients:\n")
print(coef(lm_arithmetic))

lm_arithmetic_tidy <- tidy(lm_arithmetic)
lm_arithmetic_glance <- glance(lm_arithmetic)

write_csv(
  lm_arithmetic_tidy,
  file.path(output_dir, "NNk_over_NN1_arithmetic_linear_model_coefficients.csv")
)

write_csv(
  lm_arithmetic_glance,
  file.path(output_dir, "NNk_over_NN1_arithmetic_linear_model_fit.csv")
)

# Optional: predicted line and residuals
ratio_summary_all <- ratio_summary_all %>%
  mutate(
    predicted_mean_ratio = predict(lm_arithmetic, newdata = ratio_summary_all),
    residual = mean_ratio - predicted_mean_ratio
  )

write_csv(
  ratio_summary_all,
  file.path(output_dir, "NNk_over_NN1_all_colonies_mean_SE_SD_with_lm_predictions.csv")
)

# ---------------------------------------------------------
# 5. Shared plotting theme
# ---------------------------------------------------------

theme_spacing <- theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------
# 6. Plot with SE, visually emphasized
# ---------------------------------------------------------
# Note: SE bars are tiny because n = 3559 per rank.
# This version uses thicker bars so they are visible.

okabe_ito <- c(
  vermillion = "#D55E00",
  black = "#000000"
)

p_ratio_mean_se <- ggplot(
  ratio_summary_all,
  aes(x = nn_rank_num, y = mean_ratio)
) +
  geom_errorbar(
    aes(
      ymin = mean_ratio - se_ratio,
      ymax = mean_ratio + se_ratio
    ),
    width = 0.18,
    linewidth = 1.1,
    color = okabe_ito["black"]
  ) +
  geom_line(linewidth = 1.2, color = okabe_ito["black"]) +
  geom_point(size = 3.8, color = okabe_ito["black"]) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1,
    linetype = "dashed",
    color = okabe_ito["vermillion"]
  ) +
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

print(p_ratio_mean_se)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_arithmetic_mean_SE.png"),
  plot = p_ratio_mean_se,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_arithmetic_mean_SE.pdf"),
  plot = p_ratio_mean_se,
  width = 7,
  height = 5,
  bg = "white"
)

# ---------------------------------------------------------
# 7. Plot with SD
# ---------------------------------------------------------
# This shows spread among individual polyp-level ratios,
# not uncertainty around the mean.

p_ratio_mean_sd <- ggplot(
  ratio_summary_all,
  aes(x = nn_rank_num, y = mean_ratio)
) +
  geom_errorbar(
    aes(
      ymin = mean_ratio - sd_ratio,
      ymax = mean_ratio + sd_ratio
    ),
    width = 0.18,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_point(size = 3.8, color = "black") +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1,
    linetype = "dashed",
    color = "#D55E00"
  ) +
  scale_x_continuous(
    breaks = 2:8,
    labels = paste0("NN", 2:8)
  ) +
  labs(
    title = "Mean normalized nearest-neighbor spacing",
    x = "Nearest-neighbor rank",
    y = expression(NN[k] / NN[1])
  ) +
  theme_spacing

print(p_ratio_mean_sd)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_arithmetic_mean_SD.png"),
  plot = p_ratio_mean_sd,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_arithmetic_mean_SD.pdf"),
  plot = p_ratio_mean_sd,
  width = 7,
  height = 5,
  bg = "white"
)

# ---------------------------------------------------------
# 8. Optional: residual plot for linearity check
# ---------------------------------------------------------

p_ratio_residuals <- ggplot(
  ratio_summary_all,
  aes(x = nn_rank_num, y = residual)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_line(linewidth = 1, color = "black") +
  geom_point(size = 3.5, color = "black") +
  scale_x_continuous(
    breaks = 2:8,
    labels = paste0("NN", 2:8)
  ) +
  labs(
    title = "Residuals from linear fit",
    x = "Nearest-neighbor rank",
    y = "Residual"
  ) +
  theme_spacing

print(p_ratio_residuals)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_arithmetic_linear_fit_residuals.png"),
  plot = p_ratio_residuals,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "NNk_over_NN1_arithmetic_linear_fit_residuals.pdf"),
  plot = p_ratio_residuals,
  width = 7,
  height = 5,
  bg = "white"
)

cat("\nDone. Arithmetic progression summaries, model, and plots saved.\n")
