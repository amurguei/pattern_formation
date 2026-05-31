############################################################
## Boundary-biased positioning composite figure
## Four panels with model-based letters on one line
## Boxplots: cadetblue3 + red mean dots
## Panel B: extreme ratio jitter points hidden visually only
############################################################

# ---------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

library(tidyverse)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(patchwork)
library(grid)
library(multcomp)

input_dir  <- "inputs"
output_dir <- "outputs"
plot_dir   <- "plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 1. Read data
# ---------------------------------------------------------

boundary_df <- read_csv(
  file.path(input_dir, "PERIMETER_NUBBIN_DISTANCE_RATIO_REPORT.csv"),
  show_col_types = FALSE
)

# ---------------------------------------------------------
# 2. Clean / prepare
# ---------------------------------------------------------

geno_levels <- c("1", "2", "3", "4", "6", "7", "8", "9", "10")

boundary_df <- boundary_df %>%
  rename(
    genotype    = genotype,
    nubbin_rep  = replica,
    day         = day,
    d_perimeter = dist_perimeter,
    d_nubbin    = dist_nubbin,
    ratio       = nubbin_perimeter_ratio
  ) %>%
  mutate(
    genotype = factor(as.character(genotype), levels = geno_levels),
    nubbin_id = interaction(genotype, nubbin_rep, drop = TRUE),
    day_sc = as.numeric(scale(day)),
    closer_perimeter = ifelse(d_perimeter < d_nubbin, 1, 0),
    closer_nubbin = 1 - closer_perimeter
  ) %>%
  filter(
    !is.na(genotype),
    !is.na(day),
    !is.na(d_perimeter),
    !is.na(d_nubbin),
    !is.na(ratio)
  )

cat("\nData check:\n")
print(glimpse(boundary_df))
cat("\nN observations:", nrow(boundary_df), "\n")
cat("N nubbins:", n_distinct(boundary_df$nubbin_id), "\n")

# ---------------------------------------------------------
# 3. Models used for statistics
# ---------------------------------------------------------

m_ratio_day <- lmer(
  ratio ~ genotype + day_sc + (1 | nubbin_id),
  data = boundary_df
)

m_nubbin_day <- lmer(
  d_nubbin ~ genotype * day_sc + (1 | nubbin_id),
  data = boundary_df
)

m_perimeter_day <- lmer(
  d_perimeter ~ genotype + day_sc + (1 | nubbin_id),
  data = boundary_df
)

m_binary_day <- glmer(
  closer_perimeter ~ genotype + day_sc + (1 | nubbin_id),
  data = boundary_df,
  family = binomial
)

# ---------------------------------------------------------
# 4. Print and save model outputs
# ---------------------------------------------------------

cat("\n==============================\n")
cat("RATIO LMM: ratio ~ genotype + day_sc + (1|nubbin_id)\n")
cat("==============================\n")
print(anova(m_ratio_day))
print(summary(m_ratio_day))

cat("\n==============================\n")
cat("NUBBIN DISTANCE LMM: d_nubbin ~ genotype * day_sc + (1|nubbin_id)\n")
cat("==============================\n")
print(anova(m_nubbin_day))
print(summary(m_nubbin_day))
print(emtrends(m_nubbin_day, ~ genotype, var = "day_sc"))

cat("\n==============================\n")
cat("PERIMETER DISTANCE LMM: d_perimeter ~ genotype + day_sc + (1|nubbin_id)\n")
cat("==============================\n")
print(anova(m_perimeter_day))
print(summary(m_perimeter_day))

cat("\n==============================\n")
cat("BINARY GLMM: closer_perimeter ~ genotype + day_sc + (1|nubbin_id)\n")
cat("==============================\n")
print(summary(m_binary_day))
print(car::Anova(m_binary_day, type = 3))
print(emmeans(m_binary_day, ~ genotype, type = "response"))

boundary_lmm_results <- bind_rows(
  as.data.frame(anova(m_ratio_day)) %>%
    rownames_to_column("term") %>%
    mutate(model = "ratio_lmm"),
  as.data.frame(anova(m_nubbin_day)) %>%
    rownames_to_column("term") %>%
    mutate(model = "nubbin_distance_lmm"),
  as.data.frame(anova(m_perimeter_day)) %>%
    rownames_to_column("term") %>%
    mutate(model = "perimeter_distance_lmm")
)

write_csv(
  boundary_lmm_results,
  file.path(output_dir, "boundary_positioning_lmm_results.csv")
)

binary_glmm_results <- as.data.frame(car::Anova(m_binary_day, type = 3)) %>%
  rownames_to_column("term")

write_csv(
  binary_glmm_results,
  file.path(output_dir, "boundary_positioning_binary_glmm_results.csv")
)

# ---------------------------------------------------------
# 5. Compact letters for both boxplots
# Highest estimated mean gets "a"
# ---------------------------------------------------------

emm_ratio <- emmeans(m_ratio_day, ~ genotype)

letters_ratio <- multcomp::cld(
  emm_ratio,
  adjust = "tukey",
  Letters = letters,
  sort = TRUE,
  reversed = TRUE
) %>%
  as.data.frame() %>%
  transmute(
    genotype = factor(as.character(genotype), levels = geno_levels),
    ratio_emmean = emmean,
    ratio_letter = stringr::str_trim(.group)
  )

emm_perimeter <- emmeans(m_perimeter_day, ~ genotype)

letters_perimeter <- multcomp::cld(
  emm_perimeter,
  adjust = "tukey",
  Letters = letters,
  sort = TRUE,
  reversed = TRUE
) %>%
  as.data.frame() %>%
  transmute(
    genotype = factor(as.character(genotype), levels = geno_levels),
    perimeter_emmean = emmean,
    perimeter_letter = stringr::str_trim(.group)
  )

write_csv(
  letters_ratio,
  file.path(output_dir, "boundary_ratio_model_letters.csv")
)

write_csv(
  letters_perimeter,
  file.path(output_dir, "boundary_perimeter_distance_model_letters.csv")
)

cat("\nRatio letters:\n")
print(letters_ratio)

cat("\nPerimeter-distance letters:\n")
print(letters_perimeter)

# ---------------------------------------------------------
# 6. Single-line letter positions
# ---------------------------------------------------------

ratio_label_y_const <- 15.55

ratio_letter_pos <- letters_ratio %>%
  dplyr::select(genotype, ratio_letter) %>%
  mutate(label_y = ratio_label_y_const)

perimeter_label_y_const <- 6.6

perimeter_letter_pos <- letters_perimeter %>%
  dplyr::select(genotype, perimeter_letter) %>%
  mutate(label_y = perimeter_label_y_const)

# ---------------------------------------------------------
# 7. Panel (a): binary probability + binomial-test stars
# ---------------------------------------------------------

prob_df <- boundary_df %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    successes = sum(closer_perimeter),
    p_perimeter = mean(closer_perimeter),
    p_nubbin = mean(closer_nubbin),
    se = sd(closer_perimeter) / sqrt(n),
    p_value = binom.test(successes, n, p = 0.5)$p.value,
    signif = ifelse(p_value < 0.05, "*", ""),
    star_y = pmax(p_perimeter + se, p_nubbin + se) + 0.018,
    .groups = "drop"
  )

write_csv(
  prob_df,
  file.path(output_dir, "boundary_binary_probability_by_genotype.csv")
)

prob_long <- prob_df %>%
  dplyr::select(genotype, p_perimeter, p_nubbin, se) %>%
  pivot_longer(
    cols = c(p_perimeter, p_nubbin),
    names_to = "position",
    values_to = "p"
  ) %>%
  mutate(
    position = factor(
      position,
      levels = c("p_perimeter", "p_nubbin"),
      labels = c("Closer to perimeter", "Closer to nubbin")
    )
  )

# ---------------------------------------------------------
# 8. Time summaries for panel (c): SE and SD versions
# ---------------------------------------------------------

time_summary <- boundary_df %>%
  group_by(day) %>%
  summarise(
    n = n(),
    mean_nubbin = mean(d_nubbin, na.rm = TRUE),
    sd_nubbin = sd(d_nubbin, na.rm = TRUE),
    se_nubbin = sd_nubbin / sqrt(n),
    mean_perim = mean(d_perimeter, na.rm = TRUE),
    sd_perim = sd(d_perimeter, na.rm = TRUE),
    se_perim = sd_perim / sqrt(n),
    .groups = "drop"
  )

write_csv(
  time_summary,
  file.path(output_dir, "boundary_time_summary_SE_SD.csv")
)

# ---------------------------------------------------------
# 9. Colours
# ---------------------------------------------------------

col_perimeter_bar <- "#CC79A7"
col_nubbin_bar    <- "#56B4E9"

col_nubbin_line   <- "#56B4E9"
col_perim_line    <- "black"

col_box_cadet     <- "cadetblue3"
col_mean_red      <- "red"

# ---------------------------------------------------------
# 10. Shared theme
# ---------------------------------------------------------

theme_boundary <- theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.2, "lines"),
    legend.spacing.x = unit(0.45, "cm"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
  )

# ---------------------------------------------------------
# 11. Panel (a): probability
# ---------------------------------------------------------

pA <- ggplot(prob_long, aes(x = genotype, y = p, fill = position)) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.65,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = p - se, ymax = p + se),
    position = position_dodge(width = 0.7),
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_text(
    data = prob_df,
    aes(x = genotype, y = star_y, label = signif),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Closer to perimeter" = col_perimeter_bar,
      "Closer to nubbin"    = col_nubbin_bar
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(
    ylim = c(0, 1.06),
    clip = "on"
  ) +
  labs(
    x = "Colony",
    y = "Probability",
    fill = NULL,
    title = "(a) Probability of boundary-biased\npolyp formation"
  ) +
  theme_boundary +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# ---------------------------------------------------------
# 12. Panel (b): ratio boxplot with single-line letters
# NOTE:
# - Boxplot/stat_summary/letters use full boundary_df.
# - Jitter points are visually filtered at ratio <= 16 to avoid
#   points floating outside the panel.
# ---------------------------------------------------------

pB <- ggplot(boundary_df, aes(x = genotype, y = ratio)) +
  geom_boxplot(
    fill = col_box_cadet,
    color = "black",
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    data = boundary_df %>% filter(ratio <= 16),
    aes(x = genotype, y = ratio),
    width = 0.12,
    alpha = 0.18,
    size = 1.2,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = col_mean_red,
    color = "black",
    stroke = 0.5
  ) +
  geom_text(
    data = ratio_letter_pos,
    aes(x = genotype, y = label_y, label = ratio_letter),
    inherit.aes = FALSE,
    size = 5.5,
    vjust = 0
  ) +
  scale_y_continuous(
    breaks = seq(0, 16, by = 4),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(
    ylim = c(0, 16),
    clip = "on"
  ) +
  labs(
    x = "Colony",
    y = expression(d[nubbin] / d[perimeter]),
    title = "(b) Relative positioning\nof new polyps"
  ) +
  theme_boundary +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# ---------------------------------------------------------
# 13. Panel (c): time dynamics, SE or SD
# ---------------------------------------------------------

make_time_panel <- function(error_type = c("SE", "SD")) {
  
  error_type <- match.arg(error_type)
  
  if (error_type == "SE") {
    time_df <- time_summary %>%
      mutate(
        nubbin_lower = mean_nubbin - se_nubbin,
        nubbin_upper = mean_nubbin + se_nubbin,
        perim_lower = mean_perim - se_perim,
        perim_upper = mean_perim + se_perim
      )
  } else {
    time_df <- time_summary %>%
      mutate(
        nubbin_lower = mean_nubbin - sd_nubbin,
        nubbin_upper = mean_nubbin + sd_nubbin,
        perim_lower = mean_perim - sd_perim,
        perim_upper = mean_perim + sd_perim
      )
  }
  
  ggplot(time_df, aes(x = day)) +
    geom_line(aes(y = mean_nubbin, color = "Nubbin"), linewidth = 1) +
    geom_point(aes(y = mean_nubbin, color = "Nubbin"), size = 2.8) +
    geom_errorbar(
      aes(ymin = nubbin_lower, ymax = nubbin_upper, color = "Nubbin"),
      width = 1.2,
      linewidth = 0.6,
      alpha = 0.8
    ) +
    geom_line(aes(y = mean_perim, color = "Perimeter"), linewidth = 1) +
    geom_point(aes(y = mean_perim, color = "Perimeter"), size = 2.8) +
    geom_errorbar(
      aes(ymin = perim_lower, ymax = perim_upper, color = "Perimeter"),
      width = 1.2,
      linewidth = 0.6,
      alpha = 0.8
    ) +
    scale_color_manual(
      values = c(
        "Nubbin" = col_nubbin_line,
        "Perimeter" = col_perim_line
      )
    ) +
    labs(
      x = "Day",
      y = "Mean distance (mm)",
      color = NULL,
      title = "(c) Temporal dynamics of\nboundary-biased positioning"
    ) +
    theme_boundary +
    theme(
      legend.position = "top"
    )
}

pC_se <- make_time_panel("SE")
pC_sd <- make_time_panel("SD")

# ---------------------------------------------------------
# 14. Panel (d): perimeter-distance boxplot with single-line letters
# ---------------------------------------------------------

pD <- ggplot(boundary_df, aes(x = genotype, y = d_perimeter)) +
  geom_boxplot(
    fill = col_box_cadet,
    color = "black",
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.18,
    size = 1.2,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = col_mean_red,
    color = "black",
    stroke = 0.5
  ) +
  geom_text(
    data = perimeter_letter_pos,
    aes(x = genotype, y = label_y, label = perimeter_letter),
    inherit.aes = FALSE,
    size = 5.5,
    vjust = 0
  ) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 2),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(
    ylim = c(0, 7),
    clip = "on"
  ) +
  labs(
    x = "Colony",
    y = expression(d[perimeter]~"(mm)"),
    title = "(d) Distance to colony perimeter"
  ) +
  theme_boundary +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# ---------------------------------------------------------
# 15. Composite figures: SE and SD versions
# ---------------------------------------------------------

figure_boundary_se <- (pA | pB) / (pC_se | pD) +
  plot_layout(
    heights = c(1, 1.05),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(10, 16, 10, 16)
  )

figure_boundary_sd <- (pA | pB) / (pC_sd | pD) +
  plot_layout(
    heights = c(1, 1.05),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(10, 16, 10, 16)
  )

print(figure_boundary_se)
print(figure_boundary_sd)

# ---------------------------------------------------------
# 16. Save composite figures
# ---------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SE_FINAL.png"),
  plot = figure_boundary_se,
  width = 15,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SE_FINAL.pdf"),
  plot = figure_boundary_se,
  width = 15,
  height = 10,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SD_FINAL.png"),
  plot = figure_boundary_sd,
  width = 15,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SD_FINAL.pdf"),
  plot = figure_boundary_sd,
  width = 15,
  height = 10,
  bg = "white"
)

# ---------------------------------------------------------
# 17. Save individual panels
# ---------------------------------------------------------

ggsave(
  file.path(plot_dir, "boundary_panel_a_binary_probability_FINAL.png"),
  pA,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_a_binary_probability_FINAL.pdf"),
  pA,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_b_ratio_letters_FINAL.png"),
  pB,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_b_ratio_letters_FINAL.pdf"),
  pB,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_c_time_SE_FINAL.png"),
  pC_se,
  width = 8,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_c_time_SE_FINAL.pdf"),
  pC_se,
  width = 8,
  height = 5.5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_c_time_SD_FINAL.png"),
  pC_sd,
  width = 8,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_c_time_SD_FINAL.pdf"),
  pC_sd,
  width = 8,
  height = 5.5,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_d_perimeter_letters_FINAL.png"),
  pD,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_d_perimeter_letters_FINAL.pdf"),
  pD,
  width = 7,
  height = 5,
  bg = "white"
)

cat("\nDone. Final four-panel boundary-positioning figures saved.\n")