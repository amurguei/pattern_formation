############################################################
## Boundary-biased positioning composite figure
## Cleaner titles + legends moved upward
## No letters in boxplots
############################################################

# ---------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

library(tidyverse)
library(readr)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(patchwork)
library(grid)

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

cat("\nN observations:", nrow(boundary_df), "\n")
cat("N nubbins:", n_distinct(boundary_df$nubbin_id), "\n")

# ---------------------------------------------------------
# 3. Models kept for reporting
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

ratio_anova <- as.data.frame(anova(m_ratio_day)) %>%
  rownames_to_column("term")

nubbin_anova <- as.data.frame(anova(m_nubbin_day)) %>%
  rownames_to_column("term")

perimeter_anova <- as.data.frame(anova(m_perimeter_day)) %>%
  rownames_to_column("term")

binary_anova <- as.data.frame(car::Anova(m_binary_day, type = 3)) %>%
  rownames_to_column("term")

write_csv(ratio_anova, file.path(output_dir, "boundary_ratio_lmm_anova.csv"))
write_csv(nubbin_anova, file.path(output_dir, "boundary_nubbin_lmm_anova.csv"))
write_csv(perimeter_anova, file.path(output_dir, "boundary_perimeter_lmm_anova.csv"))
write_csv(binary_anova, file.path(output_dir, "boundary_binary_glmm_anova.csv"))

# ---------------------------------------------------------
# 4. Panel A data
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
  select(genotype, p_perimeter, p_nubbin, se) %>%
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
# 5. Panel C data
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
# 6. Colors
# ---------------------------------------------------------

col_perimeter_bar <- "#CC79A7"
col_nubbin_bar    <- "#56B4E9"

col_nubbin_line   <- "#56B4E9"
col_perim_line    <- "black"

col_box_fill      <- "cadetblue3"
col_mean_red      <- "red"

# ---------------------------------------------------------
# 7. Compact shared theme
# ---------------------------------------------------------

theme_boundary <- theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.85, "lines"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-8, 0, -8, 0),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 6, r = 12, b = 6, l = 8)
  )

# ---------------------------------------------------------
# 8. Panel A: location of new polyps
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
      "Closer to nubbin" = col_nubbin_bar
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.025))
  ) +
  coord_cartesian(
    ylim = c(0, 1.05),
    clip = "on"
  ) +
  labs(
    x = "Colony",
    y = "Probability",
    fill = NULL,
    title = "(a) Location of new polyps"
  ) +
  theme_boundary +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = c(0.62, 0.995),
    legend.justification = c(0.5, 1),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.90), color = NA)
  )

# ---------------------------------------------------------
# 9. Panel B: relative positioning
# ---------------------------------------------------------

pB <- ggplot(boundary_df, aes(x = genotype, y = ratio)) +
  geom_boxplot(
    fill = col_box_fill,
    color = "black",
    width = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    data = boundary_df %>% filter(ratio <= 16),
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
  scale_y_continuous(
    breaks = seq(0, 16, by = 4),
    expand = expansion(mult = c(0.01, 0.025))
  ) +
  coord_cartesian(
    ylim = c(0, 16),
    clip = "on"
  ) +
  labs(
    x = "Colony",
    y = expression(d[nubbin] / d[perimeter]),
    title = "(b) Relative positioning of new polyps"
  ) +
  theme_boundary +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# ---------------------------------------------------------
# 10. Panel C: temporal dynamics
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
    scale_x_continuous(
      breaks = c(60, 90, 120, 150)
    ) +
    labs(
      x = "Day",
      y = "Mean distance (mm)",
      color = NULL,
      title = "(c) Temporal dynamics of new polyp formation"
    ) +
    theme_boundary +
    theme(
      legend.position = c(0.55, 0.995),
      legend.justification = c(0.5, 1),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = alpha("white", 0.90), color = NA)
    )
}

pC_se <- make_time_panel("SE")
pC_sd <- make_time_panel("SD")

# ---------------------------------------------------------
# 11. Panel D: distance to perimeter
# ---------------------------------------------------------

pD <- ggplot(boundary_df, aes(x = genotype, y = d_perimeter)) +
  geom_boxplot(
    fill = col_box_fill,
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
  scale_y_continuous(
    breaks = seq(0, 6, by = 2),
    expand = expansion(mult = c(0.01, 0.025))
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
# 12. Composite figures
# ---------------------------------------------------------

figure_boundary_se <- (pA | pB) / (pC_se | pD) +
  plot_layout(
    heights = c(1, 1),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(4, 6, 4, 6)
  )

figure_boundary_sd <- (pA | pB) / (pC_sd | pD) +
  plot_layout(
    heights = c(1, 1),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(4, 6, 4, 6)
  )

print(figure_boundary_se)
print(figure_boundary_sd)

# ---------------------------------------------------------
# 13. Save composite figures
# ---------------------------------------------------------

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SE_clean_titles.png"),
  plot = figure_boundary_se,
  width = 14.5,
  height = 8.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SE_clean_titles.pdf"),
  plot = figure_boundary_se,
  width = 14.5,
  height = 8.2,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SD_clean_titles.png"),
  plot = figure_boundary_sd,
  width = 14.5,
  height = 8.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "boundary_positioning_composite_4panel_SD_clean_titles.pdf"),
  plot = figure_boundary_sd,
  width = 14.5,
  height = 8.2,
  bg = "white"
)

# ---------------------------------------------------------
# 14. Save individual panels
# ---------------------------------------------------------

ggsave(
  file.path(plot_dir, "boundary_panel_a_location_new_polyps.png"),
  pA,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_b_relative_positioning.png"),
  pB,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_c_temporal_SE.png"),
  pC_se,
  width = 8,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_c_temporal_SD.png"),
  pC_sd,
  width = 8,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(plot_dir, "boundary_panel_d_distance_perimeter.png"),
  pD,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

cat("\nDone: clean-title compact boundary-positioning figures saved.\n")

