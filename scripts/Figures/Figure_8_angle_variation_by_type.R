############################################################
## ANGLE VARIATION ACROSS COLONIES
## Model-based and raw-data point-range plots
##
############################################################

## =========================
## 0. PACKAGES
## =========================

packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "lme4",
  "lmerTest",
  "emmeans",
  "broom.mixed"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)

## =========================
## 1. PATHS
## =========================

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

input_file <- "triangles_full_data.xlsx"

output_dir <- "outputs_model_triangles"
plot_dir   <- "plots/model3"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 2. LOAD + CLEAN DATA
## =========================

tri <- read_excel(input_file) %>%
  clean_names() %>%
  mutate(
    genotype = factor(genotype),
    replica  = factor(replica),
    nubbin_id = interaction(genotype, replica, drop = TRUE, sep = "_"),
    day = as.numeric(day),
    
    ## Manuscript-facing label: colony
    colony = paste0("SC", as.character(genotype)),
    colony = factor(
      colony,
      levels = paste0("SC", sort(as.numeric(as.character(unique(genotype)))))
    )
  ) %>%
  filter(
    dist_to_nn1 > 0,
    dist_to_nn2 > 0,
    angle_main > 0, angle_main < 180,
    angle_nn1  > 0, angle_nn1  < 180,
    angle_nn2  > 0, angle_nn2  < 180
  ) %>%
  mutate(
    angle_sum = angle_main + angle_nn1 + angle_nn2,
    angle_sum_error = angle_sum - 180
  )

cat("\nAngle-sum sanity check:\n")
print(summary(tri$angle_sum))
print(summary(tri$angle_sum_error))

## =========================
## 3. LONG FORMAT FOR ANGLES
## =========================

angles_long <- tri %>%
  dplyr::select(
    colony,
    genotype,
    replica,
    nubbin_id,
    day,
    angle_main,
    angle_nn1,
    angle_nn2
  ) %>%
  pivot_longer(
    cols = c(angle_main, angle_nn1, angle_nn2),
    names_to = "angle_type",
    values_to = "angle"
  ) %>%
  mutate(
    angle_type = factor(
      angle_type,
      levels = c("angle_main", "angle_nn1", "angle_nn2"),
      labels = c("main", "NN1", "NN2")
    ),
    colony = factor(colony),
    day_sc = as.numeric(scale(day))
  )

## =========================
## 4. COLORS + THEME
## =========================

angle_cols <- c(
  "main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

theme_angle <- theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold", color = "black"),
    axis.text  = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

## =========================
## 5. MODEL VERSION
## =========================

options(contrasts = c("contr.sum", "contr.poly"))

m_angle3_colony <- lmer(
  angle ~ angle_type * colony + angle_type * day_sc + (1 | nubbin_id),
  data = angles_long
)

cat("\n=========================\n")
cat("MODEL SUMMARY\n")
cat("=========================\n")
print(summary(m_angle3_colony))

cat("\n=========================\n")
cat("MODEL ANOVA\n")
cat("=========================\n")
anova_m_angle3_colony <- anova(m_angle3_colony)
print(anova_m_angle3_colony)

## Model-based estimated angle means by colony
emm_angle_by_colony <- emmeans(
  m_angle3_colony,
  ~ angle_type | colony
)

pairs_angle_by_colony <- pairs(
  emm_angle_by_colony,
  adjust = "tukey"
)

cat("\n=========================\n")
cat("MODEL-BASED ANGLE MEANS BY COLONY\n")
cat("=========================\n")
print(emm_angle_by_colony)

cat("\n=========================\n")
cat("PAIRWISE ANGLE CONTRASTS WITHIN COLONY\n")
cat("=========================\n")
print(pairs_angle_by_colony)

emm_angle_by_colony_df <- as.data.frame(emm_angle_by_colony)

## Helper for CI column names, depending on emmeans output
fix_emmeans_ci <- function(df) {
  if ("asymp.LCL" %in% names(df)) {
    df <- df %>%
      rename(
        lower.CL = asymp.LCL,
        upper.CL = asymp.UCL
      )
  }
  df
}

emm_angle_by_colony_df <- emm_angle_by_colony_df %>%
  fix_emmeans_ci() %>%
  mutate(
    angle_type = factor(angle_type, levels = c("main", "NN1", "NN2")),
    colony = factor(colony, levels = levels(angles_long$colony))
  )

p_angle_colony_model <- ggplot(
  emm_angle_by_colony_df,
  aes(
    x = colony,
    y = emmean,
    color = angle_type,
    group = angle_type
  )
) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.8
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(30, 90, by = 10),
    limits = c(27, 93),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Colony",
    y = "Estimated angle (degrees)",
    color = "Angle type",
    title = "Angle variation across colonies"
  ) +
  theme_angle

print(p_angle_colony_model)

## =========================
## 6. RAW-DATA VERSION
## =========================

angle_colony_raw <- angles_long %>%
  group_by(colony, angle_type) %>%
  summarise(
    n = n(),
    mean_angle = mean(angle, na.rm = TRUE),
    sd_angle = sd(angle, na.rm = TRUE),
    se_angle = sd_angle / sqrt(n),
    lower_95 = mean_angle - 1.96 * se_angle,
    upper_95 = mean_angle + 1.96 * se_angle,
    .groups = "drop"
  ) %>%
  mutate(
    angle_type = factor(angle_type, levels = c("main", "NN1", "NN2")),
    colony = factor(colony, levels = levels(angles_long$colony))
  )

cat("\n=========================\n")
cat("RAW-DATA ANGLE SUMMARY BY COLONY\n")
cat("=========================\n")
print(angle_colony_raw)

p_angle_colony_raw <- ggplot(
  angle_colony_raw,
  aes(
    x = colony,
    y = mean_angle,
    color = angle_type,
    group = angle_type
  )
) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.8
  ) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(30, 90, by = 10),
    limits = c(27, 93),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Colony",
    y = "Observed mean angle (degrees)",
    color = "Angle type",
    title = "Angle variation across colonies"
  ) +
  theme_angle

print(p_angle_colony_raw)

## =========================
## 7. EXPORT TABLES
## =========================

write.csv(
  as.data.frame(anova_m_angle3_colony) %>%
    rownames_to_column("term"),
  file.path(output_dir, "angle_model_colony_type3_anova.csv"),
  row.names = FALSE
)

write.csv(
  emm_angle_by_colony_df,
  file.path(output_dir, "angle_model_colony_emmeans.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(pairs_angle_by_colony),
  file.path(output_dir, "angle_model_colony_pairwise_angle_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  angle_colony_raw,
  file.path(output_dir, "angle_raw_colony_summary.csv"),
  row.names = FALSE
)

## =========================
## 8. SAVE PLOTS
## =========================

ggsave(
  filename = file.path(plot_dir, "fig_8_angle_colony_model_point_range.png"),
  plot = p_angle_colony_model,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "fig_8angle_colony_model_point_range.pdf"),
  plot = p_angle_colony_model,
  width = 8,
  height = 5,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "fig_8_angle_colony_raw_point_range.png"),
  plot = p_angle_colony_raw,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "fig_8_angle_colony_raw_point_range.pdf"),
  plot = p_angle_colony_raw,
  width = 8,
  height = 5,
  bg = "white"
)

cat("\nDone. Saved model and raw-data colony angle plots.\n")
