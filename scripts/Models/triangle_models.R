############################################################
## MODEL 3 ONLY: FINAL ANGLE ANALYSIS
## Keeps one model and extracts everything needed
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
  "patchwork",
  "car",
  "broom.mixed",
  "scales"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)
library(car)
library(broom.mixed)
library(scales)

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

## =========================
## 1. LOAD + CLEAN
## =========================
file_path <- "triangles_full_data.xlsx"

tri <- read_excel(file_path) %>%
  clean_names() %>%
  mutate(
    genotype = factor(genotype),
    replica  = factor(replica),
    nubbin_id = interaction(genotype, replica, drop = TRUE, sep = "_"),
    day = as.numeric(day)
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

## Optional sanity check
summary(tri$angle_sum)
summary(tri$angle_sum_error)

## =========================
## 2. LONG FORMAT
## =========================
angles_long <- tri %>%
  select(genotype, replica, nubbin_id, day,
         angle_main, angle_nn1, angle_nn2) %>%
  pivot_longer(
    cols = c(angle_main, angle_nn1, angle_nn2),
    names_to = "angle_type",
    values_to = "angle"
  ) %>%
  mutate(
    angle_type = factor(
      angle_type,
      levels = c("angle_main", "angle_nn1", "angle_nn2"),
      labels = c("Main", "NN1", "NN2")
    ),
    genotype = factor(genotype)
  )

## Set contrasts for Type III tests
options(contrasts = c("contr.sum", "contr.poly"))

## =========================
## 3. FINAL MODEL (MODEL 3)
## =========================
m_angle3 <- lmer(
  angle ~ angle_type * genotype + angle_type * scale(day) + (1 | nubbin_id),
  data = angles_long
)

cat("\n=========================\n")
cat("MODEL 3 SUMMARY\n")
cat("=========================\n")
print(summary(m_angle3))

cat("\n=========================\n")
cat("MODEL 3 TYPE III ANOVA\n")
cat("=========================\n")
anova_m3 <- anova(m_angle3)
print(anova_m3)

## =========================
## 4. OVERALL ANGLE DIFFERENCES
## (Model 1 equivalent, but from Model 3)
## =========================
emm_angle_overall <- emmeans(m_angle3, ~ angle_type)
pairs_angle_overall <- pairs(emm_angle_overall, adjust = "tukey")

cat("\n=========================\n")
cat("OVERALL ANGLE-TYPE MEANS\n")
cat("=========================\n")
print(emm_angle_overall)

cat("\n=========================\n")
cat("OVERALL ANGLE-TYPE PAIRWISE CONTRASTS\n")
cat("=========================\n")
print(pairs_angle_overall)

## =========================
## 5. ANGLE DIFFERENCES WITHIN GENOTYPE
## (Model 2 equivalent)
## =========================
emm_angle_by_genotype <- emmeans(m_angle3, ~ angle_type | genotype)
pairs_angle_by_genotype <- pairs(emm_angle_by_genotype, adjust = "tukey")

cat("\n=========================\n")
cat("ANGLE TYPE WITHIN GENOTYPE\n")
cat("=========================\n")
print(emm_angle_by_genotype)

cat("\n=========================\n")
cat("PAIRWISE ANGLE-TYPE CONTRASTS WITHIN GENOTYPE\n")
cat("=========================\n")
print(pairs_angle_by_genotype)

## =========================
## 6. GENOTYPE DIFFERENCES WITHIN ANGLE TYPE
## =========================
emm_genotype_by_angle <- emmeans(m_angle3, ~ genotype | angle_type)
pairs_genotype_by_angle <- pairs(emm_genotype_by_angle, adjust = "tukey")

cat("\n=========================\n")
cat("GENOTYPE WITHIN ANGLE TYPE\n")
cat("=========================\n")
print(emm_genotype_by_angle)

cat("\n=========================\n")
cat("PAIRWISE GENOTYPE CONTRASTS WITHIN ANGLE TYPE\n")
cat("=========================\n")
print(pairs_genotype_by_angle)

## =========================
## 7. TEMPORAL TRENDS BY ANGLE TYPE
## =========================
angle_slopes <- emtrends(m_angle3, ~ angle_type, var = "day")
pairs_angle_slopes <- pairs(angle_slopes, adjust = "tukey")

cat("\n=========================\n")
cat("TEMPORAL SLOPES BY ANGLE TYPE\n")
cat("=========================\n")
print(angle_slopes)

cat("\n=========================\n")
cat("PAIRWISE COMPARISONS OF SLOPES\n")
cat("=========================\n")
print(pairs_angle_slopes)

## =========================
## 8. ANGLE DIFFERENCES AT REPRESENTATIVE TIMES
## early / mid / late from observed day distribution
## =========================
angles_long <- angles_long %>%
  mutate(
    day_sc = as.numeric(scale(day))
  )

m_angle3 <- lmer(
  angle ~ angle_type * genotype + angle_type * day_sc + (1 | nubbin_id),
  data = angles_long
)


## =========================
## 8. ANGLE DIFFERENCES AT REPRESENTATIVE TIMES
## early / mid / late from observed day distribution
## =========================
day_mean <- mean(angles_long$day, na.rm = TRUE)
day_sd   <- sd(angles_long$day, na.rm = TRUE)

time_points_raw <- c(
  Early = day_mean - day_sd,
  Mid   = day_mean,
  Late  = day_mean + day_sd
)

time_points_scaled <- c(-1, 0, 1)

emm_angle_by_time <- emmeans(
  m_angle3,
  ~ angle_type | day_sc,
  at = list(day_sc = time_points_scaled)
)

emm_angle_by_time_df <- as.data.frame(emm_angle_by_time) %>%
  mutate(
    time_label = rep(names(time_points_raw), each = 3),
    day_value  = rep(unname(round(time_points_raw, 1)), each = 3)
  )

pairs_angle_by_time <- pairs(emm_angle_by_time, adjust = "tukey")

cat("\n=========================\n")
cat("ANGLE TYPE AT EARLY / MID / LATE TIMES\n")
cat("=========================\n")
print(emm_angle_by_time_df)

cat("\n=========================\n")
cat("PAIRWISE ANGLE CONTRASTS AT EARLY / MID / LATE TIMES\n")
cat("=========================\n")
print(pairs_angle_by_time)


## =========================
## 9. EARLY / MID / LATE WITHIN EACH ANGLE TYPE
## =========================
emm_period_by_angle <- emmeans(
  m_angle3,
  ~ day_sc | angle_type,
  at = list(day_sc = time_points_scaled)
)

emm_period_by_angle_df <- as.data.frame(emm_period_by_angle) %>%
  mutate(
    time_label = rep(names(time_points_raw), times = 3),
    day_value  = rep(unname(round(time_points_raw, 1)), times = 3)
  )

pairs_period_by_angle <- pairs(emm_period_by_angle, adjust = "tukey")

cat("\n=========================\n")
cat("EARLY / MID / LATE WITHIN EACH ANGLE TYPE\n")
cat("=========================\n")
print(emm_period_by_angle_df)

cat("\n=========================\n")
cat("PAIRWISE EARLY / MID / LATE CONTRASTS WITHIN EACH ANGLE TYPE\n")
cat("=========================\n")
print(pairs_period_by_angle)
## =========================
## 10. EXPORT TABLES
## =========================
dir.create("outputs_model_triangles", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

write.csv(as.data.frame(anova_m3),
          "outputs_model_triangles/model3_type3_anova.csv")

write.csv(as.data.frame(emm_angle_overall),
          "outputs_model_triangles/model3_emmeans_angle_overall.csv",
          row.names = FALSE)

write.csv(as.data.frame(pairs_angle_overall),
          "outputs_model_triangles/model3_pairs_angle_overall.csv",
          row.names = FALSE)

write.csv(as.data.frame(emm_angle_by_genotype),
          "outputs_model_triangles/model3_emmeans_angle_by_genotype.csv",
          row.names = FALSE)

write.csv(as.data.frame(pairs_angle_by_genotype),
          "outputs_model_triangles/model3_pairs_angle_by_genotype.csv",
          row.names = FALSE)

write.csv(as.data.frame(emm_genotype_by_angle),
          "outputs_model_triangles/model3_emmeans_genotype_by_angle.csv",
          row.names = FALSE)

write.csv(as.data.frame(pairs_genotype_by_angle),
          "outputs_model_triangles/model3_pairs_genotype_by_angle.csv",
          row.names = FALSE)

write.csv(as.data.frame(angle_slopes),
          "outputs_model_triangles/model3_angle_slopes.csv",
          row.names = FALSE)

write.csv(as.data.frame(pairs_angle_slopes),
          "outputs_model_triangles/model3_pairs_angle_slopes.csv",
          row.names = FALSE)

write.csv(emm_angle_by_time_df,
          "outputs_model_triangles/model3_emmeans_angle_by_time.csv",
          row.names = FALSE)

write.csv(as.data.frame(pairs_angle_by_time),
          "outputs_model_triangles/model3_pairs_angle_by_time.csv",
          row.names = FALSE)

write.csv(emm_period_by_angle_df,
          "outputs_model_triangles/model3_emmeans_period_by_angle.csv",
          row.names = FALSE)

write.csv(as.data.frame(pairs_period_by_angle),
          "outputs_model_triangles/model3_pairs_period_by_angle.csv",
          row.names = FALSE)

## =========================
## 11. PLOTS
## =========================

angle_cols <- c(
  "Main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

## A. overall violin
p_angle_violin <- ggplot(
  angles_long,
  aes(x = angle_type, y = angle, fill = angle_type)
) +
  geom_violin(alpha = 0.5, trim = FALSE, color = NA) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.alpha = 0.2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.8,
    color = "black"
  ) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "Angle (degrees)",
    title = "Angles differ by triangle vertex"
  )

## B. genotype-specific estimated means
emm_gen_plot <- as.data.frame(emm_angle_by_genotype)

p_angle_genotype <- ggplot(
  emm_gen_plot,
  aes(x = genotype, y = emmean, color = angle_type, group = angle_type)
) +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(width = 0.4),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Genotype",
    y = "Estimated marginal mean angle (degrees)",
    color = "Angle type",
    title = "Genotype-specific angle means from Model 3"
  )

## C. model-based time trajectories
newdat <- expand.grid(
  angle_type = factor(c("Main", "NN1", "NN2"), levels = c("Main", "NN1", "NN2")),
  genotype = levels(angles_long$genotype)[1],   # genotype main effect is negligible
  day = seq(min(angles_long$day), max(angles_long$day), length.out = 200),
  nubbin_id = angles_long$nubbin_id[1]
)

## C. model-based time trajectories
m_angle3_clean <- lmer(
  angle ~ angle_type * genotype + angle_type * day_sc + (1 | nubbin_id),
  data = angles_long
)

day_seq <- seq(min(angles_long$day), max(angles_long$day), length.out = 200)

newdat <- expand.grid(
  angle_type = factor(c("Main", "NN1", "NN2"), levels = c("Main", "NN1", "NN2")),
  genotype = levels(angles_long$genotype)[1],   # genotype main effect is negligible
  day = day_seq,
  nubbin_id = angles_long$nubbin_id[1]
) %>%
  mutate(
    day_sc = (day - mean(angles_long$day, na.rm = TRUE)) / sd(angles_long$day, na.rm = TRUE)
  )

newdat$pred <- predict(m_angle3_clean, newdata = newdat, re.form = NA)
## D. early / mid / late summary plot
p_time_slices <- ggplot(
  emm_angle_by_time_df,
  aes(x = time_label, y = emmean, color = angle_type, group = angle_type)
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Estimated angle (degrees)",
    color = "Angle type",
    title = "Model 3 angle estimates at representative times"
  )
names(emm_gen_plot)

print(p_angle_violin)
print(p_angle_genotype)
print(p_time_model)
print(p_time_slices)

ggsave("plots/model3_angle_violin.png", p_angle_violin, width = 6, height = 5, dpi = 600, bg = "white")
ggsave("plots/model3_angle_genotype.png", p_angle_genotype, width = 8, height = 5, dpi = 600, bg = "white")
ggsave("plots/model3_time_model.png", p_time_model, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("plots/model3_time_slices.png", p_time_slices, width = 7, height = 5, dpi = 600, bg = "white")

ggsave("plots/model3_angle_violin.pdf", p_angle_violin, width = 6, height = 5, bg = "white")
ggsave("plots/model3_angle_genotype.pdf", p_angle_genotype, width = 8, height = 5, bg = "white")
ggsave("plots/model3_time_model.pdf", p_time_model, width = 7, height = 5, bg = "white")
ggsave("plots/model3_time_slices.pdf", p_time_slices, width = 7, height = 5, bg = "white")



## =========================
## 0. PREP: clean time variable
## =========================
angles_long <- angles_long %>%
  mutate(day_sc = as.numeric(scale(day)))

## =========================
## 1. REFIT MODEL 3 (clean version)
## =========================
m_angle3_clean <- lmer(
  angle ~ angle_type * genotype + angle_type * day_sc + (1 | nubbin_id),
  data = angles_long
)

## =========================
## 2. Helper: fix emmeans CI names
## =========================
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

## =========================
## 3. GENOTYPE EFFECT PLOT
## =========================
emm_angle_by_genotype <- emmeans(
  m_angle3_clean,
  ~ angle_type | genotype
)

emm_gen_plot <- as.data.frame(emm_angle_by_genotype) %>%
  fix_emmeans_ci()

p_angle_genotype <- ggplot(
  emm_gen_plot,
  aes(x = genotype, y = emmean, color = angle_type, group = angle_type)
) +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(width = 0.4),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 10),
    limits = c(30, 90)
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Genotype",
    y = "Estimated angle (degrees)",
    color = "Angle type",
    title = "Angle variation across genotypes"
  )

print(p_angle_genotype)

## =========================
## 4. REPRESENTATIVE TIME POINTS (early / mid / late)
## =========================
day_mean <- mean(angles_long$day, na.rm = TRUE)
day_sd   <- sd(angles_long$day, na.rm = TRUE)

time_points_raw <- c(
  Early = day_mean - day_sd,
  Mid   = day_mean,
  Late  = day_mean + day_sd
)

time_points_scaled <- c(-1, 0, 1)

emm_angle_by_time <- emmeans(
  m_angle3_clean,
  ~ angle_type | day_sc,
  at = list(day_sc = time_points_scaled)
)

emm_time_plot <- as.data.frame(emm_angle_by_time) %>%
  fix_emmeans_ci() %>%
  mutate(
    time_label = rep(names(time_points_raw), each = 3),
    day_value  = rep(round(time_points_raw, 1), each = 3)
  )

## =========================
## 5. TIME-SLICE PLOT
## =========================
p_time_slices <- ggplot(
  emm_time_plot,
  aes(x = time_label, y = emmean, color = angle_type, group = angle_type)
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.15,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180)
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Estimated angle (degrees)",
    color = "Angle type",
    title = "Angle changes through time (Model 3)"
  )

print(p_time_slices)

## =========================
## 6. FULL TRAJECTORY (smooth curves)
## =========================
day_seq <- seq(min(angles_long$day), max(angles_long$day), length.out = 200)

newdat <- expand.grid(
  angle_type = factor(c("Main", "NN1", "NN2"), levels = c("Main", "NN1", "NN2")),
  genotype = levels(angles_long$genotype)[1],
  day = day_seq,
  nubbin_id = angles_long$nubbin_id[1]
) %>%
  mutate(
    day_sc = (day - day_mean) / day_sd
  )

newdat$pred <- predict(m_angle3_clean, newdata = newdat, re.form = NA)

p_time_model <- ggplot(
  newdat,
  aes(x = day, y = pred, color = angle_type)
) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180)
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Day",
    y = "Predicted angle (degrees)",
    color = "Angle type",
    title = "Model 3 predicted trajectories"
  )

print(p_time_model)


## =========================
## 7. SAVE PLOTS
## =========================

# create folder if it doesn't exist
dir.create("plots/model3", recursive = TRUE, showWarnings = FALSE)

## ---- 1. Genotype plot ----
ggsave(
  filename = "plots/model3/angle_genotype.png",
  plot = p_angle_genotype,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/angle_genotype.pdf",
  plot = p_angle_genotype,
  width = 7,
  height = 5,
  bg = "white"
)

## ---- 2. Time slices ----
ggsave(
  filename = "plots/model3/angle_time_slices.png",
  plot = p_time_slices,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/angle_time_slices.pdf",
  plot = p_time_slices,
  width = 7,
  height = 5,
  bg = "white"
)

## ---- 3. Trajectories ----
ggsave(
  filename = "plots/model3/angle_trajectories.png",
  plot = p_time_model,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/angle_trajectories.pdf",
  plot = p_time_model,
  width = 8,
  height = 5,
  bg = "white"
)

## =========================
## 12. OPTIONAL: COMPACT CONSOLE SUMMARY
## =========================
cat("\n====================================\n")
cat("MODEL 3: MAIN TAKE-HOMES\n")
cat("====================================\n")
cat("\n1. Type III ANOVA:\n")
print(anova_m3)

cat("\n2. Overall angle means:\n")
print(emm_angle_overall)

cat("\n3. Overall angle contrasts:\n")
print(pairs_angle_overall)

cat("\n4. Genotype effects within angle type:\n")
print(emm_genotype_by_angle)

cat("\n5. Time slopes by angle type:\n")
print(angle_slopes)

cat("\n6. Angle estimates at representative times:\n")
print(emm_angle_by_time_df)

## =========================
## REAL-DATA EARLY / MID / LATE BINS
## =========================

early_cut <- time_points_raw["Early"]
late_cut  <- time_points_raw["Late"]

angles_binned <- angles_long %>%
  mutate(
    time_bin = case_when(
      day <= early_cut ~ "Early",
      day >= late_cut  ~ "Late",
      TRUE             ~ "Mid"
    ),
    time_bin = factor(time_bin, levels = c("Early", "Mid", "Late"))
  )

## Quick summary table
angle_bin_summary <- angles_binned %>%
  group_by(time_bin, angle_type) %>%
  summarise(
    n = n(),
    mean = mean(angle, na.rm = TRUE),
    sd = sd(angle, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

print(angle_bin_summary)

p_time_boxes <- ggplot(
  angles_binned,
  aes(x = time_bin, y = angle, fill = angle_type)
) +
  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.15
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.2,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Angle (degrees)",
    title = "Angle distributions across Early, Mid, and Late stages"
  )

print(p_time_boxes)


p_time_violins <- ggplot(
  angles_binned,
  aes(x = time_bin, y = angle, fill = angle_type)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    color = NA
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.alpha = 0.15
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.2,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Angle (degrees)",
    title = "Angle distributions across Early, Mid, and Late stages"
  )

print(p_time_violins)

p_time_bars <- ggplot(
  angle_bin_summary,
  aes(x = time_bin, y = mean, fill = angle_type)
) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Mean angle (degrees)",
    title = "Mean angles across Early, Mid, and Late stages"
  )

print(p_time_bars)


ggsave(
  "plots/model3/angle_time_boxes.png",
  p_time_boxes,
  width = 9,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  "plots/model3/angle_time_boxes.pdf",
  p_time_boxes,
  width = 9,
  height = 4,
  bg = "white"
)

ggsave(
  "plots/model3/angle_time_violins.png",
  p_time_violins,
  width = 9,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  "plots/model3/angle_time_violins.pdf",
  p_time_violins,
  width = 9,
  height = 4,
  bg = "white"
)

ggsave(
  "plots/model3/angle_time_bars.png",
  p_time_bars,
  width = 9,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  "plots/model3/angle_time_bars.pdf",
  p_time_bars,
  width = 9,
  height = 4,
  bg = "white"
)

## =========================
## REAL-DATA GENOTYPE SUMMARY
## =========================
angle_gen_summary <- angles_long %>%
  group_by(genotype, angle_type) %>%
  summarise(
    n = n(),
    mean = mean(angle, na.rm = TRUE),
    sd = sd(angle, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

print(angle_gen_summary)

p_gen_boxes_real <- ggplot(
  angles_long,
  aes(x = genotype, y = angle, fill = angle_type)
) +
  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.15
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.2,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Genotype",
    y = "Angle (degrees)",
    title = "Angle distributions across genotypes"
  )

print(p_gen_boxes_real)

p_gen_violins_real <- ggplot(
  angles_long,
  aes(x = genotype, y = angle, fill = angle_type)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.45,
    color = NA
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.alpha = 0.15
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.2,
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Genotype",
    y = "Angle (degrees)",
    title = "Angle distributions across genotypes"
  )

print(p_gen_violins_real)

p_gen_bars_real <- ggplot(
  angle_gen_summary,
  aes(x = genotype, y = mean, fill = angle_type)
) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Genotype",
    y = "Mean angle (degrees)",
    title = "Mean angles across genotypes"
  )

print(p_gen_bars_real)

## =========================
## RAW-DATA EARLY / MID / LATE SUMMARY
## =========================

early_cut <- time_points_raw["Early"]
late_cut  <- time_points_raw["Late"]

angles_binned <- angles_long %>%
  mutate(
    time_bin = case_when(
      day <= early_cut ~ "Early",
      day >= late_cut  ~ "Late",
      TRUE             ~ "Mid"
    ),
    time_bin = factor(time_bin, levels = c("Early", "Mid", "Late"))
  )

angle_time_raw <- angles_binned %>%
  group_by(time_bin, angle_type) %>%
  summarise(
    n = n(),
    mean_angle = mean(angle, na.rm = TRUE),
    sd_angle = sd(angle, na.rm = TRUE),
    se_angle = sd_angle / sqrt(n),
    lower_95 = mean_angle - 1.96 * se_angle,
    upper_95 = mean_angle + 1.96 * se_angle,
    .groups = "drop"
  )

p_time_slices_raw <- ggplot(
  angle_time_raw,
  aes(x = time_bin, y = mean_angle, color = angle_type, group = angle_type)
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0.15,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Observed mean angle (degrees)",
    color = "Angle type",
    title = "Angle changes through time"
  )

print(p_time_slices_raw)

## =========================
## RAW-DATA GENOTYPE POINT-RANGE FIGURE
## =========================

## Summary table from raw data
angle_gen_raw <- angles_long %>%
  group_by(genotype, angle_type) %>%
  summarise(
    n = n(),
    mean_angle = mean(angle, na.rm = TRUE),
    sd_angle = sd(angle, na.rm = TRUE),
    se_angle = sd_angle / sqrt(n),
    lower_95 = mean_angle - 1.96 * se_angle,
    upper_95 = mean_angle + 1.96 * se_angle,
    .groups = "drop"
  )

print(angle_gen_raw)

## Plot
p_angle_genotype_raw <- ggplot(
  angle_gen_raw,
  aes(x = genotype, y = mean_angle, color = angle_type, group = angle_type)
) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 2.8
  ) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.4),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_color_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(30, 90, by = 10),
    limits = c(27, 93),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Genotype",
    y = "Observed mean angle (degrees)",
    color = "Angle type",
    title = "Angle variation across genotypes"
  )

print(p_angle_genotype_raw)

## Save
dir.create("plots/model3", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "plots/model3/angle_genotype_raw.png",
  plot = p_angle_genotype_raw,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/angle_genotype_raw.pdf",
  plot = p_angle_genotype_raw,
  width = 8,
  height = 5,
  bg = "white"
)
