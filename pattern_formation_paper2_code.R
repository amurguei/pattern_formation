############################################################
## TRIANGLE / NEW POLYP ANALYSIS
## Template for mixed models, diagnostics, figures, and
## geometric invariants
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
  "performance",
  "DHARMa",
  "broom.mixed",
  "patchwork",
  "car",
  "ggdist"
)

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(DHARMa)
library(broom.mixed)
library(patchwork)
library(car)
library(ggdist)

## =========================
## 1. LOAD DATA
## =========================
file_path <- "triangles_full_data.xlsx"

tri <- read_excel(file_path) %>%
  clean_names()

glimpse(tri)
summary(tri)

## Expected columns after clean_names():
## dist_to_nn1, dist_to_nn2, angle_main, angle_nn1, angle_nn2,
## genotype, replica, day

## =========================
## 2. BASIC CLEANING
## =========================
tri <- tri %>%
  mutate(
    genotype = factor(genotype),
    replica  = factor(replica),
    nubbin_id = interaction(genotype, replica, drop = TRUE, sep = "_"),
    day = as.numeric(day)
  )

## Check missing values
colSums(is.na(tri))

## Optional: remove impossible values if needed
tri <- tri %>%
  filter(
    dist_to_nn1 > 0,
    dist_to_nn2 > 0,
    angle_main > 0, angle_main < 180,
    angle_nn1 > 0,  angle_nn1 < 180,
    angle_nn2 > 0,  angle_nn2 < 180
  )

## Check triangle angle sum
tri <- tri %>%
  mutate(
    angle_sum = angle_main + angle_nn1 + angle_nn2,
    angle_sum_error = angle_sum - 180
  )

summary(tri$angle_sum)
summary(tri$angle_sum_error)

## =========================
## 3. DERIVED VARIABLES / GEOMETRIC INVARIANTS
## =========================
tri <- tri %>%
  mutate(
    dist_ratio = dist_to_nn2 / dist_to_nn1,
    dist_diff  = dist_to_nn2 - dist_to_nn1,
    
    ## order statistics for angles
    angle_min = pmin(angle_main, angle_nn1, angle_nn2),
    angle_max = pmax(angle_main, angle_nn1, angle_nn2),
    angle_mid = angle_sum - angle_min - angle_max,
    
    ## equilateral deviation
    dev_main_60 = angle_main - 60,
    dev_nn1_60  = angle_nn1 - 60,
    dev_nn2_60  = angle_nn2 - 60,
    
    ## absolute deviations from equilateral expectation
    abs_dev_main_60 = abs(angle_main - 60),
    abs_dev_nn1_60  = abs(angle_nn1 - 60),
    abs_dev_nn2_60  = abs(angle_nn2 - 60),
    
    ## total equilateral deviation
    total_abs_dev_60 = abs_dev_main_60 + abs_dev_nn1_60 + abs_dev_nn2_60,
    
    ## deviation from 70-70-40 motif after sorting
    dev_704070 = abs(angle_max - 70) + abs(angle_mid - 70) + abs(angle_min - 40),
    
    ## a simple shape asymmetry metric
    angle_asymmetry = angle_max - angle_min
  )

summary(tri$dist_ratio)
summary(tri$total_abs_dev_60)
summary(tri$dev_704070)
summary(tri$angle_asymmetry)

tri <- tri %>%
  mutate(
    dev_806040 =
      abs(angle_max - 80) +
      abs(angle_mid - 60) +
      abs(angle_min - 40)
  )

summary(tri$dev_806040)
## =========================
## 4. LONG FORMAT FOR ANGLES
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
    )
  )

glimpse(angles_long)

## =========================
## 5. EXPLORATORY PLOTS
## =========================

### =========================
## 5. EXPLORATORY PLOTS
## =========================

## =========================
## 5. EXPLORATORY PLOTS
## =========================

library(scales)

## Create folder for plots
if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

## Okabe-Ito palette
angle_cols <- c(
  "Main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

dist_cols <- c(
  "NN1" = "#0072B2",
  "NN2" = "#D55E00"
)

## Force factor levels to match palette names exactly
angles_long <- angles_long %>%
  mutate(angle_type = factor(angle_type, levels = c("Main", "NN1", "NN2")))

## Distances long format
dist_long <- tri %>%
  select(genotype, replica, nubbin_id, day, dist_to_nn1, dist_to_nn2) %>%
  pivot_longer(
    cols = c(dist_to_nn1, dist_to_nn2),
    names_to = "dist_type",
    values_to = "distance"
  ) %>%
  mutate(
    dist_type = factor(
      dist_type,
      levels = c("dist_to_nn1", "dist_to_nn2"),
      labels = c("NN1", "NN2")
    )
  ) %>%
  mutate(dist_type = factor(dist_type, levels = c("NN1", "NN2")))

## Sanity checks
stopifnot(all(levels(angles_long$angle_type) %in% names(angle_cols)))
stopifnot(all(levels(dist_long$dist_type) %in% names(dist_cols)))

## Mean values for density plot vertical lines
angle_means <- angles_long %>%
  group_by(angle_type) %>%
  summarise(mean_angle = mean(angle, na.rm = TRUE), .groups = "drop")

## 5a. Distributions of angles by type
## Keep this one without extra y-axis tinkering
p_angle_density <- ggplot(
  angles_long,
  aes(x = angle, fill = angle_type, color = angle_type)
) +
  geom_density(alpha = 0.30, linewidth = 0.9) +
  geom_vline(
    data = angle_means,
    aes(xintercept = mean_angle, color = angle_type),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  facet_wrap(~ angle_type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = angle_cols) +
  scale_color_manual(values = angle_cols) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Angle (degrees)",
    y = "Density",
    title = "Angle distributions by type"
  )

## 5b. Violin + boxplot of angle by type
p_angle_violin <- ggplot(
  angles_long,
  aes(x = angle_type, y = angle, fill = angle_type)
) +
  geom_violin(alpha = 0.5, trim = FALSE, color = NA) +
  geom_boxplot(
    width = 0.15,
    outlier.alpha = 0.2,
    fill = "white",
    color = "black"
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
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Angle (degrees)",
    title = "Angles differ by triangle vertex"
  )

## 5c. Distances
p_dist <- ggplot(
  dist_long,
  aes(x = dist_type, y = distance, fill = dist_type)
) +
  geom_violin(alpha = 0.5, trim = FALSE, color = NA) +
  geom_boxplot(
    width = 0.15,
    outlier.alpha = 0.2,
    fill = "white",
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.8,
    color = "black"
  ) +
  scale_fill_manual(values = dist_cols) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  ) +
  labs(
    x = NULL,
    y = "Distance (mm)",
    title = "Distances to nearest neighbors"
  )

## 5d. Time trends by nubbin
p_time <- ggplot(
  angles_long,
  aes(x = day, y = angle, group = nubbin_id, color = angle_type)
) +
  geom_line(alpha = 0.12, linewidth = 0.4) +
  stat_summary(
    aes(group = angle_type),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  stat_summary(
    aes(group = angle_type),
    fun = mean,
    geom = "point",
    size = 2
  ) +
  facet_wrap(~ angle_type, ncol = 1) +
  scale_color_manual(values = angle_cols) +
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
    x = "Day",
    y = "Angle (degrees)",
    title = "Time trends in angles"
  )

## Combined plots
p_angles_combined <- p_angle_density / p_angle_violin
p_dist_time_combined <- p_dist / p_time

## Print to viewer
print(p_angle_density)
print(p_angle_violin)
print(p_dist)
print(p_time)
print(p_angles_combined)
print(p_dist_time_combined)

## =========================
## SAVE PLOTS
## =========================

## Individual PNGs
ggsave(
  filename = "plots/p_angle_density.png",
  plot = p_angle_density,
  width = 7,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/p_angle_violin.png",
  plot = p_angle_violin,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist.png",
  plot = p_dist,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/p_time.png",
  plot = p_time,
  width = 7,
  height = 9,
  dpi = 600,
  bg = "white"
)

## Combined PNGs
ggsave(
  filename = "plots/p_angles_combined.png",
  plot = p_angles_combined,
  width = 10,
  height = 12,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist_time_combined.png",
  plot = p_dist_time_combined,
  width = 10,
  height = 12,
  dpi = 600,
  bg = "white"
)

## Individual PDFs
ggsave(
  filename = "plots/p_angle_density.pdf",
  plot = p_angle_density,
  width = 7,
  height = 9,
  bg = "white"
)

ggsave(
  filename = "plots/p_angle_violin.pdf",
  plot = p_angle_violin,
  width = 6,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist.pdf",
  plot = p_dist,
  width = 6,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "plots/p_time.pdf",
  plot = p_time,
  width = 7,
  height = 9,
  bg = "white"
)

## Combined PDFs
ggsave(
  filename = "plots/p_angles_combined.pdf",
  plot = p_angles_combined,
  width = 10,
  height = 12,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist_time_combined.pdf",
  plot = p_dist_time_combined,
  width = 10,
  height = 12,
  bg = "white"
)
## Combined plots
ggsave(
  filename = "plots/p_angles_combined.png",
  plot = p_angles_combined,
  width = 10,
  height = 12,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist_time_combined.png",
  plot = p_dist_time_combined,
  width = 10,
  height = 12,
  dpi = 600,
  bg = "white"
)

## Optional PDFs too (great for papers)
ggsave(
  filename = "plots/p_angle_density.pdf",
  plot = p_angle_density,
  width = 7,
  height = 9,
  bg = "white"
)

ggsave(
  filename = "plots/p_angle_violin.pdf",
  plot = p_angle_violin,
  width = 6,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist.pdf",
  plot = p_dist,
  width = 6,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "plots/p_time.pdf",
  plot = p_time,
  width = 7,
  height = 9,
  bg = "white"
)

ggsave(
  filename = "plots/p_angles_combined.pdf",
  plot = p_angles_combined,
  width = 10,
  height = 12,
  bg = "white"
)

ggsave(
  filename = "plots/p_dist_time_combined.pdf",
  plot = p_dist_time_combined,
  width = 10,
  height = 12,
  bg = "white"
)

## =========================
## 6. CORE MODEL 1:
##    ARE THE 3 ANGLES DIFFERENT?
## =========================
m_angle0 <- lmer(
  angle ~ angle_type + (1 | nubbin_id),
  data = angles_long
)

summary(m_angle0)
anova(m_angle0)

emmeans(m_angle0, pairwise ~ angle_type)

## =========================
## 7. CORE MODEL 2:
##    DO ANGLE DIFFERENCES DEPEND ON GENOTYPE?
## =========================
m_angle1 <- lmer(
  angle ~ angle_type * genotype + (1 | nubbin_id),
  data = angles_long
)

summary(m_angle1)
anova(m_angle1)
emmeans(m_angle1, pairwise ~ angle_type)
emmeans(m_angle1, pairwise ~ angle_type | genotype)
emmeans(m_angle1, pairwise ~ genotype | angle_type)

#CORE MODEL 2 (MODIFIED):
  ##    DO ANGLE DIFFERENCES DEPEND ON GENOTYPE?
  ## =========================

## Packages needed
library(car)
library(emmeans)
library(broom)
library(broom.mixed)

## Make sure factors are set correctly
angles_long <- angles_long %>%
  mutate(
    angle_type = factor(angle_type, levels = c("Main", "NN1", "NN2")),
    genotype   = factor(genotype)
  )

## Set sum-to-zero contrasts for Type III tests
## (important when using car::Anova(type = 3))
options(contrasts = c("contr.sum", "contr.poly"))

## -------------------------
## 7.1 MAIN MODEL: fixed-effects model
## -------------------------
m_angle1_lm <- lm(
  angle ~ angle_type * genotype,
  data = angles_long
)

summary(m_angle1_lm)

## Type III ANOVA
anova_angle1_lm <- car::Anova(m_angle1_lm, type = 3)
anova_angle1_lm

## -------------------------
## 7.2 ESTIMATED MARGINAL MEANS
## -------------------------

## Overall angle_type effect
emm_angle_type <- emmeans(m_angle1_lm, ~ angle_type)
pairs_angle_type <- pairs(emm_angle_type, adjust = "tukey")

emm_angle_type
pairs_angle_type

## Angle type within each genotype
emm_angle_by_genotype <- emmeans(m_angle1_lm, ~ angle_type | genotype)
pairs_angle_by_genotype <- pairs(emm_angle_by_genotype, adjust = "tukey")

emm_angle_by_genotype
pairs_angle_by_genotype

## Genotype within each angle type
emm_genotype_by_angle <- emmeans(m_angle1_lm, ~ genotype | angle_type)
pairs_genotype_by_angle <- pairs(emm_genotype_by_angle, adjust = "tukey")

emm_genotype_by_angle
pairs_genotype_by_angle

## -------------------------
## 7.3 SENSITIVITY CHECK: mixed model
## -------------------------
m_angle1_lmer <- lmer(
  angle ~ angle_type * genotype + (1 | nubbin_id),
  data = angles_long
)

summary(m_angle1_lmer)
anova(m_angle1_lmer)

## Optional: compare fitted values from lm and lmer
cor(
  fitted(m_angle1_lm),
  fitted(m_angle1_lmer)
)

## -------------------------
## 7.4 TIDY OUTPUTS FOR EXPORT
## -------------------------

## Fixed effects table
tidy_m_angle1_lm <- broom::tidy(m_angle1_lm, conf.int = TRUE)
write.csv(
  tidy_m_angle1_lm,
  file = "model_angle1_lm_coefficients.csv",
  row.names = FALSE
)

## Type III ANOVA table
anova_angle1_lm_df <- as.data.frame(anova_angle1_lm)
anova_angle1_lm_df$term <- rownames(anova_angle1_lm_df)
rownames(anova_angle1_lm_df) <- NULL

write.csv(
  anova_angle1_lm_df,
  file = "model_angle1_lm_type3_anova.csv",
  row.names = FALSE
)

## Estimated marginal means
emm_angle_type_df <- as.data.frame(emm_angle_type)
emm_angle_by_genotype_df <- as.data.frame(emm_angle_by_genotype)
emm_genotype_by_angle_df <- as.data.frame(emm_genotype_by_angle)

write.csv(
  emm_angle_type_df,
  file = "model_angle1_lm_emmeans_angle_type.csv",
  row.names = FALSE
)

write.csv(
  emm_angle_by_genotype_df,
  file = "model_angle1_lm_emmeans_angle_by_genotype.csv",
  row.names = FALSE
)

write.csv(
  emm_genotype_by_angle_df,
  file = "model_angle1_lm_emmeans_genotype_by_angle.csv",
  row.names = FALSE
)

## Pairwise comparisons
pairs_angle_type_df <- as.data.frame(pairs_angle_type)
pairs_angle_by_genotype_df <- as.data.frame(pairs_angle_by_genotype)
pairs_genotype_by_angle_df <- as.data.frame(pairs_genotype_by_angle)

write.csv(
  pairs_angle_type_df,
  file = "model_angle1_lm_pairs_angle_type.csv",
  row.names = FALSE
)

write.csv(
  pairs_angle_by_genotype_df,
  file = "model_angle1_lm_pairs_angle_by_genotype.csv",
  row.names = FALSE
)

write.csv(
  pairs_genotype_by_angle_df,
  file = "model_angle1_lm_pairs_genotype_by_angle.csv",
  row.names = FALSE
)

## -------------------------
## 7.5 DIAGNOSTICS
## -------------------------
par(mfrow = c(2, 2))
plot(m_angle1_lm)
par(mfrow = c(1, 1))

## Optional normality check of residuals
shapiro.test(sample(resid(m_angle1_lm), size = min(5000, length(resid(m_angle1_lm)))))

## -------------------------
## 7.6 PLOT: genotype-specific angle means
## -------------------------

## Use emmeans output for cleaner figure
emm_plot_df <- as.data.frame(emm_angle_by_genotype)

p_angle_genotype <- ggplot(
  emm_plot_df,
  aes(x = genotype, y = emmean, color = angle_type, group = angle_type)
) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 2.5
  ) +
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
    y = "Estimated marginal mean angle (°)",
    color = "Angle type",
    title = "Genotype-specific triangle geometry"
  )

print(p_angle_genotype)

ggsave(
  filename = "plots/p_angle_genotype_emmeans.png",
  plot = p_angle_genotype,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/p_angle_genotype_emmeans.pdf",
  plot = p_angle_genotype,
  width = 8,
  height = 5,
  bg = "white"
)

## -------------------------
## 7.7 OPTIONAL: compact summaries in console
## -------------------------
cat("\n=== TYPE III ANOVA (LM) ===\n")
print(anova_angle1_lm)

cat("\n=== EMMEANS: angle_type ===\n")
print(emm_angle_type)

cat("\n=== PAIRWISE: angle_type ===\n")
print(pairs_angle_type)

cat("\n=== EMMEANS: angle_type within genotype ===\n")
print(emm_angle_by_genotype)

cat("\n=== PAIRWISE: angle_type within genotype ===\n")
print(pairs_angle_by_genotype)

cat("\n=== EMMEANS: genotype within angle_type ===\n")
print(emm_genotype_by_angle)

cat("\n=== PAIRWISE: genotype within angle_type ===\n")
print(pairs_genotype_by_angle)

## =========================
## 8. CORE MODEL 3:
##    ADD TIME
## =========================
m_angle2 <- lmer(
  angle ~ angle_type * genotype + scale(day) + (1 | nubbin_id),
  data = angles_long
)

summary(m_angle2)
anova(m_angle2)

## Optional: allow nubbin-specific slopes through time
m_angle2b <- lmer(
  angle ~ angle_type * genotype + scale(day) + (1 + scale(day) | nubbin_id),
  data = angles_long,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(m_angle2b)
anova(m_angle2, m_angle2b)

## Optional: time interaction with angle type
m_angle3 <- lmer(
  angle ~ angle_type * genotype + angle_type * scale(day) + (1 | nubbin_id),
  data = angles_long
)

summary(m_angle3)
anova(m_angle3)

library(emmeans)

emtrends(m_angle3, ~ angle_type, var = "day")

## =========================
## 9. DISTANCE MODELS
## =========================
m_dist0 <- lmer(
  distance ~ dist_type + (1 | nubbin_id),
  data = dist_long
)

summary(m_dist0)
anova(m_dist0)
emmeans(m_dist0, pairwise ~ dist_type)

m_dist1 <- lmer(
  distance ~ dist_type * genotype + scale(day) + (1 | nubbin_id),
  data = dist_long
)

summary(m_dist1)
anova(m_dist1)

## =========================
## 10. DISTANCE RATIO MODEL
## =========================
m_ratio <- lmer(
  dist_ratio ~ genotype + scale(day) + (1 | nubbin_id),
  data = tri
)

summary(m_ratio)
anova(m_ratio)
emmeans(m_ratio, pairwise ~ genotype)

## =========================
## 11. TESTING AGAINST HYPOTHESES
## =========================

## 11a. Is the angle sum close to 180? (sanity check)
t.test(tri$angle_sum, mu = 180)

## 11b. Are mean angles equal to 60?
t.test(tri$angle_main, mu = 60)
t.test(tri$angle_nn1,  mu = 60)
t.test(tri$angle_nn2,  mu = 60)

## Better: mixed model on deviations from 60 in long format
angles_long_60 <- angles_long %>%
  mutate(dev_from_60 = angle - 60)

m_60 <- lmer(
  dev_from_60 ~ angle_type + genotype + scale(day) + (1 | nubbin_id),
  data = angles_long_60
)

summary(m_60)
anova(m_60)
emmeans(m_60, pairwise ~ angle_type)

## 11c. Compare fit to 60-60-60 vs 70-70-40 motif
## Lower values = closer fit
m_compare_shape <- lmer(
  cbind(total_abs_dev_60, dev_704070)[,1] ~ genotype + scale(day) + (1 | nubbin_id),
  data = tri
)

## That formulation is clunky for direct comparison, so instead:
shape_long <- tri %>%
  select(genotype, nubbin_id, day, total_abs_dev_60, dev_704070) %>%
  pivot_longer(
    cols = c(total_abs_dev_60, dev_704070),
    names_to = "shape_model",
    values_to = "deviation"
  ) %>%
  mutate(
    shape_model = factor(
      shape_model,
      levels = c("total_abs_dev_60", "dev_704070"),
      labels = c("60-60-60", "70-70-40")
    )
  )

m_shape <- lmer(
  deviation ~ shape_model + genotype + scale(day) + (1 | nubbin_id),
  data = shape_long
)

summary(m_shape)
anova(m_shape)
emmeans(m_shape, pairwise ~ shape_model)

tri %>%
  mutate(
    better_fit = case_when(
      dev_704070 < dev_806040 ~ "70-70-40",
      dev_806040 < dev_704070 ~ "80-60-40",
      TRUE ~ "tie"
    )
  ) %>%
  group_by(day, better_fit) %>%
  summarise(n = n(), .groups = "drop")


shape_long2 <- tri %>%
  select(genotype, nubbin_id, day, dev_704070, dev_806040) %>%
  pivot_longer(
    cols = c(dev_704070, dev_806040),
    names_to = "shape_model",
    values_to = "deviation"
  ) %>%
  mutate(
    shape_model = factor(
      shape_model,
      levels = c("dev_704070", "dev_806040"),
      labels = c("70-70-40", "80-60-40")
    )
  )
shape_long2 <- tri %>%
  select(genotype, nubbin_id, day, dev_704070, dev_806040) %>%
  pivot_longer(
    cols = c(dev_704070, dev_806040),
    names_to = "shape_model",
    values_to = "deviation"
  ) %>%
  mutate(
    shape_model = factor(
      shape_model,
      levels = c("dev_704070", "dev_806040"),
      labels = c("70-70-40", "80-60-40")
    )
  )

m_shape_time <- lmer(
  deviation ~ shape_model * scale(day) + genotype + (1 | nubbin_id),
  data = shape_long2
)

anova(m_shape_time)
emmeans(m_shape_time, ~ shape_model)

## =========================
## 12. MODEL DIAGNOSTICS
## =========================

## Residual diagnostics for main angle model
check_model(m_angle2)

## Basic residual plots
par(mfrow = c(2, 2))
plot(m_angle2)
qqnorm(resid(m_angle2))
qqline(resid(m_angle2))
hist(resid(m_angle2), breaks = 30)

## DHARMa-style simulated residuals
sim_m_angle2 <- simulateResiduals(m_angle2)
plot(sim_m_angle2)
testUniformity(sim_m_angle2)
testDispersion(sim_m_angle2)
testOutliers(sim_m_angle2)

## Check singularity / overfitting
check_singularity(m_angle2b)
check_collinearity(m_angle2)


## =========================
## 13. FIGURES FOR THE PAPER
## =========================

## Figure 1: angle distributions + means/CI
fig1 <- ggplot(angles_long, aes(x = angle_type, y = angle, fill = angle_type)) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.25, point_colour = NA) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.15) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.08) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Angle (degrees)", title = "Triangle angles associated with new polyp insertion")

## Figure 2: genotype-specific angle means
fig2 <- ggplot(angles_long, aes(x = genotype, y = angle, color = angle_type, group = angle_type)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.4), size = 2) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",
               position = position_dodge(width = 0.4), width = 0.2) +
  theme_classic(base_size = 12) +
  labs(x = "Genotype", y = "Angle (degrees)", title = "Genotype-specific geometry")

## Figure 3: time trends
fig3 <- ggplot(angles_long, aes(x = day, y = angle, color = angle_type)) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ genotype) +
  theme_classic(base_size = 12) +
  labs(x = "Day", y = "Angle (degrees)", title = "Temporal dynamics of triangle angles")

## Figure 4: distance ratio and shape invariant
fig4a <- ggplot(tri, aes(x = genotype, y = dist_ratio)) +
  geom_violin(fill = "grey85", alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.15) +
  theme_classic(base_size = 12) +
  labs(x = "Genotype", y = "NN2 / NN1 distance ratio", title = "Distance ratio across genotypes")

fig4b <- ggplot(tri, aes(x = genotype, y = dev_704070)) +
  geom_violin(fill = "grey85", alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.15) +
  theme_classic(base_size = 12) +
  labs(x = "Genotype", y = "Deviation from 70-70-40", title = "Closeness to 70-70-40 motif")

fig4 <- fig4a + fig4b

fig1
fig2
fig3
fig4

## Save figures
ggsave("fig1_angles_main.png", fig1, width = 6, height = 5, dpi = 300)
ggsave("fig2_genotypes.png", fig2, width = 7, height = 5, dpi = 300)
ggsave("fig3_time_trends.png", fig3, width = 9, height = 6, dpi = 300)
ggsave("fig4_invariants.png", fig4, width = 10, height = 5, dpi = 300)

## =========================
## 14. OPTIONAL: NUBBIN-LEVEL SUMMARIES
## =========================
nubbin_summary <- tri %>%
  group_by(genotype, nubbin_id) %>%
  summarise(
    n_triangles = n(),
    mean_main = mean(angle_main, na.rm = TRUE),
    mean_nn1  = mean(angle_nn1,  na.rm = TRUE),
    mean_nn2  = mean(angle_nn2,  na.rm = TRUE),
    mean_ratio = mean(dist_ratio, na.rm = TRUE),
    mean_dev_704070 = mean(dev_704070, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(nubbin_summary, "nubbin_summary.csv")

## =========================
## 15. EXPORT MODEL OUTPUTS
## =========================
tidy_m_angle2 <- tidy(m_angle2, effects = "fixed", conf.int = TRUE)
tidy_m_dist1  <- tidy(m_dist1, effects = "fixed", conf.int = TRUE)
tidy_m_ratio  <- tidy(m_ratio, effects = "fixed", conf.int = TRUE)
tidy_m_shape  <- tidy(m_shape, effects = "fixed", conf.int = TRUE)

write_csv(tidy_m_angle2, "model_angle2_fixed_effects.csv")
write_csv(tidy_m_dist1,  "model_dist1_fixed_effects.csv")
write_csv(tidy_m_ratio,  "model_ratio_fixed_effects.csv")
write_csv(tidy_m_shape,  "model_shape_fixed_effects.csv")


#Checking early vs. late

early_5360 <- angles_long %>%
  filter(day %in% c(53, 60))

late_127 <- angles_long %>%
  filter(day >= 127)

summary_5360 <- early_5360 %>%
  group_by(angle_type) %>%
  summarise(
    n    = n(),
    mean = mean(angle, na.rm = TRUE),
    sd   = sd(angle, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(period = "Days 53–60")

summary_127 <- late_127 %>%
  group_by(angle_type) %>%
  summarise(
    n    = n(),
    mean = mean(angle, na.rm = TRUE),
    sd   = sd(angle, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(period = "Day ≥ 127")

summary_window <- bind_rows(summary_5360, summary_127)
summary_window

library(tidyr)

summary_window %>%
  select(period, angle_type, mean, sd) %>%
  pivot_wider(
    names_from = period,
    values_from = c(mean, sd)
  )

q1 <- quantile(angles_long$day, 0.25)
q3 <- quantile(angles_long$day, 0.75)

early_Q1 <- angles_long %>%
  filter(day <= q1)

late_Q4 <- angles_long %>%
  filter(day >= q3)

summary_Q1 <- early_Q1 %>%
  group_by(angle_type) %>%
  summarise(
    n    = n(),
    mean = mean(angle, na.rm = TRUE),
    sd   = sd(angle, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(period = "First quartile")

summary_Q4 <- late_Q4 %>%
  group_by(angle_type) %>%
  summarise(
    n    = n(),
    mean = mean(angle, na.rm = TRUE),
    sd   = sd(angle, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(period = "Last quartile")

summary_quartiles <- bind_rows(summary_Q1, summary_Q4)
summary_quartiles

summary_quartiles %>%
  select(period, angle_type, mean, sd) %>%
  pivot_wider(
    names_from = period,
    values_from = c(mean, sd)
  )

angles_long %>%
  mutate(period = case_when(
    day %in% c(53, 60) ~ "Early",
    day >= 127 ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(angle_type) %>%
  summarise(
    t_test = list(t.test(angle ~ period)),
    .groups = "drop"
  )


## =========================================================
## 12. EARLY vs LATE ANGLE COMPARISONS
## Two approaches:
##   A) Biological windows: days 53–60 vs day >= 127
##   B) Quartiles: first quartile vs last quartile
## =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

## Make sure angle_type has the right order
angles_long <- angles_long %>%
  mutate(
    angle_type = factor(angle_type, levels = c("Main", "NN1", "NN2"))
  )

## Optional: create plots folder
if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

## Use same colors as before
angle_cols <- c(
  "Main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

period_cols <- c(
  "Early" = "#0072B2",
  "Late"  = "#D55E00",
  "Early (Q1)" = "#0072B2",
  "Late (Q4)"  = "#D55E00"
)

## =========================================================
## A) BIOLOGICAL WINDOWS
## Early = days 53 and 60
## Late  = day >= 127
## =========================================================

early_late_window <- angles_long %>%
  mutate(period = case_when(
    day %in% c(53, 60) ~ "Early",
    day >= 127 ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  mutate(period = factor(period, levels = c("Early", "Late")))

## -------------------------
## A1. Summary table
## -------------------------
summary_window <- early_late_window %>%
  group_by(period, angle_type) %>%
  summarise(
    n    = n(),
    mean = mean(angle, na.rm = TRUE),
    sd   = sd(angle, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

print(summary_window)

summary_window_wide <- summary_window %>%
  select(period, angle_type, mean, sd, se, n) %>%
  pivot_wider(
    names_from = period,
    values_from = c(mean, sd, se, n)
  )

print(summary_window_wide)

write.csv(summary_window, "window_summary_long.csv", row.names = FALSE)
write.csv(summary_window_wide, "window_summary_wide.csv", row.names = FALSE)

## -------------------------
## A2. Mixed model
## -------------------------
m_window <- lmer(
  angle ~ angle_type * period + (1 | nubbin_id),
  data = early_late_window
)

cat("\n=== WINDOW MODEL ANOVA ===\n")
print(anova(m_window))

cat("\n=== WINDOW MODEL SUMMARY ===\n")
print(summary(m_window))

## -------------------------
## A3. emmeans
## -------------------------
emm_window_angle_by_period <- emmeans(m_window, ~ angle_type | period)
pairs_window_angle_by_period <- pairs(emm_window_angle_by_period, adjust = "tukey")

emm_window_period_by_angle <- emmeans(m_window, ~ period | angle_type)
pairs_window_period_by_angle <- pairs(emm_window_period_by_angle, adjust = "tukey")

cat("\n=== WINDOW: angle_type within period ===\n")
print(emm_window_angle_by_period)
print(pairs_window_angle_by_period)

cat("\n=== WINDOW: period within angle_type ===\n")
print(emm_window_period_by_angle)
print(pairs_window_period_by_angle)

write.csv(as.data.frame(emm_window_angle_by_period),
          "window_emmeans_angle_by_period.csv", row.names = FALSE)
write.csv(as.data.frame(pairs_window_angle_by_period),
          "window_pairs_angle_by_period.csv", row.names = FALSE)
write.csv(as.data.frame(emm_window_period_by_angle),
          "window_emmeans_period_by_angle.csv", row.names = FALSE)
write.csv(as.data.frame(pairs_window_period_by_angle),
          "window_pairs_period_by_angle.csv", row.names = FALSE)

## -------------------------
## A4. Violin plot
## -------------------------
p_window_violin <- ggplot(
  early_late_window,
  aes(x = period, y = angle, fill = period)
) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = period_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Angle (degrees)",
    title = "Angle distributions: Early vs Late (days 53–60 vs day ≥ 127)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_window_violin)

ggsave(
  filename = "plots/window_violin_angles.png",
  plot = p_window_violin,
  width = 10,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/window_violin_angles.pdf",
  plot = p_window_violin,
  width = 10,
  height = 4,
  bg = "white"
)

## -------------------------
## A5. Mean ± SE plot
## -------------------------
p_window_means <- ggplot(
  summary_window,
  aes(x = angle_type, y = mean, fill = period)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_fill_manual(values = period_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Angle type",
    y = "Mean angle (degrees)",
    title = "Mean angle values: Early vs Late (days 53–60 vs day ≥ 127)"
  ) +
  theme_classic(base_size = 12)

print(p_window_means)

ggsave(
  filename = "plots/window_mean_angles.png",
  plot = p_window_means,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/window_mean_angles.pdf",
  plot = p_window_means,
  width = 7,
  height = 5,
  bg = "white"
)

## =========================================================
## B) QUARTILE APPROACH
## Early = first quartile of day
## Late  = last quartile of day
## =========================================================

q1 <- quantile(angles_long$day, 0.25, na.rm = TRUE)
q3 <- quantile(angles_long$day, 0.75, na.rm = TRUE)

cat("\n=== QUARTILE CUTPOINTS ===\n")
print(q1)
print(q3)

quartile_data <- angles_long %>%
  mutate(period = case_when(
    day <= q1 ~ "Early (Q1)",
    day >= q3 ~ "Late (Q4)",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  mutate(period = factor(period, levels = c("Early (Q1)", "Late (Q4)")))

## -------------------------
## B1. Summary table
## -------------------------
summary_quartiles <- quartile_data %>%
  group_by(period, angle_type) %>%
  summarise(
    n    = n(),
    mean = mean(angle, na.rm = TRUE),
    sd   = sd(angle, na.rm = TRUE),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

print(summary_quartiles)

summary_quartiles_wide <- summary_quartiles %>%
  select(period, angle_type, mean, sd, se, n) %>%
  pivot_wider(
    names_from = period,
    values_from = c(mean, sd, se, n)
  )

print(summary_quartiles_wide)

write.csv(summary_quartiles, "quartile_summary_long.csv", row.names = FALSE)
write.csv(summary_quartiles_wide, "quartile_summary_wide.csv", row.names = FALSE)

## -------------------------
## B2. Mixed model
## -------------------------
m_quartile <- lmer(
  angle ~ angle_type * period + (1 | nubbin_id),
  data = quartile_data
)

cat("\n=== QUARTILE MODEL ANOVA ===\n")
print(anova(m_quartile))

cat("\n=== QUARTILE MODEL SUMMARY ===\n")
print(summary(m_quartile))

## -------------------------
## B3. emmeans
## -------------------------
emm_quartile_angle_by_period <- emmeans(m_quartile, ~ angle_type | period)
pairs_quartile_angle_by_period <- pairs(emm_quartile_angle_by_period, adjust = "tukey")

emm_quartile_period_by_angle <- emmeans(m_quartile, ~ period | angle_type)
pairs_quartile_period_by_angle <- pairs(emm_quartile_period_by_angle, adjust = "tukey")

cat("\n=== QUARTILE: angle_type within period ===\n")
print(emm_quartile_angle_by_period)
print(pairs_quartile_angle_by_period)

cat("\n=== QUARTILE: period within angle_type ===\n")
print(emm_quartile_period_by_angle)
print(pairs_quartile_period_by_angle)

write.csv(as.data.frame(emm_quartile_angle_by_period),
          "quartile_emmeans_angle_by_period.csv", row.names = FALSE)
write.csv(as.data.frame(pairs_quartile_angle_by_period),
          "quartile_pairs_angle_by_period.csv", row.names = FALSE)
write.csv(as.data.frame(emm_quartile_period_by_angle),
          "quartile_emmeans_period_by_angle.csv", row.names = FALSE)
write.csv(as.data.frame(pairs_quartile_period_by_angle),
          "quartile_pairs_period_by_angle.csv", row.names = FALSE)

## -------------------------
## B4. Violin plot
## -------------------------
p_quartile_violin <- ggplot(
  quartile_data,
  aes(x = period, y = angle, fill = period)
) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = period_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Angle (degrees)",
    title = "Angle distributions: First quartile vs Last quartile"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_quartile_violin)

ggsave(
  filename = "plots/quartile_violin_angles.png",
  plot = p_quartile_violin,
  width = 10,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/quartile_violin_angles.pdf",
  plot = p_quartile_violin,
  width = 10,
  height = 4,
  bg = "white"
)

## -------------------------
## B5. Mean ± SE plot
## -------------------------
p_quartile_means <- ggplot(
  summary_quartiles,
  aes(x = angle_type, y = mean, fill = period)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_fill_manual(values = period_cols) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Angle type",
    y = "Mean angle (degrees)",
    title = "Mean angle values: First quartile vs Last quartile"
  ) +
  theme_classic(base_size = 12)

print(p_quartile_means)

ggsave(
  filename = "plots/quartile_mean_angles.png",
  plot = p_quartile_means,
  width = 7,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/quartile_mean_angles.pdf",
  plot = p_quartile_means,
  width = 7,
  height = 5,
  bg = "white"
)



## =========================================================
## 13. NN2 / NN1 DISTANCE RATIO PANEL
## A) Overall distribution
## B) By genotype
## C) Through time
## =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

## Create folder for plots
if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

## Use tri directly
ratio_data <- tri %>%
  mutate(
    genotype = factor(genotype),
    day = as.numeric(day)
  ) %>%
  filter(
    is.finite(dist_ratio),
    !is.na(dist_ratio),
    dist_ratio > 0
  )

## Optional: inspect the ratio quickly
summary(ratio_data$dist_ratio)

## Overall mean
ratio_mean <- mean(ratio_data$dist_ratio, na.rm = TRUE)

## ---------------------------------------------------------
## PANEL A — Overall distribution
## ---------------------------------------------------------
p_ratio_A <- ggplot(ratio_data, aes(x = "All data", y = dist_ratio)) +
  geom_violin(fill = "#56B4E9", alpha = 0.7, color = NA, trim = FALSE) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "NN2 / NN1 distance ratio",
    title = "A. Overall distribution"
  ) +
  theme_classic(base_size = 12)

## ---------------------------------------------------------
## PANEL B — By genotype
## ---------------------------------------------------------
p_ratio_B <- ggplot(ratio_data, aes(x = genotype, y = dist_ratio)) +
  geom_violin(fill = "#009E73", alpha = 0.7, color = NA, trim = FALSE) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = "Genotype",
    y = "NN2 / NN1 distance ratio",
    title = "B. By genotype"
  ) +
  theme_classic(base_size = 12)

## ---------------------------------------------------------
## PANEL C — Through time
## ---------------------------------------------------------
p_ratio_C <- ggplot(ratio_data, aes(x = day, y = dist_ratio)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  labs(
    x = "Day",
    y = "NN2 / NN1 distance ratio",
    title = "C. Through time"
  ) +
  theme_classic(base_size = 12)

## ---------------------------------------------------------
## COMBINED PANEL
## ---------------------------------------------------------
p_ratio_panel <- p_ratio_A | p_ratio_B | p_ratio_C

print(p_ratio_panel)

## Save high-res versions
ggsave(
  filename = "plots/nn_ratio_panel.png",
  plot = p_ratio_panel,
  width = 15,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel.pdf",
  plot = p_ratio_panel,
  width = 15,
  height = 4.8,
  bg = "white"
)

## =========================================================
## OPTIONAL: EARLY vs LATE RATIO VIOLIN
## =========================================================

ratio_window <- ratio_data %>%
  mutate(period = case_when(
    day %in% c(53, 60) ~ "Early",
    day >= 127 ~ "Late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  mutate(period = factor(period, levels = c("Early", "Late")))

p_ratio_window <- ggplot(ratio_window, aes(x = period, y = dist_ratio, fill = period)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "NN2 / NN1 distance ratio",
    title = "NN distance ratio: Early vs Late"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

print(p_ratio_window)

ggsave(
  filename = "plots/nn_ratio_early_late.png",
  plot = p_ratio_window,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_early_late.pdf",
  plot = p_ratio_window,
  width = 6,
  height = 5,
  bg = "white"
)

## =========================================================
## OPTIONAL: SUMMARY TABLES
## =========================================================

ratio_summary_genotype <- ratio_data %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    mean = mean(dist_ratio, na.rm = TRUE),
    sd = sd(dist_ratio, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

ratio_summary_overall <- ratio_data %>%
  summarise(
    n = n(),
    mean = mean(dist_ratio, na.rm = TRUE),
    sd = sd(dist_ratio, na.rm = TRUE),
    se = sd / sqrt(n)
  )

write.csv(ratio_summary_genotype, "ratio_summary_by_genotype.csv", row.names = FALSE)
write.csv(ratio_summary_overall, "ratio_summary_overall.csv", row.names = FALSE)
############################################################
## END
############################################################


## =========================================================
## 13. NN2 / NN1 DISTANCE RATIO PANEL
## A) Overall distribution
## B) By genotype
## C) Through time
## =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

## Create folder for plots
if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

## ---------------------------------------------------------
## PREP DATA
## ---------------------------------------------------------
ratio_data <- tri %>%
  mutate(
    genotype = factor(genotype),
    day = as.numeric(day)
  ) %>%
  filter(
    is.finite(dist_ratio),
    !is.na(dist_ratio),
    dist_ratio > 0
  )

## Optional quick inspection
summary(ratio_data$dist_ratio)

## Set y zoom for all panels
ratio_ylim <- c(0.95, 2.5)

## ---------------------------------------------------------
## PANEL A — Overall distribution
## ---------------------------------------------------------
p_ratio_A <- ggplot(ratio_data, aes(x = "All data", y = dist_ratio)) +
  geom_violin(
    fill = "#56B4E9",
    alpha = 0.7,
    color = NA,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = NULL,
    y = "NN2 / NN1 distance ratio",
    title = "A. Overall distribution"
  ) +
  theme_classic(base_size = 12)

## ---------------------------------------------------------
## PANEL B — By genotype
## ---------------------------------------------------------
p_ratio_B <- ggplot(ratio_data, aes(x = genotype, y = dist_ratio)) +
  geom_violin(
    fill = "#009E73",
    alpha = 0.7,
    color = NA,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2,
    color = "black"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = "Genotype",
    y = "NN2 / NN1 distance ratio",
    title = "B. By genotype"
  ) +
  theme_classic(base_size = 12)

## ---------------------------------------------------------
## PANEL C — Through time
## ---------------------------------------------------------
p_ratio_C <- ggplot(ratio_data, aes(x = day, y = dist_ratio)) +
  geom_point(
    alpha = 0.15,
    size = 0.8
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = "Day",
    y = "NN2 / NN1 distance ratio",
    title = "C. Through time"
  ) +
  theme_classic(base_size = 12)

## ---------------------------------------------------------
## COMBINED PANEL
## ---------------------------------------------------------
p_ratio_panel <- p_ratio_A | p_ratio_B | p_ratio_C

print(p_ratio_panel)

## Save high-res versions
ggsave(
  filename = "plots/nn_ratio_panel.png",
  plot = p_ratio_panel,
  width = 15,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel.pdf",
  plot = p_ratio_panel,
  width = 15,
  height = 4.8,
  bg = "white"
)

## =========================================================
## OPTIONAL: EARLY vs LATE RATIO VIOLIN
## =========================================================

ratio_window <- ratio_data %>%
  mutate(
    period = case_when(
      day %in% c(53, 60) ~ "Early",
      day >= 127 ~ "Late",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period)) %>%
  mutate(period = factor(period, levels = c("Early", "Late")))

p_ratio_window <- ggplot(ratio_window, aes(x = period, y = dist_ratio, fill = period)) +
  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    color = NA
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = NULL,
    y = "NN2 / NN1 distance ratio",
    title = "NN distance ratio: Early vs Late"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

print(p_ratio_window)

ggsave(
  filename = "plots/nn_ratio_early_late.png",
  plot = p_ratio_window,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_early_late.pdf",
  plot = p_ratio_window,
  width = 6,
  height = 5,
  bg = "white"
)

## =========================================================
## OPTIONAL: SUMMARY TABLES
## =========================================================

ratio_summary_genotype <- ratio_data %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    mean = mean(dist_ratio, na.rm = TRUE),
    sd = sd(dist_ratio, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

ratio_summary_overall <- ratio_data %>%
  summarise(
    n = n(),
    mean = mean(dist_ratio, na.rm = TRUE),
    sd = sd(dist_ratio, na.rm = TRUE),
    se = sd / sqrt(n)
  )

write.csv(ratio_summary_genotype, "ratio_summary_by_genotype.csv", row.names = FALSE)
write.csv(ratio_summary_overall, "ratio_summary_overall.csv", row.names = FALSE)


## =========================================================
## 13. NN2 / NN1 DISTANCE RATIO PANEL
## 2 + 1 LAYOUT
## A) Overall distribution
## B) By genotype
## C) Through time
## =========================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

## Create folder for plots
if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

## ---------------------------------------------------------
## PREP DATA
## ---------------------------------------------------------
ratio_data <- tri %>%
  mutate(
    genotype = factor(genotype),
    day = as.numeric(day)
  ) %>%
  filter(
    is.finite(dist_ratio),
    !is.na(dist_ratio),
    dist_ratio > 0
  )

summary(ratio_data$dist_ratio)

## Slightly tighter zoom for visual balance
ratio_ylim <- c(0.95, 2.3)

## Shared theme tweaks
theme_ratio <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "plain"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(8, 10, 8, 10)
  )

## ---------------------------------------------------------
## PANEL A — Overall distribution
## ---------------------------------------------------------
p_ratio_A <- ggplot(ratio_data, aes(x = "All data", y = dist_ratio)) +
  geom_violin(
    fill = "#56B4E9",
    alpha = 0.7,
    color = NA,
    trim = TRUE,
    width = 0.8
  ) +
  geom_boxplot(
    width = 0.10,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.4,
    color = "black"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = NULL,
    y = "NN2 / NN1 distance ratio",
    title = "A. Overall distribution"
  ) +
  theme_ratio

## ---------------------------------------------------------
## PANEL B — By genotype
## ---------------------------------------------------------
p_ratio_B <- ggplot(ratio_data, aes(x = genotype, y = dist_ratio)) +
  geom_violin(
    fill = "#009E73",
    alpha = 0.7,
    color = NA,
    trim = TRUE,
    width = 0.85
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.1,
    color = "black"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = "Genotype",
    y = "NN2 / NN1 distance ratio",
    title = "B. By genotype"
  ) +
  theme_ratio

## ---------------------------------------------------------
## PANEL C — Through time
## ---------------------------------------------------------
p_ratio_C <- ggplot(ratio_data, aes(x = day, y = dist_ratio)) +
  geom_point(
    alpha = 0.12,
    size = 0.8,
    color = "black"
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 1,
    color = "#2C5EFF"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = "Day",
    y = "NN2 / NN1 distance ratio",
    title = "C. Through time"
  ) +
  theme_ratio

## ---------------------------------------------------------
## 2 + 1 COMBINED LAYOUT
## ---------------------------------------------------------
p_ratio_panel <- (p_ratio_A | p_ratio_B) / p_ratio_C +
  plot_layout(heights = c(1, 1.05))

print(p_ratio_panel)

## Save
ggsave(
  filename = "plots/nn_ratio_panel_2plus1.png",
  plot = p_ratio_panel,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_2plus1.pdf",
  plot = p_ratio_panel,
  width = 11,
  height = 8.5,
  bg = "white"
)
