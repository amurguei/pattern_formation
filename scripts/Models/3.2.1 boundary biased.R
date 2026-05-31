## ============================================================
## Boundary-biased positioning: full self-contained analysis
## ============================================================

# ----------------------------
# 0. Working directory
# ----------------------------
setwd("/Users/amalia/Documents/GitHub/pattern_formation")

# ----------------------------
# 1. Packages
# ----------------------------
library(tidyverse)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

# ----------------------------
# 2. Read data
# ----------------------------
boundary_df <- read_csv("inputs/PERIMETER_NUBBIN_DISTANCE_RATIO_REPORT.csv")
slope_df    <- read_csv("inputs/NUBBIN_DISTANCE_LINEAR_FIT_BY_GENOTYPE_DAY45_TO_127.csv")

# ----------------------------
# 3. Clean / prepare
# ----------------------------
boundary_df <- boundary_df %>%
  rename(
    genotype   = genotype,
    nubbin_id  = replica,
    day        = day,
    d_perimeter = dist_perimeter,
    d_nubbin    = dist_nubbin,
    ratio       = nubbin_perimeter_ratio
  ) %>%
  mutate(
    genotype  = factor(genotype),
    nubbin_id = interaction(genotype, nubbin_id, drop = TRUE),
    day_sc    = as.numeric(scale(day))
  )

# Optional quick checks
cat("\n--- Data structure ---\n")
print(glimpse(boundary_df))

cat("\n--- Missing values ---\n")
print(colSums(is.na(boundary_df)))

cat("\n--- N observations ---\n")
print(nrow(boundary_df))

cat("\n--- N nubbins ---\n")
print(n_distinct(boundary_df$nubbin_id))

# ----------------------------
# 4. Main model 1: ratio
# ----------------------------
m_ratio <- lmer(
  ratio ~ genotype + day_sc + (1 | nubbin_id),
  data = boundary_df
)

anova_ratio <- anova(m_ratio)
emm_ratio   <- emmeans(m_ratio, pairwise ~ genotype)

cat("\n==============================\n")
cat("RATIO MODEL\n")
cat("==============================\n")
print(summary(m_ratio))
print(anova_ratio)
print(emm_ratio)

# Clean genotype table for manuscript / supplement
ratio_table <- as.data.frame(emm_ratio$emmeans) %>%
  transmute(
    Genotype = as.character(genotype),
    `Estimated mean ratio (d_nubbin / d_perimeter)` = round(emmean, 2),
    SE = round(SE, 3),
    lower_CL = round(asymp.LCL, 2),
    upper_CL = round(asymp.UCL, 2)
  )

print(ratio_table)

cat("\n--- Ratio estimated means table ---\n")
print(ratio_table)

write_csv(ratio_table, "boundary_ratio_estimated_means.csv")

# ----------------------------
# 5. Main model 2: distance to nubbin
# ----------------------------
# This is the core temporal model for outward displacement
m_nubbin <- lmer(
  d_nubbin ~ genotype * day_sc + (1 | nubbin_id),
  data = boundary_df
)

anova_nubbin <- anova(m_nubbin)

cat("\n==============================\n")
cat("NUBBIN DISTANCE MODEL\n")
cat("==============================\n")
print(summary(m_nubbin))
print(anova_nubbin)

# Optional marginal means by genotype
emm_nubbin <- emmeans(m_nubbin, ~ genotype)
cat("\n--- Nubbin distance estimated means by genotype ---\n")
print(emm_nubbin)

# ----------------------------
# 6. Main model 3: distance to perimeter
# ----------------------------
m_perimeter <- lmer(
  d_perimeter ~ genotype + day_sc + (1 | nubbin_id),
  data = boundary_df
)

anova_perimeter <- anova(m_perimeter)
emm_perimeter   <- emmeans(m_perimeter, pairwise ~ genotype)

cat("\n==============================\n")
cat("PERIMETER DISTANCE MODEL\n")
cat("==============================\n")
print(summary(m_perimeter))
print(anova_perimeter)
print(emm_perimeter)

# ----------------------------
# 7. Compare temporal slopes more explicitly
# ----------------------------
# If you want genotype-specific slopes from the mixed model:
slopes_nubbin <- emtrends(m_nubbin, ~ genotype, var = "day_sc")

cat("\n==============================\n")
cat("GENOTYPE-SPECIFIC SLOPES (scaled day)\n")
cat("==============================\n")
print(slopes_nubbin)

# Pairwise slope differences, if needed
slopes_nubbin_pairs <- pairs(slopes_nubbin)
cat("\n--- Pairwise comparisons of nubbin-distance slopes ---\n")
print(slopes_nubbin_pairs)

# ----------------------------
# 8. Variability comparison
# ----------------------------
# This supports the statement that perimeter distance is more constrained
variability_by_genotype <- boundary_df %>%
  group_by(genotype) %>%
  summarise(
    mean_nubbin = mean(d_nubbin, na.rm = TRUE),
    sd_nubbin   = sd(d_nubbin, na.rm = TRUE),
    mean_perim  = mean(d_perimeter, na.rm = TRUE),
    sd_perim    = sd(d_perimeter, na.rm = TRUE),
    mean_ratio  = mean(ratio, na.rm = TRUE),
    sd_ratio    = sd(ratio, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n==============================\n")
cat("VARIABILITY BY GENOTYPE\n")
cat("==============================\n")
print(variability_by_genotype)

write_csv(variability_by_genotype, "boundary_variability_by_genotype.csv")

# Overall SD comparison
overall_variability <- tibble(
  metric = c("d_nubbin", "d_perimeter", "ratio"),
  mean   = c(mean(boundary_df$d_nubbin, na.rm = TRUE),
             mean(boundary_df$d_perimeter, na.rm = TRUE),
             mean(boundary_df$ratio, na.rm = TRUE)),
  sd     = c(sd(boundary_df$d_nubbin, na.rm = TRUE),
             sd(boundary_df$d_perimeter, na.rm = TRUE),
             sd(boundary_df$ratio, na.rm = TRUE))
)

cat("\n==============================\n")
cat("OVERALL VARIABILITY\n")
cat("==============================\n")
print(overall_variability)

write_csv(overall_variability, "boundary_overall_variability.csv")

# ----------------------------
# 9. Optional: inspect the precomputed slope file
# ----------------------------
cat("\n==============================\n")
cat("PRECOMPUTED LINEAR SLOPES (day 45-127)\n")
cat("==============================\n")
print(slope_df)

# ----------------------------
# 10. Simple figures
# ----------------------------

# Figure A: ratio by genotype
p_ratio <- ggplot(boundary_df, aes(x = genotype, y = ratio)) +
  geom_boxplot(fill = "cadetblue3", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 1) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = "red",
    color = "black"
  ) +
  labs(
    x = "Genotype",
    y = "d_nubbin / d_perimeter"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Figure B: mean nubbin / perimeter distances through time
time_summary <- boundary_df %>%
  group_by(day) %>%
  summarise(
    mean_nubbin = mean(d_nubbin, na.rm = TRUE),
    sd_nubbin   = sd(d_nubbin, na.rm = TRUE),
    mean_perim  = mean(d_perimeter, na.rm = TRUE),
    sd_perim    = sd(d_perimeter, na.rm = TRUE),
    .groups = "drop"
  )

p_time <- ggplot(time_summary, aes(x = day)) +
  geom_line(aes(y = mean_nubbin, color = "Nubbin"), linewidth = 0.9) +
  geom_point(aes(y = mean_nubbin, color = "Nubbin"), size = 2) +
  geom_errorbar(
    aes(ymin = mean_nubbin - sd_nubbin, ymax = mean_nubbin + sd_nubbin, color = "Nubbin"),
    width = 1.5, alpha = 0.7
  ) +
  geom_line(aes(y = mean_perim, color = "Perimeter"), linewidth = 0.9) +
  geom_point(aes(y = mean_perim, color = "Perimeter"), size = 2) +
  geom_errorbar(
    aes(ymin = mean_perim - sd_perim, ymax = mean_perim + sd_perim, color = "Perimeter"),
    width = 1.5, alpha = 0.7
  ) +
  scale_color_manual(values = c("Nubbin" = "blue", "Perimeter" = "black")) +
  labs(
    x = "Day",
    y = "Mean distance (mm)",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

print(p_ratio)
print(p_time)

# Save figures
ggsave("boundary_ratio_boxplot.png", p_ratio, width = 150, height = 110, units = "mm", dpi = 600, bg = "white")
ggsave("boundary_time_dynamics.png", p_time, width = 170, height = 110, units = "mm", dpi = 600, bg = "white")

# ----------------------------
# 11. Compact stats table for manuscript
# ----------------------------
results_table <- tibble(
  model = c("ratio", "ratio", "nubbin_distance", "nubbin_distance", "nubbin_distance", "perimeter_distance", "perimeter_distance"),
  term  = c("genotype", "day_sc", "genotype", "day_sc", "genotype:day_sc", "genotype", "day_sc"),
  F_value = c(
    anova_ratio["genotype", "F value"],
    anova_ratio["day_sc", "F value"],
    anova_nubbin["genotype", "F value"],
    anova_nubbin["day_sc", "F value"],
    anova_nubbin["genotype:day_sc", "F value"],
    anova_perimeter["genotype", "F value"],
    anova_perimeter["day_sc", "F value"]
  ),
  p_value = c(
    anova_ratio["genotype", "Pr(>F)"],
    anova_ratio["day_sc", "Pr(>F)"],
    anova_nubbin["genotype", "Pr(>F)"],
    anova_nubbin["day_sc", "Pr(>F)"],
    anova_nubbin["genotype:day_sc", "Pr(>F)"],
    anova_perimeter["genotype", "Pr(>F)"],
    anova_perimeter["day_sc", "Pr(>F)"]
  )
)

print(results_table)
write_csv(results_table, "boundary_positioning_model_results.csv")

# ----------------------------
# Binary variable
# ----------------------------
boundary_df <- boundary_df %>%
  mutate(
    closer_perimeter = ifelse(d_perimeter < d_nubbin, 1, 0)
  )

# Quick check
table(boundary_df$closer_perimeter)

# ----------------------------
# GLMM (logistic mixed model)
# ----------------------------
library(lme4)

m_binary <- glmer(
  closer_perimeter ~ genotype + day_sc + (1 | nubbin_id),
  data = boundary_df,
  family = binomial
)

summary(m_binary)

# ----------------------------
# ANOVA-style test
# ----------------------------
library(car)
Anova(m_binary, type = 3)

# ----------------------------
# Estimated probabilities
# ----------------------------
library(emmeans)

emm_binary <- emmeans(m_binary, ~ genotype, type = "response")

print(emm_binary)


# ---------------------------------
# Predicted probability over time
# ---------------------------------
library(emmeans)

emm_time <- emmeans(m_binary, ~ day_sc, type = "response")

emm_time_df <- as.data.frame(emm_time)

ggplot(emm_time_df, aes(x = day_sc, y = prob)) +
  geom_line(linewidth = 1.2, color = "black") +
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL),
              alpha = 0.2, fill = "grey70") +
  labs(
    x = "Scaled time",
    y = "Probability of forming closer to colony perimeter"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black")
  )


ggplot(boundary_df, aes(x = genotype, y = closer_perimeter)) +
  geom_boxplot(fill = "cadetblue3", outlier.shape = NA) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = "red",
    color = "black"
  ) +
  labs(
    x = "Colony",
    y = "Proportion of polyps closer to perimeter"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black")
  )


# ---------------------------------
# Binary summary by genotype
# ---------------------------------

boundary_df <- boundary_df %>%
  mutate(
    closer_perimeter = ifelse(d_perimeter < d_nubbin, 1, 0),
    closer_nubbin    = 1 - closer_perimeter
  )

prob_df <- boundary_df %>%
  group_by(genotype) %>%
  summarise(
    p_perimeter = mean(closer_perimeter, na.rm = TRUE),
    sd_perimeter = sd(closer_perimeter, na.rm = TRUE),
    p_nubbin = mean(closer_nubbin, na.rm = TRUE),
    sd_nubbin = sd(closer_nubbin, na.rm = TRUE),
    n = n(),
    se_perimeter = sd_perimeter / sqrt(n),
    se_nubbin = sd_nubbin / sqrt(n),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(p_perimeter, p_nubbin, se_perimeter, se_nubbin),
    names_to = c(".value", "position"),
    names_pattern = "(p|se)_(perimeter|nubbin)"
  ) %>%
  mutate(
    position = factor(position, levels = c("perimeter", "nubbin"),
                      labels = c("Closer to perimeter", "Closer to nubbin"))
  )

print(prob_df)

# ---------------------------------
# Plot
# ---------------------------------

ggplot(prob_df, aes(x = genotype, y = p, group = position)) +
  geom_point(aes(shape = position, color = position),
             position = position_dodge(width = 0.35),
             size = 3) +
  geom_errorbar(
    aes(ymin = p - se, ymax = p + se, color = position),
    position = position_dodge(width = 0.35),
    width = 0.15,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c("Closer to perimeter" = "black",
               "Closer to nubbin" = "blue")
  ) +
  scale_shape_manual(
    values = c("Closer to perimeter" = 16,
               "Closer to nubbin" = 17)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "Colony",
    y = "Probability of new polyp position",
    color = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggplot(prob_df, aes(x = genotype, y = p, fill = position)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_errorbar(
    aes(ymin = p - se, ymax = p + se),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  scale_fill_manual(
    values = c("Closer to perimeter" = "black",
               "Closer to nubbin" = "blue")
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "Colony",
    y = "Probability of new polyp position",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )


# ---------------------------------
# Prepare data
# ---------------------------------
boundary_df <- boundary_df %>%
  mutate(
    closer_perimeter = ifelse(d_perimeter < d_nubbin, 1, 0),
    closer_nubbin    = 1 - closer_perimeter
  )

# ---------------------------------
# Summary + binomial test
# ---------------------------------
prob_df <- boundary_df %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    successes = sum(closer_perimeter),
    p_perimeter = mean(closer_perimeter),
    p_nubbin    = mean(closer_nubbin),
    se = sd(closer_perimeter) / sqrt(n),
    p_value = binom.test(successes, n, p = 0.5)$p.value,
    signif = ifelse(p_value < 0.05, "*", "ns"),
    .groups = "drop"
  )

# Long format
prob_long <- prob_df %>%
  select(genotype, p_perimeter, p_nubbin, se) %>%
  pivot_longer(
    cols = c(p_perimeter, p_nubbin),
    names_to = "position",
    values_to = "p"
  ) %>%
  mutate(
    position = factor(position,
                      levels = c("p_perimeter", "p_nubbin"),
                      labels = c("Closer to perimeter", "Closer to nubbin"))
  )

# ---------------------------------
# Plot
# ---------------------------------
ggplot(prob_long, aes(x = genotype, y = p, fill = position)) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             linewidth = 0.6,
             color = "black") +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.65,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = p - se, ymax = p + se),
    position = position_dodge(width = 0.7),
    width = 0.15
  ) +
  geom_text(
    data = prob_df,
    aes(x = genotype, y = 1.03, label = signif),
    inherit.aes = FALSE,
    size = 5
  ) +
  scale_fill_manual(
    values = c("Closer to perimeter" = "#CC79A7",
               "Closer to nubbin" = "#0072B2")
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    x = "Colony",
    y = "Probability of new polyp position",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )


m_ratio_day <- lmer(
  ratio ~ genotype + day + (1 | nubbin_id),
  data = boundary_df
)

m_nubbin_day <- lmer(
  d_nubbin ~ genotype * day + (1 | nubbin_id),
  data = boundary_df
)

m_perimeter_day <- lmer(
  d_perimeter ~ genotype + day + (1 | nubbin_id),
  data = boundary_df
)

m_binary_day <- glmer(
  closer_perimeter ~ genotype + day + (1 | nubbin_id),
  data = boundary_df,
  family = binomial
)

anova(m_ratio_day)
anova(m_nubbin_day)
anova(m_perimeter_day)
car::Anova(m_binary_day, type = 3)


# ---------------------------------
# Refit models using real day
# ---------------------------------
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

boundary_df <- boundary_df %>%
  mutate(
    closer_perimeter = ifelse(d_perimeter < d_nubbin, 1, 0),
    closer_nubbin    = 1 - closer_perimeter
  )

m_ratio_day <- lmer(
  ratio ~ genotype + day + (1 | nubbin_id),
  data = boundary_df
)

m_nubbin_day <- lmer(
  d_nubbin ~ genotype * day + (1 | nubbin_id),
  data = boundary_df
)

m_perimeter_day <- lmer(
  d_perimeter ~ genotype + day + (1 | nubbin_id),
  data = boundary_df
)

m_binary_day <- glmer(
  closer_perimeter ~ genotype + day + (1 | nubbin_id),
  data = boundary_df,
  family = binomial
)

cat("\n=== RATIO MODEL ===\n")
print(anova(m_ratio_day))
print(emmeans(m_ratio_day, ~ genotype))

cat("\n=== NUBBIN DISTANCE MODEL ===\n")
print(anova(m_nubbin_day))
print(emtrends(m_nubbin_day, ~ genotype, var = "day"))

cat("\n=== PERIMETER DISTANCE MODEL ===\n")
print(anova(m_perimeter_day))
print(emmeans(m_perimeter_day, ~ genotype))

cat("\n=== BINARY MODEL ===\n")
print(summary(m_binary_day))
print(car::Anova(m_binary_day, type = 3))
print(emmeans(m_binary_day, ~ genotype, type = "response"))





# =========================================================
# Boundary-biased positioning composite figure
# =========================================================

library(tidyverse)
library(patchwork)
library(grid)

# ---------------------------------------------------------
# 0. Global theme
# ---------------------------------------------------------
theme_set(
  theme_minimal(base_size = 16) +
    theme(
      axis.title = element_text(size = 18, face = "bold"),
      axis.text  = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text  = element_text(size = 15),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.x = unit(0.45, "cm"),
      panel.grid.minor = element_blank()
    )
)

# ---------------------------------------------------------
# 1. Prepare data
# ---------------------------------------------------------
geno_levels <- c("1","2","3","4","6","7","8","9","10")

boundary_df <- boundary_df %>%
  mutate(
    genotype = factor(as.character(genotype), levels = geno_levels),
    closer_perimeter = ifelse(d_perimeter < d_nubbin, 1, 0),
    closer_nubbin    = 1 - closer_perimeter
  )

# ---------------------------------------------------------
# 2. Panel A: binary probability by genotype
# ---------------------------------------------------------
prob_df <- boundary_df %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    successes = sum(closer_perimeter),
    p_perimeter = mean(closer_perimeter),
    p_nubbin    = mean(closer_nubbin),
    se = sd(closer_perimeter) / sqrt(n),
    p_value = binom.test(successes, n, p = 0.5)$p.value,
    star_y = pmax(p_perimeter + se, p_nubbin + se) + 0.018,
    .groups = "drop"
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
  geom_point(
    data = prob_df,
    aes(x = genotype, y = star_y),
    inherit.aes = FALSE,
    shape = 8,
    size = 2,      # <- your requested asterisk size
    stroke = 0.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Closer to perimeter" = "#CC79A7",
      "Closer to nubbin"    = "#56B4E9"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1.06),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Colony",
    y = "Probability",
    fill = NULL,
    title = "A. Probability of boundary-biased polyp formation"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# ---------------------------------------------------------
# 3. Panel B: ratio by genotype
# ---------------------------------------------------------
pB <- ggplot(boundary_df, aes(x = genotype, y = ratio)) +
  geom_boxplot(
    fill = "#56B4E9",
    alpha = 0.75,
    outlier.shape = NA,
    width = 0.65,
    color = "black"
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.18,
    size = 1.2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = "red",
    color = "black",
    stroke = 0.5
  ) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(
    x = "Colony",
    y = expression(d[nubbin] / d[perimeter]),
    title = "B. Relative positioning of new polyps"
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# ---------------------------------------------------------
# 4. Panel C: time dynamics
# ---------------------------------------------------------
time_summary <- boundary_df %>%
  group_by(day) %>%
  summarise(
    mean_nubbin = mean(d_nubbin, na.rm = TRUE),
    sd_nubbin   = sd(d_nubbin, na.rm = TRUE),
    mean_perim  = mean(d_perimeter, na.rm = TRUE),
    sd_perim    = sd(d_perimeter, na.rm = TRUE),
    .groups = "drop"
  )

pC <- ggplot(time_summary, aes(x = day)) +
  geom_line(aes(y = mean_nubbin, color = "Nubbin"), linewidth = 1) +
  geom_point(aes(y = mean_nubbin, color = "Nubbin"), size = 2.8) +
  geom_errorbar(
    aes(
      ymin = mean_nubbin - sd_nubbin,
      ymax = mean_nubbin + sd_nubbin,
      color = "Nubbin"
    ),
    width = 1.2,
    linewidth = 0.6,
    alpha = 0.8
  ) +
  geom_line(aes(y = mean_perim, color = "Perimeter"), linewidth = 1) +
  geom_point(aes(y = mean_perim, color = "Perimeter"), size = 2.8) +
  geom_errorbar(
    aes(
      ymin = mean_perim - sd_perim,
      ymax = mean_perim + sd_perim,
      color = "Perimeter"
    ),
    width = 1.2,
    linewidth = 0.6,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Nubbin" = "#56B4E9",
      "Perimeter" = "black"
    )
  ) +
  labs(
    x = "Day",
    y = "Mean distance (mm)",
    color = NULL,
    title = "C. Temporal dynamics of boundary-biased positioning"
  ) +
  theme(
    legend.position = "top"
  )

# ---------------------------------------------------------
# 5. Composite
# ---------------------------------------------------------
figure_boundary <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1.05))

print(figure_boundary)

# ---------------------------------------------------------
# 6. Save composite
# ---------------------------------------------------------
ggsave(
  filename = "boundary_positioning_composite_final.png",
  plot = figure_boundary,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "boundary_positioning_composite_final.pdf",
  plot = figure_boundary,
  width = 11,
  height = 8.5,
  bg = "white"
)

# ---------------------------------------------------------
# 7. Save individual panels
# ---------------------------------------------------------
ggsave(
  filename = "boundary_panel_A_binary.png",
  plot = pA,
  width = 7,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "boundary_panel_A_binary.pdf",
  plot = pA,
  width = 7,
  height = 4.8,
  bg = "white"
)

ggsave(
  filename = "boundary_panel_B_ratio.png",
  plot = pB,
  width = 7,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "boundary_panel_B_ratio.pdf",
  plot = pB,
  width = 7,
  height = 4.8,
  bg = "white"
)

ggsave(
  filename = "boundary_panel_C_time.png",
  plot = pC,
  width = 8.5,
  height = 5.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "boundary_panel_C_time.pdf",
  plot = pC,
  width = 8.5,
  height = 5.2,
  bg = "white"
)




pA <- pA + labs(title = "Probability of boundary-biased polyp formation")

pB <- pB + labs(title = "Relative positioning of new polyps")

pC <- pC + labs(title = "Temporal dynamics of boundary-biased positioning")


ggsave("Fig_boundary_A_binary.png", pA, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("Fig_boundary_A_binary.pdf", pA, width = 7, height = 5, bg = "white")

ggsave("Fig_boundary_B_ratio.png", pB, width = 7, height = 5, dpi = 600, bg = "white")
ggsave("Fig_boundary_B_ratio.pdf", pB, width = 7, height = 5, bg = "white")

ggsave("Fig_boundary_C_time.png", pC, width = 8, height = 5.5, dpi = 600, bg = "white")
ggsave("Fig_boundary_C_time.pdf", pC, width = 8, height = 5.5, bg = "white")
