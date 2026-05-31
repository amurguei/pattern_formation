############################################################
## FIGURE 10: Early, mid, late angle plot
## Bars = observed mean angle
## Error bars = SE
## Bins = Early <= mean(day) - 1 SD
##        Mid   > mean(day) - 1 SD and < mean(day) + 1 SD
##        Late  >= mean(day) + 1 SD
############################################################

## =========================
## 0. PACKAGES
## =========================
packages <- c("tidyverse", "readxl", "janitor", "scales")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(scales)

## =========================
## 1. LOAD + CLEAN DATA
## =========================

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

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
    angle_main > 0, angle_main < 180,
    angle_nn1  > 0, angle_nn1  < 180,
    angle_nn2  > 0, angle_nn2  < 180
  )

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
    )
  )

## =========================
## 3. DEFINE EARLY / MID / LATE BINS
## =========================

day_mean <- mean(angles_long$day, na.rm = TRUE)
day_sd   <- sd(angles_long$day, na.rm = TRUE)

early_cut <- day_mean - day_sd
late_cut  <- day_mean + day_sd

## Print cutoffs so you can report exact values
cat("\nEarly cutoff:", round(early_cut, 1), "days\n")
cat("Mid range:", round(early_cut, 1), "to", round(late_cut, 1), "days\n")
cat("Late cutoff:", round(late_cut, 1), "days\n")

angles_binned <- angles_long %>%
  mutate(
    time_bin = case_when(
      day <= early_cut ~ "Early",
      day >= late_cut  ~ "Late",
      TRUE             ~ "Mid"
    ),
    time_bin = factor(time_bin, levels = c("Early", "Mid", "Late")),
    
    ## Labels with days for x-axis
    time_label = case_when(
      time_bin == "Early" ~ paste0("Early\n≤", round(early_cut), " d"),
      time_bin == "Mid"   ~ paste0("Mid\n>", round(early_cut), "–<", round(late_cut), " d"),
      time_bin == "Late"  ~ paste0("Late\n≥", round(late_cut), " d")
    ),
    time_label = factor(
      time_label,
      levels = c(
        paste0("Early\n≤", round(early_cut), " d"),
        paste0("Mid\n>", round(early_cut), "–<", round(late_cut), " d"),
        paste0("Late\n≥", round(late_cut), " d")
      )
    )
  )

## =========================
## 4. SUMMARY TABLE
## =========================

angle_bin_summary <- angles_binned %>%
  group_by(angle_type, time_bin, time_label) %>%
  summarise(
    n = n(),
    mean = mean(angle, na.rm = TRUE),
    sd = sd(angle, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

print(angle_bin_summary)

## Optional export
dir.create("outputs_model_triangles", showWarnings = FALSE)

write.csv(
  angle_bin_summary,
  "outputs_model_triangles/figure10_early_mid_late_summary.csv",
  row.names = FALSE
)

## =========================
## 5. COLORS
## =========================

angle_cols <- c(
  "Main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

## =========================
## 6. FIGURE 10 BARPLOT
## =========================

p_fig10 <- ggplot(
  angle_bin_summary,
  aes(x = time_label, y = mean, fill = angle_type)
) +
  geom_col(
    width = 0.7,
    color = NA
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.18,
    linewidth = 0.55,
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Mean angle (degrees)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.title = element_blank()
  )

print(p_fig10)

## =========================
## 7. SAVE FIGURE
## =========================

dir.create("plots/model3", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "plots/model3/figure10_early_mid_late_angles.png",
  plot = p_fig10,
  width = 8,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/figure10_early_mid_late_angles.pdf",
  plot = p_fig10,
  width = 8,
  height = 4,
  bg = "white"
)

## =========================================================
## FIGURE 10: EARLY / MID / LATE ANGLE BARPLOT
## Main version + genotype-panel supplementary version
## =========================================================

library(tidyverse)

## -------------------------
## 1. Colors
## -------------------------
angle_cols <- c(
  "Main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

## -------------------------
## 2. Define Early / Mid / Late stages
## Based on mean(day) ± 1 SD
## -------------------------
day_mean <- mean(angles_long$day, na.rm = TRUE)
day_sd   <- sd(angles_long$day, na.rm = TRUE)

early_cut <- day_mean - day_sd
late_cut  <- day_mean + day_sd

cat("\nStage definitions:\n")
cat("Early: day <=", round(early_cut, 1), "\n")
cat("Mid: day >", round(early_cut, 1), "and day <", round(late_cut, 1), "\n")
cat("Late: day >=", round(late_cut, 1), "\n")

angles_binned <- angles_long %>%
  mutate(
    time_bin = case_when(
      day <= early_cut ~ "Early",
      day >= late_cut  ~ "Late",
      TRUE             ~ "Mid"
    ),
    time_bin = factor(time_bin, levels = c("Early", "Mid", "Late")),
    angle_type = factor(angle_type, levels = c("Main", "NN1", "NN2")),
    genotype = factor(genotype)
  )

## -------------------------
## 3. Main summary table
## -------------------------
angle_bin_summary <- angles_binned %>%
  group_by(angle_type, time_bin) %>%
  summarise(
    n = n(),
    mean = mean(angle, na.rm = TRUE),
    sd = sd(angle, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

print(angle_bin_summary)

## -------------------------
## 4. Main Figure 10
## No days in x-axis labels
## -------------------------
p_fig10 <- ggplot(
  angle_bin_summary,
  aes(x = time_bin, y = mean, fill = angle_type)
) +
  geom_col(
    width = 0.7,
    color = NA
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.18,
    linewidth = 0.55,
    color = "black"
  ) +
  facet_wrap(~ angle_type, nrow = 1) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Mean angle (degrees)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.title = element_blank()
  )

print(p_fig10)

## -------------------------
## 5. Supplementary summary table by genotype
## -------------------------
angle_bin_genotype_summary <- angles_binned %>%
  group_by(genotype, angle_type, time_bin) %>%
  summarise(
    n = n(),
    mean = mean(angle, na.rm = TRUE),
    sd = sd(angle, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

print(angle_bin_genotype_summary)

## -------------------------
## 6. Supplementary genotype-panel figure
## Faceted by genotype
## -------------------------
p_fig10_genotype <- ggplot(
  angle_bin_genotype_summary,
  aes(x = time_bin, y = mean, fill = angle_type)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = NA
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.75),
    width = 0.18,
    linewidth = 0.45,
    color = "black"
  ) +
  facet_wrap(~ genotype, ncol = 4) +
  scale_fill_manual(values = angle_cols) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Mean angle (degrees)",
    fill = "Angle type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    plot.title = element_blank()
  )

print(p_fig10_genotype)

## -------------------------
## 7. Save outputs
## -------------------------
dir.create("plots/model3", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs_model_triangles", recursive = TRUE, showWarnings = FALSE)

write.csv(
  angle_bin_summary,
  "outputs_model_triangles/figure10_early_mid_late_summary.csv",
  row.names = FALSE
)

write.csv(
  angle_bin_genotype_summary,
  "outputs_model_triangles/figure10_early_mid_late_by_genotype_summary.csv",
  row.names = FALSE
)

ggsave(
  filename = "plots/model3/figure10_early_mid_late_angles.png",
  plot = p_fig10,
  width = 8,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/figure10_early_mid_late_angles.pdf",
  plot = p_fig10,
  width = 8,
  height = 4,
  bg = "white"
)

ggsave(
  filename = "plots/model3/figureS_early_mid_late_angles_by_genotype.png",
  plot = p_fig10_genotype,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/model3/figureS_early_mid_late_angles_by_genotype.pdf",
  plot = p_fig10_genotype,
  width = 10,
  height = 6,
  bg = "white"
)

############################################################
## FIGURE 10: Pooled early, mid, late angle plot
## Bars = observed mean angle
## Error bars = ±1 SE
##
## Bins:
## Early <= mean(day) - 1 SD
## Mid   > mean(day) - 1 SD and < mean(day) + 1 SD
## Late  >= mean(day) + 1 SD
##
## Edit:
## - "Main" changed to "main"
############################################################

## =========================
## 0. PACKAGES
## =========================

packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "scales"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(scales)

## =========================
## 1. PATHS
## =========================

setwd("/Users/amalia/Documents/GitHub/pattern_formation")

file_path <- "triangles_full_data.xlsx"

plot_dir <- "plots/model3"
output_dir <- "outputs_model_triangles"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 2. LOAD + CLEAN DATA
## =========================

tri <- read_excel(file_path) %>%
  clean_names() %>%
  mutate(
    genotype = factor(genotype),
    replica  = factor(replica),
    nubbin_id = interaction(genotype, replica, drop = TRUE, sep = "_"),
    day = as.numeric(day)
  ) %>%
  filter(
    angle_main > 0, angle_main < 180,
    angle_nn1  > 0, angle_nn1  < 180,
    angle_nn2  > 0, angle_nn2  < 180
  )

## =========================
## 3. LONG FORMAT
## =========================

angles_long <- tri %>%
  dplyr::select(
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
    )
  )

## =========================
## 4. DEFINE EARLY / MID / LATE BINS
## =========================

day_mean <- mean(angles_long$day, na.rm = TRUE)
day_sd   <- sd(angles_long$day, na.rm = TRUE)

early_cut <- day_mean - day_sd
late_cut  <- day_mean + day_sd

cat("\nStage definitions:\n")
cat("Early: day <=", round(early_cut, 1), "\n")
cat("Mid: day >", round(early_cut, 1), "and day <", round(late_cut, 1), "\n")
cat("Late: day >=", round(late_cut, 1), "\n")

angles_binned <- angles_long %>%
  mutate(
    time_bin = case_when(
      day <= early_cut ~ "Early",
      day >= late_cut  ~ "Late",
      TRUE             ~ "Mid"
    ),
    time_bin = factor(
      time_bin,
      levels = c("Early", "Mid", "Late")
    ),
    angle_type = factor(
      angle_type,
      levels = c("main", "NN1", "NN2")
    )
  )

## =========================
## 5. SUMMARY TABLE
## =========================

angle_bin_summary <- angles_binned %>%
  group_by(angle_type, time_bin) %>%
  summarise(
    n = n(),
    mean = mean(angle, na.rm = TRUE),
    sd = sd(angle, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

print(angle_bin_summary)

write.csv(
  angle_bin_summary,
  file.path(output_dir, "figure10_early_mid_late_summary_main_lowercase.csv"),
  row.names = FALSE
)

## =========================
## 6. COLORS
## =========================

angle_cols <- c(
  "main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

## =========================
## 7. FIGURE 10: POOLED BARPLOT
## =========================

p_fig10 <- ggplot(
  angle_bin_summary,
  aes(x = time_bin, y = mean, fill = angle_type)
) +
  geom_col(
    width = 0.7,
    color = NA
  ) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.18,
    linewidth = 0.55,
    color = "black"
  ) +
  facet_wrap(
    ~ angle_type,
    nrow = 1
  ) +
  scale_fill_manual(
    values = angle_cols,
    breaks = c("main", "NN1", "NN2")
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Mean angle (degrees)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.title = element_blank()
  )

print(p_fig10)

## =========================
## 8. SAVE FIGURE
## =========================

ggsave(
  filename = file.path(plot_dir, "figure10_early_mid_late_angles_main_lowercase.png"),
  plot = p_fig10,
  width = 8,
  height = 4,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "figure10_early_mid_late_angles_main_lowercase.pdf"),
  plot = p_fig10,
  width = 8,
  height = 4,
  bg = "white"
)

cat("\nDone. Saved pooled Figure 10 with 'main' lowercase.\n")

getwd()
)
