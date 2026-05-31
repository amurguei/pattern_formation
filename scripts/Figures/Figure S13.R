############################################################
## TIME TRENDS IN ANGLES
## main / NN1 / NN2
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
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
    ),
    angle_type = factor(angle_type, levels = c("main", "NN1", "NN2"))
  )

## =========================
## 4. COLORS
## =========================

angle_cols <- c(
  "main" = "#000000",
  "NN1"  = "#0072B2",
  "NN2"  = "#D55E00"
)

## =========================
## 5. TIME TREND PLOT
## =========================

p_time <- ggplot(
  angles_long,
  aes(x = day, y = angle, group = nubbin_id, color = angle_type)
) +
  geom_line(
    alpha = 0.12,
    linewidth = 0.4
  ) +
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
  facet_wrap(
    ~ angle_type,
    ncol = 1
  ) +
  scale_color_manual(
    values = angle_cols,
    breaks = c("main", "NN1", "NN2")
  ) +
  scale_y_continuous(
    breaks = seq(0, 180, by = 20),
    limits = c(0, 180),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8)
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
  ) +
  labs(
    x = "Day",
    y = "Angle (degrees)",
    title = "Time trends in angles"
  )

print(p_time)

## =========================
## 6. SAVE FIGURE
## =========================

ggsave(
  filename = file.path(plot_dir, "figS_13_angle_time_trends_main_lowercase.png"),
  plot = p_time,
  width = 7,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(plot_dir, "figS_13_angle_time_trends_main_lowercase.pdf"),
  plot = p_time,
  width = 7,
  height = 9,
  bg = "white"
)

cat("\nDone. Saved time trend graph with 'main' lowercase.\n")
