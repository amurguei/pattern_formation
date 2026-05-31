############################################################
## NN2 / NN1 DISTANCE RATIO PANEL
############################################################

## =========================
## 0. PACKAGES
## =========================
packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "patchwork",
  "scales"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(patchwork)
library(scales)

## =========================
## 1. WORKING DIRECTORY
## =========================
setwd("/Users/amalia/Documents/GitHub/pattern_formation")

## Create output folder
if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)

## =========================
## 2. LOAD + CLEAN DATA
## =========================
file_path <- "triangles_full_data.xlsx"

tri <- read_excel(file_path) %>%
  clean_names() %>%
  mutate(
    genotype  = factor(genotype),
    replica   = factor(replica),
    nubbin_id = interaction(genotype, replica, drop = TRUE, sep = "_"),
    day       = as.numeric(day)
  ) %>%
  filter(
    !is.na(dist_to_nn1),
    !is.na(dist_to_nn2),
    dist_to_nn1 > 0,
    dist_to_nn2 > 0
  ) %>%
  mutate(
    dist_ratio = dist_to_nn2 / dist_to_nn1
  )

## Optional sanity check
summary(tri$dist_ratio)

## =========================
## 3. DATA FOR RATIO PANEL
## =========================
ratio_data <- tri %>%
  filter(
    is.finite(dist_ratio),
    !is.na(dist_ratio),
    dist_ratio > 0
  ) %>%
  mutate(
    genotype = factor(genotype)
  )

## Visual limits
ratio_ylim <- c(0.95, 2.3)

## =========================
## 4. STYLE
## =========================
col_overall <- "#56B4E9"     # light blue
col_colony  <- "cadetblue3"  # by-colony violins
col_mean    <- "red"         # red mean dots
col_smooth  <- "#2C5EFF"     # smooth line
col_points  <- "black"

theme_ratio <- theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold", colour = "black"),
    axis.text  = element_text(size = 13, colour = "black"),
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(8, 10, 8, 10)
  )

## =========================
## 5. PANEL (a): OVERALL DISTRIBUTION
## =========================
p_ratio_A <- ggplot(ratio_data, aes(x = "All data", y = dist_ratio)) +
  geom_violin(
    fill = col_overall,
    alpha = 0.75,
    color = NA,
    trim = TRUE,
    width = 0.85
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = col_mean,
    color = "black",
    stroke = 0.4
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = NULL,
    y = expression(NN[2] / NN[1]~"distance ratio"),
    title = "(a) Overall distribution"
  ) +
  theme_ratio

## =========================
## 6. PANEL (b): BY COLONY
## =========================
p_ratio_B <- ggplot(ratio_data, aes(x = genotype, y = dist_ratio)) +
  geom_violin(
    fill = col_colony,
    alpha = 0.75,
    color = NA,
    trim = TRUE,
    width = 0.85
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.8,
    fill = col_mean,
    color = "black",
    stroke = 0.4
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = "Colony",
    y = expression(NN[2] / NN[1]~"distance ratio"),
    title = "(b) By colony"
  ) +
  theme_ratio

## =========================
## 7. PANEL (c): THROUGH TIME
## =========================
p_ratio_C <- ggplot(ratio_data, aes(x = day, y = dist_ratio)) +
  geom_point(
    alpha = 0.12,
    size = 0.8,
    color = col_points
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 1.1,
    color = col_smooth,
    fill = "grey75"
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
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
    y = expression(NN[2] / NN[1]~"distance ratio"),
    title = "(c) Through time"
  ) +
  theme_ratio +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85")
  )

## =========================
## 8. COMBINED LAYOUT
## =========================
p_ratio_panel <- (p_ratio_A | p_ratio_B) / p_ratio_C +
  plot_layout(
    heights = c(1, 1.05),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(6, 8, 6, 8)
  )

print(p_ratio_panel)

## =========================
## 9. SAVE FIGURE
## =========================
ggsave(
  filename = "plots/nn_ratio_panel_2plus1_lowercase.png",
  plot = p_ratio_panel,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_2plus1_lowercase.pdf",
  plot = p_ratio_panel,
  width = 11,
  height = 8.5,
  bg = "white"
)

## Optional: save individual panels too
ggsave(
  filename = "plots/nn_ratio_panel_a_overall.png",
  plot = p_ratio_A,
  width = 5,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_b_colony.png",
  plot = p_ratio_B,
  width = 6,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_c_time.png",
  plot = p_ratio_C,
  width = 10,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

############################################################
## END
############################################################


############################################################
## NN2 / NN1 DISTANCE RATIO PANEL
## Full self-contained script
## Saves SE and SD versions
############################################################

## =========================
## 0. PACKAGES
## =========================
packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "patchwork",
  "scales"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readxl)
library(janitor)
library(patchwork)
library(scales)

## =========================
## 1. WORKING DIRECTORY
## =========================
setwd("/Users/amalia/Documents/GitHub/pattern_formation")

if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

## =========================
## 2. LOAD + CLEAN DATA
## =========================
file_path <- "triangles_full_data.xlsx"

tri <- read_excel(file_path) %>%
  clean_names() %>%
  mutate(
    genotype  = factor(genotype),
    replica   = factor(replica),
    nubbin_id = interaction(genotype, replica, drop = TRUE, sep = "_"),
    day       = as.numeric(day)
  ) %>%
  filter(
    !is.na(dist_to_nn1),
    !is.na(dist_to_nn2),
    dist_to_nn1 > 0,
    dist_to_nn2 > 0
  ) %>%
  mutate(
    dist_ratio = dist_to_nn2 / dist_to_nn1
  )

ratio_data <- tri %>%
  filter(
    is.finite(dist_ratio),
    !is.na(dist_ratio),
    dist_ratio > 0
  )

cat("\nSummary of NN2/NN1 ratio:\n")
print(summary(ratio_data$dist_ratio))

## =========================
## 3. SUMMARY TABLES
## =========================

ratio_summary_overall <- ratio_data %>%
  summarise(
    n = n(),
    mean_ratio = mean(dist_ratio, na.rm = TRUE),
    sd_ratio = sd(dist_ratio, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n)
  )

ratio_summary_colony <- ratio_data %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    mean_ratio = mean(dist_ratio, na.rm = TRUE),
    sd_ratio = sd(dist_ratio, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    .groups = "drop"
  )

ratio_summary_day <- ratio_data %>%
  group_by(day) %>%
  summarise(
    n = n(),
    mean_ratio = mean(dist_ratio, na.rm = TRUE),
    sd_ratio = sd(dist_ratio, na.rm = TRUE),
    se_ratio = sd_ratio / sqrt(n),
    .groups = "drop"
  )

write_csv(ratio_summary_overall, "outputs/nn2_nn1_ratio_summary_overall.csv")
write_csv(ratio_summary_colony, "outputs/nn2_nn1_ratio_summary_by_colony.csv")
write_csv(ratio_summary_day, "outputs/nn2_nn1_ratio_summary_by_day.csv")

print(ratio_summary_overall)
print(ratio_summary_colony)
print(ratio_summary_day)

## =========================
## 4. STYLE
## =========================

ratio_ylim <- c(0.95, 2.3)

col_overall <- "#56B4E9"
col_colony  <- "cadetblue3"
col_mean    <- "red"
col_time    <- "#56B4E9"
col_points  <- "black"

theme_ratio <- theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold", colour = "black"),
    axis.text  = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 17, face = "bold", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(8, 10, 8, 10)
  )

## =========================
## 5. PANEL (a): OVERALL DISTRIBUTION
## =========================

p_ratio_A <- ggplot(ratio_data, aes(x = "All data", y = dist_ratio)) +
  geom_violin(
    fill = col_overall,
    alpha = 0.75,
    color = NA,
    trim = TRUE,
    width = 0.85
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3,
    fill = col_mean,
    color = "black",
    stroke = 0.4
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = NULL,
    y = expression(NN[2] / NN[1]~"distance ratio"),
    title = "(a) Overall distribution"
  ) +
  theme_ratio

## =========================
## 6. PANEL (b): BY COLONY
## =========================

p_ratio_B <- ggplot(ratio_data, aes(x = genotype, y = dist_ratio)) +
  geom_violin(
    fill = col_colony,
    alpha = 0.75,
    color = NA,
    trim = TRUE,
    width = 0.85
  ) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    color = "black",
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 2.8,
    fill = col_mean,
    color = "black",
    stroke = 0.4
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  coord_cartesian(ylim = ratio_ylim) +
  labs(
    x = "Colony",
    y = expression(NN[2] / NN[1]~"distance ratio"),
    title = "(b) By colony"
  ) +
  theme_ratio

## =========================
## 7. PANEL (c): THROUGH TIME FUNCTION
## =========================

make_ratio_time_panel <- function(error_type = c("SE", "SD")) {
  
  error_type <- match.arg(error_type)
  
  if (error_type == "SE") {
    time_df <- ratio_summary_day %>%
      mutate(
        lower = mean_ratio - se_ratio,
        upper = mean_ratio + se_ratio,
        subtitle_text = "Mean ± 1 SE"
      )
  } else {
    time_df <- ratio_summary_day %>%
      mutate(
        lower = mean_ratio - sd_ratio,
        upper = mean_ratio + sd_ratio,
        subtitle_text = "Mean ± 1 SD"
      )
  }
  
  ggplot(time_df, aes(x = day, y = mean_ratio)) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.6,
      color = "black"
    ) +
    geom_line(
      linewidth = 1.2,
      color = col_time
    ) +
    geom_point(
      size = 3,
      shape = 21,
      fill = col_time,
      color = "black",
      stroke = 0.4
    ) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 1.5,
      linewidth = 0.7,
      color = col_time
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
      y = expression(NN[2] / NN[1]~"distance ratio"),
      title = "(c) Through time"
    ) +
    theme_ratio +
    theme(
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85")
    )
}

p_ratio_C_se <- make_ratio_time_panel("SE")
p_ratio_C_sd <- make_ratio_time_panel("SD")

## =========================
## 8. COMBINED LAYOUTS
## =========================

p_ratio_panel_se <- (p_ratio_A | p_ratio_B) / p_ratio_C_se +
  plot_layout(
    heights = c(1, 1.05),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(6, 8, 6, 8)
  )

p_ratio_panel_sd <- (p_ratio_A | p_ratio_B) / p_ratio_C_sd +
  plot_layout(
    heights = c(1, 1.05),
    widths = c(1, 1)
  ) &
  theme(
    plot.margin = margin(6, 8, 6, 8)
  )

print(p_ratio_panel_se)
print(p_ratio_panel_sd)

## =========================
## 9. SAVE FIGURES
## =========================

ggsave(
  filename = "plots/nn_ratio_panel_2plus1_SE.png",
  plot = p_ratio_panel_se,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_2plus1_SE.pdf",
  plot = p_ratio_panel_se,
  width = 11,
  height = 8.5,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_2plus1_SD.png",
  plot = p_ratio_panel_sd,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_2plus1_SD.pdf",
  plot = p_ratio_panel_sd,
  width = 11,
  height = 8.5,
  bg = "white"
)

## Individual panel C versions
ggsave(
  filename = "plots/nn_ratio_panel_c_time_SE.png",
  plot = p_ratio_C_se,
  width = 10,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/nn_ratio_panel_c_time_SD.png",
  plot = p_ratio_C_sd,
  width = 10,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

cat("\nDone. NN2/NN1 ratio panels with SE and SD saved.\n")

