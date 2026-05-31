############################################################
## 2x2 figure: cross-process Gompertz parameter correlations
## + t0 lag boxplot
##
## Same stats as before:
## - panels (a-c): Spearman rank correlations
## - panel (d): lag against zero by one-sample t-test or
##   Wilcoxon depending on Shapiro
## - panel (d): among-colony lag by ANOVA+Tukey or
##   Kruskal-Wallis+Dunn depending on assumptions
##
## Visual fixes only:
## - lowercase panel labels: (a), (b), (c), (d)
## - t0 formatted with subscript
## - all axis titles bold, including plotmath labels
## - no delta symbol in asymptotic area label
## - remove stray dashed line outside panel d
############################################################

# ----------------------------
# 0. Working directory
# ----------------------------

setwd("/Users/amalia/Documents/GitHub/pattern_formation/gompertz")

# ----------------------------
# 1. Packages
# ----------------------------

packages <- c(
  "tidyverse",
  "readr",
  "ggplot2",
  "cowplot",
  "car",
  "FSA",
  "multcompView"
)

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

library(tidyverse)
library(readr)
library(ggplot2)
library(cowplot)
library(car)
library(FSA)
library(multcompView)

# ----------------------------
# 2. Output folders
# ----------------------------

dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 3. Load data
# ----------------------------

polyps_raw <- read_csv(
  "gompertz_individual_fit_parameters_polyps.csv",
  show_col_types = FALSE
)

area_raw <- read_csv(
  "gompertz_individual_fit_parameters_area.csv",
  show_col_types = FALSE
)

# ----------------------------
# 4. Prepare safe IDs and display labels
# ----------------------------

prepare_gompertz_data <- function(df) {
  
  geno_order <- sort(unique(df$Genotype))
  
  key <- tibble::tibble(
    Genotype      = geno_order,
    group_id      = sprintf("G%02d", seq_along(geno_order)),
    genotype_plot = paste0("SC", geno_order)
  )
  
  df %>%
    dplyr::left_join(key, by = "Genotype") %>%
    dplyr::mutate(
      group_id      = factor(group_id, levels = key$group_id),
      genotype_plot = factor(genotype_plot, levels = key$genotype_plot)
    )
}

polyps_raw <- prepare_gompertz_data(polyps_raw)
area_raw   <- prepare_gompertz_data(area_raw)

# ----------------------------
# 5. Check required columns
# ----------------------------

required_cols <- c("genotype_plot", "Replicate", "L", "k", "x0")

if (!all(required_cols %in% colnames(polyps_raw))) {
  stop(
    "Polyps file is missing required columns: ",
    paste(setdiff(required_cols, colnames(polyps_raw)), collapse = ", ")
  )
}

if (!all(required_cols %in% colnames(area_raw))) {
  stop(
    "Area file is missing required columns: ",
    paste(setdiff(required_cols, colnames(area_raw)), collapse = ", ")
  )
}

# ----------------------------
# 6. Rename x0 -> t0 in working objects
# ----------------------------

polyps <- polyps_raw %>%
  dplyr::rename(t0 = x0)

area <- area_raw %>%
  dplyr::rename(t0 = x0)

# ----------------------------
# 7. Merge area + polyp parameters by colony and replicate
# ----------------------------

params_poly <- polyps %>%
  dplyr::select(
    genotype_plot,
    Replicate,
    L_poly  = L,
    k_poly  = k,
    t0_poly = t0
  )

params_area <- area %>%
  dplyr::select(
    genotype_plot,
    Replicate,
    L_area  = L,
    k_area  = k,
    t0_area = t0
  )

params_combined <- dplyr::left_join(
  params_area,
  params_poly,
  by = c("genotype_plot", "Replicate")
) %>%
  dplyr::mutate(
    t0_lag = t0_poly - t0_area
  )

cat("\n--- Combined rows ---\n")
print(nrow(params_combined))

cat("\n--- Number of incomplete rows after join ---\n")
print(sum(!complete.cases(params_combined)))

params_combined <- params_combined %>%
  dplyr::filter(complete.cases(.))

cat("\n--- Rows retained after removing incomplete rows ---\n")
print(nrow(params_combined))

write_csv(
  params_combined,
  "outputs/gompertz_area_polyp_parameters_combined.csv"
)

# ----------------------------
# 8. Reusable genotype colors
# ----------------------------

okabe_ito_base <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#000000"  # black
)

fixed_cols <- c(
  "SC1" = "#0072B2",
  "SC2" = "#E69F00",
  "SC5" = "#009E73"
)

geno_levels <- levels(polyps$genotype_plot)

remaining_levels <- setdiff(geno_levels, names(fixed_cols))
remaining_pool   <- setdiff(okabe_ito_base, unname(fixed_cols))

if (length(remaining_levels) > length(remaining_pool)) {
  remaining_pool <- c(
    remaining_pool,
    "#56B4E9",
    "#D55E00",
    "#CC79A7",
    "#999999"
  )
}

other_cols  <- setNames(remaining_pool[seq_along(remaining_levels)], remaining_levels)
geno_colors <- c(fixed_cols, other_cols)
geno_colors <- geno_colors[geno_levels]

# ----------------------------
# 9. Bold plotmath labels
# ----------------------------

lab_x_L  <- expression(bold("Asymptotic area (mm"^2*")"))
lab_y_L  <- expression(bold("Asymptotic polyp number ("*L*")"))

lab_x_k  <- expression(bold("Area growth rate ("*k*")"))
lab_y_k  <- expression(bold("Polyp growth rate ("*k*")"))

lab_x_t0 <- expression(bold("Area inflection day ("*t[0]*")"))
lab_y_t0 <- expression(bold("Polyp inflection day ("*t[0]*")"))

lab_x_lag <- "Colony"
lab_y_lag <- expression(bold(t[0]~"lag (polyp - area, days)"))

# ----------------------------
# 10. Themes
# ----------------------------

base_size <- 14

common_scatter_theme <- theme_minimal(base_size = base_size) +
  theme(
    axis.title.x     = element_text(size = base_size + 2, face = "bold", colour = "black"),
    axis.title.y     = element_text(size = base_size + 2, face = "bold", colour = "black"),
    axis.text        = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    legend.title     = element_text(size = base_size, face = "bold"),
    legend.text      = element_text(size = base_size - 1),
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(8, 8, 8, 8)
  )

lag_theme <- theme_minimal(base_size = base_size) +
  theme(
    axis.title.x       = element_text(size = base_size + 2, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 2, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(8, 8, 8, 8)
  )

# ----------------------------
# 11. Cross-process Spearman correlations
# ----------------------------

cor_L <- cor.test(
  params_combined$L_area,
  params_combined$L_poly,
  method = "spearman",
  exact = FALSE
)

cor_k <- cor.test(
  params_combined$k_area,
  params_combined$k_poly,
  method = "spearman",
  exact = FALSE
)

cor_t0 <- cor.test(
  params_combined$t0_area,
  params_combined$t0_poly,
  method = "spearman",
  exact = FALSE
)

make_cor_label <- function(ct) {
  paste0(
    "Spearman \u03C1 = ", round(unname(ct$estimate), 2),
    "\n",
    "p = ", format.pval(ct$p.value, digits = 2, eps = 1e-4)
  )
}

cor_labels <- list(
  L  = make_cor_label(cor_L),
  k  = make_cor_label(cor_k),
  t0 = make_cor_label(cor_t0)
)

cor_summary <- tibble::tibble(
  comparison = c("L_area vs L_poly", "k_area vs k_poly", "t0_area vs t0_poly"),
  method = "Spearman rank correlation",
  rho = c(
    unname(cor_L$estimate),
    unname(cor_k$estimate),
    unname(cor_t0$estimate)
  ),
  statistic = c(
    unname(cor_L$statistic),
    unname(cor_k$statistic),
    unname(cor_t0$statistic)
  ),
  p_value = c(
    cor_L$p.value,
    cor_k$p.value,
    cor_t0$p.value
  )
)

cat("\n==================================================\n")
cat("SPEARMAN CORRELATIONS\n")
cat("==================================================\n")
print(cor_summary)

write_csv(
  cor_summary,
  "outputs/gompertz_area_polyp_spearman_correlations.csv"
)

# ----------------------------
# 12. Scatter plots
# ----------------------------

p_L <- ggplot(
  params_combined,
  aes(x = L_area, y = L_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = cor_labels$L,
    hjust = -0.1,
    vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = lab_x_L,
    y = lab_y_L,
    color = "Colony"
  ) +
  common_scatter_theme

p_k <- ggplot(
  params_combined,
  aes(x = k_area, y = k_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = cor_labels$k,
    hjust = -0.1,
    vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = lab_x_k,
    y = lab_y_k,
    color = "Colony"
  ) +
  common_scatter_theme

p_t0 <- ggplot(
  params_combined,
  aes(x = t0_area, y = t0_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.8,
    color = "black"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = cor_labels$t0,
    hjust = -0.1,
    vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = lab_x_t0,
    y = lab_y_t0,
    color = "Colony"
  ) +
  common_scatter_theme

# ----------------------------
# 13. Lag statistics
# ----------------------------

# 13a. Test whether lag differs from zero
lag_shapiro <- shapiro.test(params_combined$t0_lag)

if (lag_shapiro$p.value > 0.05) {
  
  lag_zero_test <- t.test(params_combined$t0_lag, mu = 0)
  
  lag_zero_method <- "One-sample t-test on paired differences"
  lag_zero_stat_type <- "t"
  lag_zero_stat <- unname(lag_zero_test$statistic)
  lag_zero_df <- unname(lag_zero_test$parameter)
  lag_zero_p <- lag_zero_test$p.value
  
} else {
  
  lag_zero_test <- wilcox.test(
    params_combined$t0_lag,
    mu = 0,
    exact = FALSE
  )
  
  lag_zero_method <- "One-sample Wilcoxon signed-rank test on paired differences"
  lag_zero_stat_type <- "V"
  lag_zero_stat <- unname(lag_zero_test$statistic)
  lag_zero_df <- NA_real_
  lag_zero_p <- lag_zero_test$p.value
}

lag_zero_summary <- tibble::tibble(
  test_question = "Does t0_lag = t0_poly - t0_area differ from zero?",
  method = lag_zero_method,
  statistic_type = lag_zero_stat_type,
  statistic = lag_zero_stat,
  df = lag_zero_df,
  p_value = lag_zero_p,
  shapiro_p_lag = lag_shapiro$p.value,
  normality_assumption_met = lag_shapiro$p.value > 0.05
)

# 13b. Test colony effect on lag
lag_lm <- lm(t0_lag ~ genotype_plot, data = params_combined)

lag_resid_shapiro <- shapiro.test(residuals(lag_lm))
lag_levene <- car::leveneTest(t0_lag ~ genotype_plot, data = params_combined)

lag_use_anova <- lag_resid_shapiro$p.value > 0.05 &&
  lag_levene[["Pr(>F)"]][1] > 0.05

if (lag_use_anova) {
  
  lag_model <- aov(t0_lag ~ genotype_plot, data = params_combined)
  lag_posthoc <- TukeyHSD(lag_model)
  lag_letters <- multcompLetters4(lag_model, lag_posthoc)
  
  lag_letters_df <- data.frame(
    genotype_plot = names(lag_letters$genotype_plot$Letters),
    Letter = as.character(lag_letters$genotype_plot$Letters),
    stringsAsFactors = FALSE
  )
  
  lag_aov_tab <- summary(lag_model)[[1]]
  
  lag_group_method <- "One-way ANOVA + Tukey HSD"
  lag_group_stat_type <- "F"
  lag_group_stat <- unname(as.numeric(lag_aov_tab[1, "F value"]))
  lag_group_df1 <- unname(as.numeric(lag_aov_tab[1, "Df"]))
  lag_group_df2 <- unname(as.numeric(lag_aov_tab[2, "Df"]))
  lag_group_p <- unname(as.numeric(lag_aov_tab[1, "Pr(>F)"]))
  
  lag_pairwise <- as.data.frame(lag_posthoc$genotype_plot)
  lag_pairwise$comparison <- rownames(lag_pairwise)
  rownames(lag_pairwise) <- NULL
  
  lag_pairwise <- lag_pairwise %>%
    as_tibble() %>%
    dplyr::rename(
      diff      = diff,
      conf_low  = lwr,
      conf_high = upr,
      p_adj     = `p adj`
    ) %>%
    dplyr::mutate(
      method = "Tukey HSD",
      .before = 1
    )
  
} else {
  
  lag_model <- kruskal.test(t0_lag ~ genotype_plot, data = params_combined)
  
  lag_posthoc <- FSA::dunnTest(
    t0_lag ~ genotype_plot,
    data = params_combined,
    method = "bh"
  )
  
  pvals_lag <- lag_posthoc$res$P.adj
  comp_names <- lag_posthoc$res$Comparison
  comp_names <- gsub(" - ", "-", comp_names)
  comp_names <- gsub(" ", "", comp_names)
  names(pvals_lag) <- comp_names
  
  lag_letters <- multcompView::multcompLetters(pvals_lag)
  
  lag_letters_df <- data.frame(
    genotype_plot = names(lag_letters$Letters),
    Letter = as.character(lag_letters$Letters),
    stringsAsFactors = FALSE
  )
  
  lag_group_method <- "Kruskal-Wallis + Dunn post hoc test with Benjamini-Hochberg correction"
  lag_group_stat_type <- "chi_squared"
  lag_group_stat <- unname(as.numeric(lag_model$statistic))
  lag_group_df1 <- unname(as.numeric(lag_model$parameter))
  lag_group_df2 <- NA_real_
  lag_group_p <- lag_model$p.value
  
  lag_pairwise <- lag_posthoc$res %>%
    as_tibble() %>%
    dplyr::rename(
      comparison = Comparison,
      statistic  = Z,
      p_unadj    = P.unadj,
      p_adj      = P.adj
    ) %>%
    dplyr::mutate(
      method = "Dunn BH-adjusted",
      .before = 1
    )
}

lag_group_summary <- tibble::tibble(
  test_question = "Does t0_lag differ among colonies?",
  method = lag_group_method,
  statistic_type = lag_group_stat_type,
  statistic = lag_group_stat,
  df1 = lag_group_df1,
  df2 = lag_group_df2,
  p_value = lag_group_p,
  shapiro_p_residuals = lag_resid_shapiro$p.value,
  levene_p = lag_levene[["Pr(>F)"]][1],
  assumptions_met_for_anova = lag_use_anova
)

cat("\n==================================================\n")
cat("LAG AGAINST ZERO\n")
cat("==================================================\n")
print(lag_zero_summary)

cat("\n==================================================\n")
cat("AMONG-COLONY LAG COMPARISON\n")
cat("==================================================\n")
print(lag_group_summary)

cat("\n==================================================\n")
cat("POST HOC LETTERS FOR PANEL (d)\n")
cat("==================================================\n")
print(lag_letters_df)

cat("\n==================================================\n")
cat("PAIRWISE POST HOC COMPARISONS FOR PANEL (d)\n")
cat("==================================================\n")
print(lag_pairwise)

write_csv(
  lag_zero_summary,
  "outputs/t0_lag_zero_test_summary.csv"
)

write_csv(
  lag_group_summary,
  "outputs/t0_lag_among_colony_test_summary.csv"
)

write_csv(
  lag_letters_df,
  "outputs/t0_lag_posthoc_letters.csv"
)

write_csv(
  lag_pairwise,
  "outputs/t0_lag_pairwise_posthoc.csv"
)

# ----------------------------
# 14. Label positions for lag letters
# ----------------------------

lag_span <- diff(range(params_combined$t0_lag, na.rm = TRUE))
if (lag_span == 0) lag_span <- 1

label_lag <- params_combined %>%
  dplyr::group_by(genotype_plot) %>%
  dplyr::summarise(
    y_max = max(t0_lag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(lag_letters_df, by = "genotype_plot") %>%
  dplyr::mutate(
    label_y = max(params_combined$t0_lag, na.rm = TRUE) + 0.08 * lag_span
  )

# ----------------------------
# 15. Lag boxplot
# ----------------------------

p_lag <- ggplot(
  params_combined,
  aes(x = genotype_plot, y = t0_lag)
) +
  geom_boxplot(
    fill = "cadetblue3",
    color = "black",
    width = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.7,
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.5,
    fill = "red",
    color = "black",
    stroke = 0.5
  ) +
  geom_text(
    data = label_lag,
    aes(x = genotype_plot, y = label_y, label = Letter),
    inherit.aes = FALSE,
    size = 4.2,
    vjust = 0
  ) +
  labs(
    x = lab_x_lag,
    y = lab_y_lag
  ) +
  coord_cartesian(
    ylim = c(
      min(params_combined$t0_lag, na.rm = TRUE) - 0.05 * lag_span,
      max(label_lag$label_y, na.rm = TRUE) + 0.08 * lag_span
    ),
    clip = "on"
  ) +
  lag_theme

# ----------------------------
# 16. Remove legends from panel plots
# ----------------------------

p_L_noleg  <- p_L  + theme(legend.position = "none")
p_k_noleg  <- p_k  + theme(legend.position = "none")
p_t0_noleg <- p_t0 + theme(legend.position = "none")

# ----------------------------
# 17. Extract shared legend
# ----------------------------

legend_plot <- p_L +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold")
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 3)
    )
  )

shared_legend <- cowplot::get_legend(legend_plot)

# ----------------------------
# 18. Assemble 2x2 figure with lowercase labels
# ----------------------------

figure_body <- cowplot::plot_grid(
  cowplot::plot_grid(
    p_L_noleg,
    p_k_noleg,
    labels = c("(a)", "(b)"),
    label_size = base_size + 3,
    label_fontface = "bold",
    label_x = 0.01,
    label_y = 0.99,
    hjust = 0,
    vjust = 1,
    ncol = 2,
    align = "hv"
  ),
  cowplot::plot_grid(
    p_t0_noleg,
    p_lag,
    labels = c("(c)", "(d)"),
    label_size = base_size + 3,
    label_fontface = "bold",
    label_x = 0.01,
    label_y = 0.99,
    hjust = 0,
    vjust = 1,
    ncol = 2,
    align = "hv"
  ),
  ncol = 1,
  align = "v"
)

figure_main_2x2 <- cowplot::plot_grid(
  figure_body,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.13)
)

print(figure_main_2x2)

# ----------------------------
# 19. Save figure
# ----------------------------

ggsave(
  filename = "plots/gompertz_area_polyp_parameter_correlations_t0_lag_FINAL.png",
  plot = figure_main_2x2,
  width = 9.5,
  height = 9.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "plots/gompertz_area_polyp_parameter_correlations_t0_lag_FINAL.pdf",
  plot = figure_main_2x2,
  width = 9.5,
  height = 9.5,
  bg = "white"
)

# ----------------------------
# 20. Compact console summary for manuscript writing
# ----------------------------

cat("\n\n==================================================\n")
cat("MANUSCRIPT STATS SUMMARY\n")
cat("==================================================\n")

cat("\nPanels (a-c): Spearman rank correlations\n")
print(cor_summary)

cat("\nPanel (d): Lag against zero\n")
print(lag_zero_summary)

cat("\nPanel (d): Among-colony lag comparison\n")
print(lag_group_summary)

cat("\nPanel (d): Letters\n")
print(lag_letters_df)

cat("\nSaved outputs:\n")
cat("- outputs/gompertz_area_polyp_spearman_correlations.csv\n")
cat("- outputs/t0_lag_zero_test_summary.csv\n")
cat("- outputs/t0_lag_among_colony_test_summary.csv\n")
cat("- outputs/t0_lag_posthoc_letters.csv\n")
cat("- outputs/t0_lag_pairwise_posthoc.csv\n")
cat("- plots/gompertz_area_polyp_parameter_correlations_t0_lag_FINAL.png/pdf\n")

############################################################
## END
############################################################