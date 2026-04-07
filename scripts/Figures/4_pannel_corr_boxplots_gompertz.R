## ============================================================
## SELF-CONTAINED SCRIPT
## 2x2 figure: cross-process correlations + lag boxplot
## Uses t0 instead of x0
## ============================================================

# ----------------------------
# 0. Working directory
# ----------------------------
setwd("/Users/amalia/Documents/GitHub/pattern_formation/gompertz")

# ----------------------------
# 1. Packages
# ----------------------------
library(tidyverse)
library(readr)
library(ggplot2)
library(cowplot)
library(car)
library(FSA)
library(multcompView)

# ----------------------------
# 2. Load data
# ----------------------------
polyps_raw <- read_csv("gompertz_individual_fit_parameters_polyps.csv")
area_raw   <- read_csv("gompertz_individual_fit_parameters_area.csv")

# ----------------------------
# 3. Prepare safe IDs and display labels
# ----------------------------
prepare_gompertz_data <- function(df) {
  geno_order <- sort(unique(df$Genotype))
  
  key <- tibble(
    Genotype      = geno_order,
    group_id      = sprintf("G%02d", seq_along(geno_order)),
    genotype_plot = paste0("SC", geno_order)
  )
  
  df %>%
    left_join(key, by = "Genotype") %>%
    mutate(
      group_id      = factor(group_id, levels = key$group_id),
      genotype_plot = factor(genotype_plot, levels = key$genotype_plot)
    )
}

polyps_raw <- prepare_gompertz_data(polyps_raw)
area_raw   <- prepare_gompertz_data(area_raw)

# ----------------------------
# 4. Check required columns
# ----------------------------
required_cols <- c("genotype_plot", "Replicate", "L", "k", "x0")

if (!all(required_cols %in% colnames(polyps_raw))) {
  stop("Polyps file is missing required columns: ",
       paste(setdiff(required_cols, colnames(polyps_raw)), collapse = ", "))
}

if (!all(required_cols %in% colnames(area_raw))) {
  stop("Area file is missing required columns: ",
       paste(setdiff(required_cols, colnames(area_raw)), collapse = ", "))
}

# ----------------------------
# 5. Rename x0 -> t0 in working objects
# ----------------------------
polyps <- polyps_raw %>%
  rename(t0 = x0)

area <- area_raw %>%
  rename(t0 = x0)

# ----------------------------
# 6. Merge area + polyp parameters by colony and replicate
# ----------------------------
params_poly <- polyps %>%
  select(genotype_plot, Replicate,
         L_poly = L, k_poly = k, t0_poly = t0)

params_area <- area %>%
  select(genotype_plot, Replicate,
         L_area = L, k_area = k, t0_area = t0)

params_combined <- left_join(
  params_area,
  params_poly,
  by = c("genotype_plot", "Replicate")
) %>%
  mutate(
    t0_lag = t0_poly - t0_area
  )

cat("\n--- Combined rows ---\n")
print(nrow(params_combined))

cat("\n--- Number of incomplete rows after join ---\n")
print(sum(!complete.cases(params_combined)))

# ----------------------------
# 7. Reusable genotype colors (Okabe–Ito style)
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
# 8. Themes
# ----------------------------
base_size <- 14

common_scatter_theme <- theme_minimal(base_size = base_size) +
  theme(
    plot.title         = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    axis.title.x       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(linewidth = 0.25),
    legend.title       = element_text(size = base_size, face = "bold"),
    legend.text        = element_text(size = base_size - 1),
    plot.background    = element_rect(fill = "white", color = NA)
  )

lag_theme <- theme_minimal(base_size = base_size) +
  theme(
    axis.title.x       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.background    = element_rect(fill = "white", color = NA)
  )

# ----------------------------
# 9. Cross-process Spearman correlations
# ----------------------------
cor_L  <- cor.test(params_combined$L_area,  params_combined$L_poly,
                   method = "spearman", exact = FALSE)

cor_k  <- cor.test(params_combined$k_area,  params_combined$k_poly,
                   method = "spearman", exact = FALSE)

cor_t0 <- cor.test(params_combined$t0_area, params_combined$t0_poly,
                   method = "spearman", exact = FALSE)

make_cor_label <- function(ct) {
  paste0(
    "Spearman \u03C1 = ", round(unname(ct$estimate), 2),
    "\np = ", format.pval(ct$p.value, digits = 2, eps = 1e-4)
  )
}

cor_labels <- list(
  L  = make_cor_label(cor_L),
  k  = make_cor_label(cor_k),
  t0 = make_cor_label(cor_t0)
)

# ----------------------------
# 10. Cross-process scatter plots
# ----------------------------
p_L <- ggplot(
  params_combined,
  aes(x = L_area, y = L_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels$L,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Asymptotic \u0394 area (mm\u00B2)",
    y = "Asymptotic polyp number (L)",
    color = "Colony"
  ) +
  common_scatter_theme

p_k <- ggplot(
  params_combined,
  aes(x = k_area, y = k_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels$k,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Area growth rate (k)",
    y = "Polyp growth rate (k)",
    color = "Colony"
  ) +
  common_scatter_theme

p_t0 <- ggplot(
  params_combined,
  aes(x = t0_area, y = t0_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels$t0,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Area inflection day (t0)",
    y = "Polyp inflection day (t0)",
    color = "Colony"
  ) +
  common_scatter_theme

# ----------------------------
# 11. Lag statistics
# ----------------------------

# 11a. Test whether lag differs from 0
lag_shapiro <- shapiro.test(params_combined$t0_lag)

if (lag_shapiro$p.value > 0.05) {
  lag_zero_test <- t.test(params_combined$t0_lag, mu = 0)
  lag_zero_method <- "One-sample t-test"
  lag_zero_stat_type <- "t"
  lag_zero_stat <- unname(lag_zero_test$statistic)
  lag_zero_df <- unname(lag_zero_test$parameter)
  lag_zero_p <- lag_zero_test$p.value
} else {
  lag_zero_test <- wilcox.test(params_combined$t0_lag, mu = 0, exact = FALSE)
  lag_zero_method <- "One-sample Wilcoxon signed-rank test"
  lag_zero_stat_type <- "V"
  lag_zero_stat <- unname(lag_zero_test$statistic)
  lag_zero_df <- NA_real_
  lag_zero_p <- lag_zero_test$p.value
}

# 11b. Test colony effect on lag
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
  
  lag_group_method <- "ANOVA + Tukey"
  lag_group_stat_type <- "F"
  lag_group_stat <- unname(as.numeric(lag_aov_tab[1, "F value"]))
  lag_group_df1 <- unname(as.numeric(lag_aov_tab[1, "Df"]))
  lag_group_df2 <- unname(as.numeric(lag_aov_tab[2, "Df"]))
  lag_group_p <- unname(as.numeric(lag_aov_tab[1, "Pr(>F)"]))
  
  lag_pairwise <- as.data.frame(lag_posthoc$genotype_plot)
  lag_pairwise$comparison <- rownames(lag_pairwise)
  rownames(lag_pairwise) <- NULL
  
  lag_pairwise <- as_tibble(lag_pairwise) %>%
    rename(
      diff      = diff,
      conf_low  = lwr,
      conf_high = upr,
      p_adj     = `p adj`
    ) %>%
    mutate(
      method = "Tukey",
      .before = 1
    )
  
} else {
  lag_model <- kruskal.test(t0_lag ~ genotype_plot, data = params_combined)
  lag_posthoc <- dunnTest(t0_lag ~ genotype_plot, data = params_combined, method = "bh")
  
  pvals_lag <- lag_posthoc$res$P.adj
  comp_names <- lag_posthoc$res$Comparison
  comp_names <- gsub(" - ", "-", comp_names)
  comp_names <- gsub(" ", "", comp_names)
  names(pvals_lag) <- comp_names
  
  lag_letters <- multcompLetters(pvals_lag)
  
  lag_letters_df <- data.frame(
    genotype_plot = names(lag_letters$Letters),
    Letter = as.character(lag_letters$Letters),
    stringsAsFactors = FALSE
  )
  
  lag_group_method <- "Kruskal-Wallis + Dunn"
  lag_group_stat_type <- "chi_squared"
  lag_group_stat <- unname(as.numeric(lag_model$statistic))
  lag_group_df1 <- unname(as.numeric(lag_model$parameter))
  lag_group_df2 <- NA_real_
  lag_group_p <- lag_model$p.value
  
  lag_pairwise <- lag_posthoc$res %>%
    as_tibble() %>%
    rename(
      comparison = Comparison,
      statistic  = Z,
      p_unadj    = P.unadj,
      p_adj      = P.adj
    ) %>%
    mutate(
      method = "Dunn",
      .before = 1
    )
}

# ----------------------------
# 12. Label positions for lag letters
# ----------------------------
lag_span <- diff(range(params_combined$t0_lag, na.rm = TRUE))
if (lag_span == 0) lag_span <- 1

label_lag <- params_combined %>%
  group_by(genotype_plot) %>%
  summarise(
    y_max = max(t0_lag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lag_letters_df, by = "genotype_plot") %>%
  mutate(
    label_y = max(params_combined$t0_lag, na.rm = TRUE) + 0.08 * lag_span
  )

# ----------------------------
# 13. Lag boxplot
# ----------------------------
p_lag <- ggplot(
  params_combined,
  aes(x = genotype_plot, y = t0_lag)
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.6,
    linetype = "dashed",
    color = "black"
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
    x = "Colony",
    y = "Polyp \u2212 Area inflection lag (days)"
  ) +
  coord_cartesian(
    ylim = c(
      min(params_combined$t0_lag, na.rm = TRUE) - 0.05 * lag_span,
      max(label_lag$label_y, na.rm = TRUE) + 0.08 * lag_span
    ),
    clip = "off"
  ) +
  lag_theme

# ----------------------------
# 14. Remove legends from top panels
# ----------------------------
p_L_noleg  <- p_L  + theme(legend.position = "none")
p_k_noleg  <- p_k  + theme(legend.position = "none")
p_t0_noleg <- p_t0 + theme(legend.position = "none")

# ----------------------------
# 15. Assemble 2x2 figure
# ----------------------------
figure_main_2x2 <- plot_grid(
  plot_grid(
    p_L_noleg, p_k_noleg,
    labels = c("A", "B"),
    label_size = base_size + 3,
    label_fontface = "bold",
    ncol = 2,
    align = "h"
  ),
  plot_grid(
    p_t0_noleg, p_lag,
    labels = c("C", "D"),
    label_size = base_size + 3,
    label_fontface = "bold",
    ncol = 2,
    align = "h"
  ),
  ncol = 1,
  align = "v"
)

print(figure_main_2x2)

# ----------------------------
# 16. Save figure
# ----------------------------
ggsave(
  filename = "Fig_main_2x2_coupling_lag_t0.png",
  plot     = figure_main_2x2,
  width    = 250,
  height   = 230,
  units    = "mm",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename    = "Fig_main_2x2_coupling_lag_t0.tif",
  plot        = figure_main_2x2,
  width       = 300,
  height      = 300,
  units       = "mm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

# ----------------------------
# 17. Export statistics
# ----------------------------
cross_process_correlations <- tibble(
  parameter = c("L", "k", "t0"),
  rho       = c(unname(cor_L$estimate), unname(cor_k$estimate), unname(cor_t0$estimate)),
  p_value   = c(cor_L$p.value, cor_k$p.value, cor_t0$p.value)
)

lag_omnibus <- tibble(
  analysis       = c("lag_vs_zero", "lag_among_colonies"),
  method         = c(lag_zero_method, lag_group_method),
  statistic_type = c(lag_zero_stat_type, lag_group_stat_type),
  statistic      = c(lag_zero_stat, lag_group_stat),
  df1            = c(lag_zero_df, lag_group_df1),
  df2            = c(NA_real_, lag_group_df2),
  p_value        = c(lag_zero_p, lag_group_p),
  shapiro_p      = c(lag_shapiro$p.value, lag_resid_shapiro$p.value),
  levene_p       = c(NA_real_, lag_levene[["Pr(>F)"]][1])
)

write_csv(cross_process_correlations, "Cross_process_spearman_correlations_t0.csv")
write_csv(lag_omnibus, "Lag_omnibus_tests_t0.csv")
write_csv(lag_pairwise, "Lag_pairwise_tests_t0.csv")

# ----------------------------
# 18. Console output so you can check the letters source
# ----------------------------
cat("\n==============================\n")
cat("CROSS-PROCESS CORRELATIONS\n")
cat("==============================\n")
print(cross_process_correlations)

cat("\n==============================\n")
cat("LAG VS ZERO\n")
cat("==============================\n")
cat("Method used:", lag_zero_method, "\n")
print(lag_zero_test)

cat("\n==============================\n")
cat("LAG AMONG COLONIES\n")
cat("==============================\n")
cat("Method used for letters:", lag_group_method, "\n")
cat("This means the boxplot letters come from:\n")
if (lag_use_anova) {
  cat("  -> Tukey post hoc after ANOVA\n")
} else {
  cat("  -> Dunn post hoc after Kruskal-Wallis\n")
}

cat("\nAssumption checks for colony effect on lag:\n")
cat("Residual Shapiro p =", lag_resid_shapiro$p.value, "\n")
cat("Levene p =", lag_levene[["Pr(>F)"]][1], "\n\n")

if (lag_use_anova) {
  print(summary(lag_model))
  print(lag_posthoc)
} else {
  print(lag_model)
  print(lag_posthoc$res)
}

cat("\n==============================\n")
cat("LAG LETTERS ACTUALLY PLOTTED\n")
cat("==============================\n")
print(lag_letters_df)


k_diff <- params_combined$k_area - params_combined$k_poly
shapiro.test(k_diff)

# robust choice (safe)
wilcox.test(params_combined$k_area,
            params_combined$k_poly,
            paired = TRUE)
