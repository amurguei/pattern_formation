## ============================================================
## Gompertz parameter comparisons among genotypes
## - safe internal IDs (G01, G02, ...)
## - display labels kept as SC1, SC2, ...
## - ANOVA + Tukey when assumptions are met
## - Kruskal-Wallis + Dunn (BH) otherwise
## - boxplots + jitter + means + significance letters
## - t0 displayed as t0 subscript in labels
## - delta removed from area labels
## ============================================================

# 0. Working directory ------------------------------------------------------
setwd("/Users/amalia/Documents/GitHub/pattern_formation/gompertz")

# 1. Packages ---------------------------------------------------------------
library(tidyverse)
library(readr)
library(car)           # leveneTest
library(multcompView)  # multcompLetters, multcompLetters4
library(FSA)           # dunnTest
library(cowplot)

# 2. Load data --------------------------------------------------------------
polyps <- read_csv("gompertz_individual_fit_parameters_polyps.csv")
area   <- read_csv("gompertz_individual_fit_parameters_area.csv")

# 3. Prepare safe IDs + display labels -------------------------------------
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

polyps <- prepare_gompertz_data(polyps)
area   <- prepare_gompertz_data(area)

# 4. Theme and colors -------------------------------------------------------
base_size <- 13

common_box_theme <- theme_minimal(base_size = base_size) +
  theme(
    plot.title         = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    axis.title.x       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.background    = element_rect(fill = "white", color = NA),
    plot.margin        = margin(t = 14, r = 10, b = 8, l = 10)
  )

fill_cadet <- function(levels_vec) {
  scale_fill_manual(
    values = setNames(rep("cadetblue3", length(levels_vec)), levels_vec),
    drop = FALSE
  )
}

# 5. Assumption checks ------------------------------------------------------
check_assumptions <- function(df, response, group = "group_id") {
  formula_obj <- as.formula(paste(response, "~", group))
  lm_obj      <- lm(formula_obj, data = df)
  
  shapiro_p <- tryCatch(
    shapiro.test(residuals(lm_obj))$p.value,
    error = function(e) NA_real_
  )
  
  levene_p <- tryCatch(
    car::leveneTest(formula_obj, data = df)[["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  
  tibble(
    parameter   = response,
    shapiro_p   = shapiro_p,
    levene_p    = levene_p,
    normal_ok   = !is.na(shapiro_p) && shapiro_p > 0.05,
    variance_ok = !is.na(levene_p)  && levene_p  > 0.05,
    use_anova   = normal_ok && variance_ok
  )
}

# 6. Statistics + compact letters ------------------------------------------
run_group_test <- function(df, response, group = "group_id") {
  assump <- check_assumptions(df, response, group)
  formula_obj <- as.formula(paste(response, "~", group))
  
  if (assump$use_anova) {
    model_obj <- aov(formula_obj, data = df)
    tuk       <- TukeyHSD(model_obj)
    let       <- multcompLetters4(model_obj, tuk)
    
    letters_df <- tibble(
      group_id = names(let[[group]]$Letters),
      Letter   = as.character(let[[group]]$Letters)
    )
    
    omnibus_p <- summary(model_obj)[[1]][["Pr(>F)"]][1]
    
    list(
      method      = "ANOVA + Tukey",
      omnibus_p   = omnibus_p,
      model       = model_obj,
      posthoc     = tuk,
      letters_df  = letters_df,
      assumptions = assump
    )
    
  } else {
    kw   <- kruskal.test(formula_obj, data = df)
    dunn <- dunnTest(formula_obj, data = df, method = "bh")
    
    pvals <- dunn$res$P.adj
    comp_names <- dunn$res$Comparison
    
    comp_names <- gsub(" - ", "-", comp_names)
    comp_names <- gsub(" ", "", comp_names)
    names(pvals) <- comp_names
    
    let <- multcompLetters(pvals)
    
    letters_df <- tibble(
      group_id = names(let$Letters),
      Letter   = as.character(let$Letters)
    )
    
    list(
      method      = "Kruskal-Wallis + Dunn",
      omnibus_p   = kw$p.value,
      model       = kw,
      posthoc     = dunn,
      letters_df  = letters_df,
      assumptions = assump
    )
  }
}

# 7. Label positions --------------------------------------------------------
make_label_data <- function(df, letters_df, value_col) {
  
  ydata <- df %>%
    group_by(group_id, genotype_plot) %>%
    summarise(
      y_max = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  out <- ydata %>%
    left_join(letters_df, by = "group_id")
  
  out
}

# 8. Plot function ----------------------------------------------------------
make_gompertz_plot <- function(df, response, ylab, title = NULL,
                               point_fill = "red") {
  
  test_res <- run_group_test(df, response, group = "group_id")
  
  label_df <- make_label_data(
    df = df,
    letters_df = test_res$letters_df,
    value_col = response
  )
  
  y_min_data <- min(df[[response]], na.rm = TRUE)
  y_max_data <- max(df[[response]], na.rm = TRUE)
  y_span     <- y_max_data - y_min_data
  
  if (y_span == 0) y_span <- 1
  
  # Put all significance letters on one aligned top line
  base_label_y <- y_max_data + 0.10 * y_span
  
  label_df <- label_df %>%
    mutate(label_y = base_label_y)
  
  y_upper <- base_label_y + 0.08 * y_span
  y_lower <- y_min_data - 0.05 * y_span
  
  p <- ggplot(df, aes(x = genotype_plot, y = .data[[response]], fill = genotype_plot)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.5) +
    geom_jitter(width = 0.12, alpha = 0.7, size = 2, color = "black") +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 21,
      size = 3.2,
      fill = point_fill,
      color = "black",
      stroke = 0.5
    ) +
    geom_text(
      data = label_df,
      aes(x = genotype_plot, y = label_y, label = Letter),
      inherit.aes = FALSE,
      vjust = 0,
      size = 4
    ) +
    labs(
      title = title,
      x = "Colony",
      y = ylab
    ) +
    fill_cadet(levels(df$genotype_plot)) +
    coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
    common_box_theme
  
  list(
    plot   = p,
    test   = test_res,
    labels = label_df
  )
}

## 10. Build plots: POLYPS --------------------------------------------------

poly_L <- make_gompertz_plot(
  polyps,
  "L",
  ylab  = "Asymptotic polyp number (L)",
  title = "Polyp Gompertz parameter: L"
)

poly_k <- make_gompertz_plot(
  polyps,
  "k",
  ylab  = "Growth rate (k)",
  title = "Polyp Gompertz parameter: k"
)

poly_x0 <- make_gompertz_plot(
  polyps,
  "x0",
  ylab  = expression(bold("Inflection day (" * t[0] * ")")),
  title = "Polyp Gompertz parameter: t0"
)

## 11. Build plots: AREA ----------------------------------------------------

area_L <- make_gompertz_plot(
  area,
  "L",
  ylab  = "Asymptotic area (mm²)",
  title = "Area Gompertz parameter: L"
)

area_k <- make_gompertz_plot(
  area,
  "k",
  ylab  = "Growth rate (k)",
  title = "Area Gompertz parameter: k"
)

area_x0 <- make_gompertz_plot(
  area,
  "x0",
  ylab  = expression(bold("Inflection day (" * t[0] * ")")),
  title = "Area Gompertz parameter: t0"
)

## ============================================================
## Export Gompertz statistical results to CSV
## ============================================================

extract_test_info <- function(obj, dataset_name, param_name) {
  
  test <- obj$test
  
  tibble(
    dataset   = dataset_name,
    parameter = param_name,
    method    = test$method,
    p_value   = test$omnibus_p,
    
    statistic = if (!is.null(test$model$statistic))
      as.numeric(test$model$statistic)
    else NA_real_,
    
    df = if (!is.null(test$model$parameter))
      as.numeric(test$model$parameter)
    else NA_real_
  )
}

gompertz_results <- bind_rows(
  extract_test_info(poly_L,  "polyps", "L"),
  extract_test_info(poly_k,  "polyps", "k"),
  extract_test_info(poly_x0, "polyps", "x0"),
  extract_test_info(area_L,  "area", "L"),
  extract_test_info(area_k,  "area", "k"),
  extract_test_info(area_x0, "area", "x0")
)

print(gompertz_results)

write_csv(
  gompertz_results,
  "Gompertz_parameter_statistics.csv"
)

# 12. Arrange figures -------------------------------------------------------
figure_polyps <- plot_grid(
  poly_L$plot, poly_k$plot, poly_x0$plot,
  labels = c("A", "B", "C"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h"
)

figure_area <- plot_grid(
  area_L$plot, area_k$plot, area_x0$plot,
  labels = c("A", "B", "C"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h"
)

print(figure_polyps)
print(figure_area)

# 13. Save ------------------------------------------------------------------
out_dir <- file.path(getwd(), "gompertz_parameter_figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = file.path(out_dir, "Fig_polyps_gompertz_parameters.png"),
  plot = figure_polyps,
  width = 300, height = 110, units = "mm",
  dpi = 600, bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig_polyps_gompertz_parameters.tif"),
  plot = figure_polyps,
  width = 300, height = 110, units = "mm",
  dpi = 600, compression = "lzw", bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig_area_gompertz_parameters.png"),
  plot = figure_area,
  width = 300, height = 110, units = "mm",
  dpi = 600, bg = "white"
)

ggsave(
  filename = file.path(out_dir, "Fig_area_gompertz_parameters.tif"),
  plot = figure_area,
  width = 300, height = 110, units = "mm",
  dpi = 600, compression = "lzw", bg = "white"
)

# 14. Summary table ---------------------------------------------------------
extract_summary_line <- function(res, dataset_name, parameter_name) {
  tibble(
    dataset   = dataset_name,
    parameter = parameter_name,
    method    = res$method,
    omnibus_p = res$omnibus_p,
    shapiro_p = res$assumptions$shapiro_p,
    levene_p  = res$assumptions$levene_p
  )
}

stats_summary <- bind_rows(
  extract_summary_line(poly_L$test,  "polyps", "L"),
  extract_summary_line(poly_k$test,  "polyps", "k"),
  extract_summary_line(poly_x0$test, "polyps", "x0"),
  extract_summary_line(area_L$test,  "area",   "L"),
  extract_summary_line(area_k$test,  "area",   "k"),
  extract_summary_line(area_x0$test, "area",   "x0")
)

print(stats_summary)
write_csv(stats_summary, file.path(out_dir, "gompertz_parameter_stats_summary.csv"))

# 15. Quick checks ----------------------------------------------------------
cat("\n--- Methods chosen ---\n")
cat("Polyps L :", poly_L$test$method, "\n")
cat("Polyps k :", poly_k$test$method, "\n")
cat("Polyps x0:", poly_x0$test$method, "\n")
cat("Area L   :", area_L$test$method, "\n")
cat("Area k   :", area_k$test$method, "\n")
cat("Area x0  :", area_x0$test$method, "\n")

cat("\n--- Letters for polyp L ---\n")
print(poly_L$test$letters_df)
print(poly_L$labels)

## ============================================================
## Paneled boxplots: AREA (top) + POLYPS (bottom)
## ============================================================

# --- AREA row (top) ---
area_row <- plot_grid(
  area_L$plot,
  area_k$plot,
  area_x0$plot,
  labels = c("(a)", "(b)", "(c)"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h"
)

# --- POLYPS row (bottom) ---
poly_row <- plot_grid(
  poly_L$plot,
  poly_k$plot,
  poly_x0$plot,
  labels = c("(d)", "(e)", "(f)"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h"
)

# --- Combine into one figure ---
figure_combined <- plot_grid(
  area_row,
  poly_row,
  ncol = 1,
  align = "v"
)

print(figure_combined)

## Save combined boxplot figure ---------------------------------------------

ggsave(
  filename = "Gompertz_boxplots_area_top_polyps_bottom_mod.png",
  plot     = figure_combined,
  width    = 330,
  height   = 200,
  units    = "mm",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename    = "Gompertz_boxplots_area_top_polyps_bottom_mod.tif",
  plot        = figure_combined,
  width       = 330,
  height      = 200,
  units       = "mm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

getwd()
