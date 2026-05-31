## ============================================================
## Gompertz parameter comparisons among genotypes
## - safe internal IDs (G01, G02, ...)
## - display labels kept as SC1, SC2, ...
## - ANOVA + Tukey when assumptions are met
## - Kruskal-Wallis + Dunn (BH) otherwise
## - boxplots + jitter + means + significance letters
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
  ylab  = "Inflection day (t0)",
  title = "Polyp Gompertz parameter: t0"
)

## 11. Build plots: AREA ----------------------------------------------------

area_L <- make_gompertz_plot(
  area,
  "L",
  ylab  = "Asymptotic Δ area (mm²)",
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
  ylab  = "Inflection day (t0)",
  title = "Area Gompertz parameter: t0"
)


## ============================================================
## Export Gompertz statistical results to CSV
## ============================================================

library(tidyverse)

# ---- Helper function to extract info safely ------------------
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

# ---- Build results table ------------------------------------
gompertz_results <- bind_rows(
  
  extract_test_info(poly_L,  "polyps", "L"),
  extract_test_info(poly_k,  "polyps", "k"),
  extract_test_info(poly_x0, "polyps", "x0"),
  
  extract_test_info(area_L,  "area", "L"),
  extract_test_info(area_k,  "area", "k"),
  extract_test_info(area_x0, "area", "x0")
)

# ---- View in console ----------------------------------------
print(gompertz_results)

# ---- Save to CSV --------------------------------------------
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

library(cowplot)

# --- AREA row (top) ---
area_row <- plot_grid(
  area_L$plot,
  area_k$plot,
  area_x0$plot,
  labels = c("A", "B", "C"),
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
  labels = c("D", "E", "F"),
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
  filename = "Gompertz_boxplots_area_top_polyps_bottom.png",
  plot     = figure_combined,
  width    = 330,
  height   = 200,
  units    = "mm",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename    = "Gompertz_boxplots_area_top_polyps_bottom.tif",
  plot        = figure_combined,
  width       = 200,
  height      = 120,
  units       = "mm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

getwd()
## ============================================================
## Export omnibus + pairwise Gompertz statistics to CSV
## ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Omnibus extractor
# ------------------------------------------------------------
extract_omnibus <- function(obj, dataset_name, param_name) {
  
  test <- obj$test
  
  if (grepl("Kruskal", test$method)) {
    
    tibble(
      dataset        = dataset_name,
      parameter      = param_name,
      method         = test$method,
      statistic_type = "chi_squared",
      statistic      = unname(as.numeric(test$model$statistic)),
      df1            = unname(as.numeric(test$model$parameter)),
      df2            = NA_real_,
      p_value        = test$omnibus_p
    )
    
  } else if (grepl("ANOVA", test$method)) {
    
    aov_tab <- summary(test$model)[[1]]
    
    tibble(
      dataset        = dataset_name,
      parameter      = param_name,
      method         = test$method,
      statistic_type = "F",
      statistic      = unname(as.numeric(aov_tab[1, "F value"])),
      df1            = unname(as.numeric(aov_tab[1, "Df"])),
      df2            = unname(as.numeric(aov_tab[2, "Df"])),
      p_value        = unname(as.numeric(aov_tab[1, "Pr(>F)"]))
    )
    
  } else {
    
    tibble(
      dataset        = dataset_name,
      parameter      = param_name,
      method         = test$method,
      statistic_type = NA_character_,
      statistic      = NA_real_,
      df1            = NA_real_,
      df2            = NA_real_,
      p_value        = test$omnibus_p
    )
  }
}

# ------------------------------------------------------------
# 2. Pairwise extractor
# ------------------------------------------------------------
extract_pairwise <- function(obj, dataset_name, param_name) {
  
  test <- obj$test
  
  if (grepl("Kruskal", test$method)) {
    
    test$posthoc$res %>%
      as_tibble() %>%
      rename(
        comparison = Comparison,
        statistic  = Z,
        p_unadj    = P.unadj,
        p_adj      = P.adj
      ) %>%
      mutate(
        dataset   = dataset_name,
        parameter = param_name,
        method    = "Dunn",
        .before = 1
      )
    
  } else if (grepl("ANOVA", test$method)) {
    
    tuk_list <- test$posthoc
    term_name <- names(tuk_list)[1]
    
    tuk_df <- as.data.frame(tuk_list[[term_name]])
    tuk_df$comparison <- rownames(tuk_df)
    rownames(tuk_df) <- NULL
    
    as_tibble(tuk_df) %>%
      rename(
        diff      = diff,
        conf_low  = lwr,
        conf_high = upr,
        p_adj     = `p adj`
      ) %>%
      mutate(
        dataset   = dataset_name,
        parameter = param_name,
        method    = "Tukey",
        .before = 1
      )
    
  } else {
    
    tibble(
      dataset   = dataset_name,
      parameter = param_name,
      method    = NA_character_
    )
  }
}

# ------------------------------------------------------------
# 3. Build omnibus table
# ------------------------------------------------------------
gompertz_omnibus <- bind_rows(
  extract_omnibus(poly_L,  "polyps", "L"),
  extract_omnibus(poly_k,  "polyps", "k"),
  extract_omnibus(poly_x0, "polyps", "x0"),
  extract_omnibus(area_L,  "area",   "L"),
  extract_omnibus(area_k,  "area",   "k"),
  extract_omnibus(area_x0, "area",   "x0")
)

print(gompertz_omnibus)

# ------------------------------------------------------------
# 4. Build pairwise table
# ------------------------------------------------------------
gompertz_pairwise <- bind_rows(
  extract_pairwise(poly_L,  "polyps", "L"),
  extract_pairwise(poly_k,  "polyps", "k"),
  extract_pairwise(poly_x0, "polyps", "x0"),
  extract_pairwise(area_L,  "area",   "L"),
  extract_pairwise(area_k,  "area",   "k"),
  extract_pairwise(area_x0, "area",   "x0")
)

print(gompertz_pairwise)

# ------------------------------------------------------------
# 5. Save to CSV
# ------------------------------------------------------------
write_csv(gompertz_omnibus,  "Gompertz_parameter_omnibus_tests.csv")
write_csv(gompertz_pairwise, "Gompertz_parameter_pairwise_tests.csv")
getwd()

colnames(polyps)
colnames(area)

params_poly <- polyps %>%
  select(genotype_plot, Replicate, L_poly = L, k_poly = k, x0_poly = x0)

params_area <- area %>%
  select(genotype_plot, Replicate, L_area = L, k_area = k, x0_area = x0)

params_combined <- left_join(
  params_area,
  params_poly,
  by = c("genotype_plot", "Replicate")
)

print(params_combined)

polyps %>% count(genotype_plot, Replicate)
area %>% count(genotype_plot, Replicate)

anti_join(
  area %>% select(genotype_plot, Replicate),
  polyps %>% select(genotype_plot, Replicate),
  by = c("genotype_plot", "Replicate")
)

anti_join(
  polyps %>% select(genotype_plot, Replicate),
  area %>% select(genotype_plot, Replicate),
  by = c("genotype_plot", "Replicate")
)


## ============================================================
## Cross-process parameter relationships (Area ↔ Polyps)
## ============================================================

library(ggplot2)
library(cowplot)

base_size <- 14

p_L <- ggplot(params_combined,
              aes(x = L_area, y = L_poly)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    x = "Area asymptote (L)",
    y = "Polyp asymptote (L)"
  ) +
  theme_minimal(base_size = base_size)

p_k <- ggplot(params_combined,
              aes(x = k_area, y = k_poly)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    x = "Area growth rate (k)",
    y = "Polyp growth rate (k)"
  ) +
  theme_minimal(base_size = base_size)

p_x0 <- ggplot(params_combined,
               aes(x = x0_area, y = x0_poly)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    x = "Area inflection day (x0)",
    y = "Polyp inflection day (x0)"
  ) +
  theme_minimal(base_size = base_size)

figure_relations <- plot_grid(
  p_L, p_k, p_x0,
  labels = c("A", "B", "C"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h",
)

print(figure_relations)

## ============================================================
## Reusable genotype colors (Okabe–Ito style)
## ============================================================

# Core Okabe-Ito palette
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

# Fixed assignments by colony label
# Adjust these if you want to preserve an older exact mapping
fixed_cols <- c(
  "SC1" = "#0072B2", # blue
  "SC2" = "#E69F00", # orange
  "SC5" = "#009E73"  # bluish green
)

# Build full mapping from the levels already present in your data
geno_levels <- levels(polyps$genotype_plot)

remaining_levels <- setdiff(geno_levels, names(fixed_cols))
remaining_pool   <- setdiff(okabe_ito_base, unname(fixed_cols))

# If needed, extend palette safely
if (length(remaining_levels) > length(remaining_pool)) {
  remaining_pool <- c(
    remaining_pool,
    "#56B4E9", # reuse if necessary
    "#D55E00",
    "#CC79A7",
    "#999999"
  )
}

other_cols <- setNames(remaining_pool[seq_along(remaining_levels)], remaining_levels)

geno_colors <- c(fixed_cols, other_cols)
geno_colors <- geno_colors[geno_levels]  # preserve plotting order

print(geno_colors)

## ============================================================
## Cross-process parameter relationships (Area ↔ Polyps)
## Colored by colony using consistent Okabe-Ito mapping
## ============================================================

library(ggplot2)
library(cowplot)

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

p_L <- ggplot(
  params_combined,
  aes(x = L_area, y = L_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Asymptotic Δ area (L, mm²)",
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
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Area growth rate (k)",
    y = "Polyp growth rate (k)",
    color = "Colony"
  ) +
  common_scatter_theme

p_x0 <- ggplot(
  params_combined,
  aes(x = x0_area, y = x0_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Area inflection day (x0)",
    y = "Polyp inflection day (x0)",
    color = "Colony"
  ) +
  common_scatter_theme

legend_shared <- get_legend(
  p_L + theme(legend.position = "bottom")
)

p_L_noleg  <- p_L  + theme(legend.position = "none")
p_k_noleg  <- p_k  + theme(legend.position = "none")
p_x0_noleg <- p_x0 + theme(legend.position = "none")

figure_relations_colored <- plot_grid(
  plot_grid(
    p_L_noleg, p_k_noleg, p_x0_noleg,
    labels = c("A", "B", "C"),
    label_size = base_size + 3,
    label_fontface = "bold",
    ncol = 3,
    align = "h"
  ),
  legend_shared,
  ncol = 1,
  rel_heights = c(1, 0.14)
)

print(figure_relations_colored)


## ============================================================
## Correlation summaries for Area ↔ Polyps parameter relationships
## ============================================================

cor_L  <- cor.test(params_combined$L_area,  params_combined$L_poly,  method = "spearman", exact = FALSE)
cor_k  <- cor.test(params_combined$k_area,  params_combined$k_poly,  method = "spearman", exact = FALSE)
cor_x0 <- cor.test(params_combined$x0_area, params_combined$x0_poly, method = "spearman", exact = FALSE)

cor_labels <- list(
  L  = paste0("Spearman \u03C1 = ", round(unname(cor_L$estimate), 2),
              "\np = ", format.pval(cor_L$p.value, digits = 2, eps = 1e-4)),
  k  = paste0("Spearman \u03C1 = ", round(unname(cor_k$estimate), 2),
              "\np = ", format.pval(cor_k$p.value, digits = 2, eps = 1e-4)),
  x0 = paste0("Spearman \u03C1 = ", round(unname(cor_x0$estimate), 2),
              "\np = ", format.pval(cor_x0$p.value, digits = 2, eps = 1e-4))
)

base_size <- 14

common_scatter_theme <- theme_minimal(base_size = base_size) +
  theme(
    plot.title         = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    axis.title.x       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(size = base_size, face = "bold"),
    legend.text        = element_text(size = base_size - 1),
    plot.background    = element_rect(fill = "white", color = NA)
  )

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

p_x0 <- ggplot(
  params_combined,
  aes(x = x0_area, y = x0_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels$x0,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Area inflection day (x0)",
    y = "Polyp inflection day (x0)",
    color = "Colony"
  ) +
  common_scatter_theme

legend_shared <- cowplot::get_legend(
  p_L + theme(legend.position = "bottom")
)

p_L_noleg  <- p_L  + theme(legend.position = "none")
p_k_noleg  <- p_k  + theme(legend.position = "none")
p_x0_noleg <- p_x0 + theme(legend.position = "none")

figure_relations_colored <- plot_grid(
  plot_grid(
    p_L_noleg, p_k_noleg, p_x0_noleg,
    labels = c("A", "B", "C"),
    label_size = base_size + 3,
    label_fontface = "bold",
    ncol = 3,
    align = "h"
  ),
  legend_shared,
  ncol = 1,
  rel_heights = c(1, 0.14)
)

print(figure_relations_colored)

ggsave(
  filename = "Gompertz_parameter_relationships_okabeito.png",
  plot     = figure_relations_colored,
  width    = 280,
  height   = 140,
  units    = "mm",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename    = "Gompertz_parameter_relationships_okabeito.tif",
  plot        = figure_relations_colored,
  width       = 280,
  height      = 140,
  units       = "mm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)


## ============================================================
## Inflection lag boxplot: x0_poly - x0_area
## Positive values = area inflects earlier than polyps
## ============================================================

params_combined <- params_combined %>%
  mutate(x0_lag = x0_poly - x0_area)

base_size <- 14

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

p_lag <- ggplot(
  params_combined,
  aes(x = genotype_plot, y = x0_lag, fill = genotype_plot)
) +
  geom_hline(yintercept = 0, linewidth = 0.6, linetype = "dashed", color = "black") +
  geom_boxplot(outlier.shape = NA, width = 0.65, linewidth = 0.5) +
  geom_jitter(aes(color = genotype_plot), width = 0.12, alpha = 0.8, size = 2) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.2,
    fill = "white",
    color = "black",
    stroke = 0.7
  ) +
  scale_fill_manual(values = geno_colors, drop = FALSE) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Colony",
    y = "Polyp − Area inflection lag (days)"
  ) +
  lag_theme

print(p_lag)

ggplot(params_combined,
       aes(x = genotype_plot, y = x0_lag)) +
  
  # Boxplots
  geom_boxplot(
    fill = "cadetblue3",
    color = "black",
    width = 0.7,
    outlier.shape = NA
  ) +
  
  # Jittered individual colonies (optional but recommended)
  geom_jitter(
    width = 0.12,
    alpha = 0.7,
    size = 2
  ) +
  
  # Mean points (RED)
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 21,
    size = 3.5,
    fill = "red",
    color = "black",
    stroke = 0.5
  ) +
  
  labs(
    x = "Colony",
    y = "Polyp − Area inflection lag (days)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


## ============================================================
## Within-process parameter relationships
## Area (top):    L vs k ; x0 vs k
## Polyps (bottom): L vs k ; x0 vs k
## ============================================================

base_size <- 14

common_scatter_theme <- theme_minimal(base_size = base_size) +
  theme(
    plot.title         = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    axis.title.x       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(size = base_size, face = "bold"),
    legend.text        = element_text(size = base_size - 1),
    plot.background    = element_rect(fill = "white", color = NA)
  )

p_area_Lk <- ggplot(
  params_combined,
  aes(x = L_area, y = k_area, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels_tradeoff$area_Lk,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Asymptotic Δ area L (mm²)",
    y = "Area growth rate (k)",
    color = "Colony"
  ) +
  common_scatter_theme

p_area_x0k <- ggplot(
  params_combined,
  aes(x = x0_area, y = k_area, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels_tradeoff$area_x0k,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Area inflection day (x0)",
    y = "Area growth rate (k)",
    color = "Colony"
  ) +
  common_scatter_theme

p_poly_Lk <- ggplot(
  params_combined,
  aes(x = L_poly, y = k_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels_tradeoff$poly_Lk,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Asymptotic polyp number (L)",
    y = "Polyp growth rate (k)",
    color = "Colony"
  ) +
  common_scatter_theme

p_poly_x0k <- ggplot(
  params_combined,
  aes(x = x0_poly, y = k_poly, color = genotype_plot)
) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = cor_labels_tradeoff$poly_x0k,
    hjust = -0.1, vjust = 1.1,
    size = 4.2
  ) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(
    x = "Polyp inflection day (x0)",
    y = "Polyp growth rate (k)",
    color = "Colony"
  ) +
  common_scatter_theme

legend_tradeoff <- cowplot::get_legend(
  p_area_Lk + theme(legend.position = "bottom")
)

p_area_Lk_noleg  <- p_area_Lk  + theme(legend.position = "none")
p_area_x0k_noleg <- p_area_x0k + theme(legend.position = "none")
p_poly_Lk_noleg  <- p_poly_Lk  + theme(legend.position = "none")
p_poly_x0k_noleg <- p_poly_x0k + theme(legend.position = "none")

top_row_tradeoff <- plot_grid(
  p_area_Lk_noleg, p_area_x0k_noleg,
  labels = c("A", "B"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 2,
  align = "h"
)

bottom_row_tradeoff <- plot_grid(
  p_poly_Lk_noleg, p_poly_x0k_noleg,
  labels = c("C", "D"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 2,
  align = "h"
)

figure_tradeoff <- plot_grid(
  plot_grid(top_row_tradeoff, bottom_row_tradeoff, ncol = 1, align = "v"),
  legend_tradeoff,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

print(figure_tradeoff)

# --- Save within-process trade-off figure ---

out_width_mm  <- 200   # good for single-column or compact 2-column layout
out_height_mm <- 200
out_dpi       <- 600

ggsave(
  filename = "Gompertz_parameter_tradeoffs_okabeito.png",
  plot     = figure_tradeoff,
  width    = out_width_mm,
  height   = out_height_mm,
  units    = "mm",
  dpi      = out_dpi,
  bg       = "white"
)

ggsave(
  filename    = "Gompertz_parameter_tradeoffs_okabeito.tif",
  plot        = figure_tradeoff,
  width       = out_width_mm,
  height      = out_height_mm,
  units       = "mm",
  dpi         = out_dpi,
  compression = "lzw",
  bg          = "white"
)

## ============================================================
## Within-process parameter relationships including L vs x0
## Area (top):    L vs k ; x0 vs k ; x0 vs L
## Polyps (bottom): L vs k ; x0 vs k ; x0 vs L
## Spearman correlations shown on panels
## Uses: params_combined, geno_colors
## ============================================================

library(tidyverse)
library(ggplot2)
library(cowplot)

base_size <- 14

# ------------------------------------------------------------
# 1. Spearman correlations
# ------------------------------------------------------------
cor_area_Lk  <- cor.test(params_combined$L_area,  params_combined$k_area,
                         method = "spearman", exact = FALSE)

cor_area_x0k <- cor.test(params_combined$x0_area, params_combined$k_area,
                         method = "spearman", exact = FALSE)

cor_area_Lx0 <- cor.test(params_combined$L_area,  params_combined$x0_area,
                         method = "spearman", exact = FALSE)

cor_poly_Lk  <- cor.test(params_combined$L_poly,  params_combined$k_poly,
                         method = "spearman", exact = FALSE)

cor_poly_x0k <- cor.test(params_combined$x0_poly, params_combined$k_poly,
                         method = "spearman", exact = FALSE)

cor_poly_Lx0 <- cor.test(params_combined$L_poly,  params_combined$x0_poly,
                         method = "spearman", exact = FALSE)

make_cor_label <- function(ct) {
  paste0(
    "Spearman \u03C1 = ",
    round(unname(ct$estimate), 2),
    "\np = ",
    format.pval(ct$p.value, digits = 2, eps = 1e-4)
  )
}

cor_labels_tradeoff <- list(
  area_Lk  = make_cor_label(cor_area_Lk),
  area_x0k = make_cor_label(cor_area_x0k),
  area_Lx0 = make_cor_label(cor_area_Lx0),
  poly_Lk  = make_cor_label(cor_poly_Lk),
  poly_x0k = make_cor_label(cor_poly_x0k),
  poly_Lx0 = make_cor_label(cor_poly_Lx0)
)

# ------------------------------------------------------------
# 2. Theme
# ------------------------------------------------------------
common_scatter_theme <- theme_minimal(base_size = base_size) +
  theme(
    plot.title         = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    axis.title.x       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.title.y       = element_text(size = base_size + 1, face = "bold", colour = "black"),
    axis.text          = element_text(size = base_size - 1, colour = "black"),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(size = base_size, face = "bold"),
    legend.text        = element_text(size = base_size - 1),
    plot.background    = element_rect(fill = "white", color = NA)
  )

# ------------------------------------------------------------
# 3. Plot helper
# ------------------------------------------------------------
make_tradeoff_plot <- function(df, xvar, yvar, xlab, ylab, label_text) {
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]], color = genotype_plot)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = label_text,
      hjust = -0.1, vjust = 1.1,
      size = 4.2
    ) +
    scale_color_manual(values = geno_colors, drop = FALSE) +
    labs(
      x = xlab,
      y = ylab,
      color = "Colony"
    ) +
    common_scatter_theme
}

# ------------------------------------------------------------
# 4. Area panels
# ------------------------------------------------------------
p_area_Lk <- make_tradeoff_plot(
  params_combined,
  xvar = "L_area",
  yvar = "k_area",
  xlab = "Asymptotic \u0394 area (mm\u00B2) L",
  ylab = "Area growth rate (k)",
  label_text = cor_labels_tradeoff$area_Lk
)

p_area_x0k <- make_tradeoff_plot(
  params_combined,
  xvar = "x0_area",
  yvar = "k_area",
  xlab = "Area inflection day (x0)",
  ylab = "Area growth rate (k)",
  label_text = cor_labels_tradeoff$area_x0k
)

p_area_Lx0 <- make_tradeoff_plot(
  params_combined,
  xvar = "x0_area",
  yvar = "L_area",
  xlab = "Area inflection day (x0)",
  ylab = "Asymptotic \u0394 area (mm\u00B2) L",
  label_text = cor_labels_tradeoff$area_Lx0
)

# ------------------------------------------------------------
# 5. Polyp panels
# ------------------------------------------------------------
p_poly_Lk <- make_tradeoff_plot(
  params_combined,
  xvar = "L_poly",
  yvar = "k_poly",
  xlab = "Asymptotic polyp number (L)",
  ylab = "Polyp growth rate (k)",
  label_text = cor_labels_tradeoff$poly_Lk
)

p_poly_x0k <- make_tradeoff_plot(
  params_combined,
  xvar = "x0_poly",
  yvar = "k_poly",
  xlab = "Polyp inflection day (x0)",
  ylab = "Polyp growth rate (k)",
  label_text = cor_labels_tradeoff$poly_x0k
)

p_poly_Lx0 <- make_tradeoff_plot(
  params_combined,
  xvar = "x0_poly",
  yvar = "L_poly",
  xlab = "Polyp inflection day (x0)",
  ylab = "Asymptotic polyp number (L)",
  label_text = cor_labels_tradeoff$poly_Lx0
)

# ------------------------------------------------------------
# 6. Shared legend
# ------------------------------------------------------------
legend_tradeoff <- cowplot::get_legend(
  p_area_Lk + theme(legend.position = "bottom")
)

p_area_Lk_noleg  <- p_area_Lk  + theme(legend.position = "none")
p_area_x0k_noleg <- p_area_x0k + theme(legend.position = "none")
p_area_Lx0_noleg <- p_area_Lx0 + theme(legend.position = "none")

p_poly_Lk_noleg  <- p_poly_Lk  + theme(legend.position = "none")
p_poly_x0k_noleg <- p_poly_x0k + theme(legend.position = "none")
p_poly_Lx0_noleg <- p_poly_Lx0 + theme(legend.position = "none")

# ------------------------------------------------------------
# 7. Assemble figure
# ------------------------------------------------------------
top_row_tradeoff <- plot_grid(
  p_area_Lk_noleg, p_area_x0k_noleg, p_area_Lx0_noleg,
  labels = c("A", "B", "C"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h"
)

bottom_row_tradeoff <- plot_grid(
  p_poly_Lk_noleg, p_poly_x0k_noleg, p_poly_Lx0_noleg,
  labels = c("D", "E", "F"),
  label_size = base_size + 3,
  label_fontface = "bold",
  ncol = 3,
  align = "h"
)

figure_tradeoff_full <- plot_grid(
  plot_grid(top_row_tradeoff, bottom_row_tradeoff, ncol = 1, align = "v"),
  legend_tradeoff,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

print(figure_tradeoff_full)

# ------------------------------------------------------------
# 8. Save figure
# ------------------------------------------------------------
out_width_mm  <- 280
out_height_mm <- 230
out_dpi       <- 600

ggsave(
  filename = "Gompertz_within_process_relationships_full_okabeito.png",
  plot     = figure_tradeoff_full,
  width    = out_width_mm,
  height   = out_height_mm,
  units    = "mm",
  dpi      = out_dpi,
  bg       = "white"
)

ggsave(
  filename    = "Gompertz_within_process_relationships_full_okabeito.tif",
  plot        = figure_tradeoff_full,
  width       = out_width_mm,
  height      = out_height_mm,
  units       = "mm",
  dpi         = out_dpi,
  compression = "lzw",
  bg          = "white"
)

# ------------------------------------------------------------
# 9. Export Spearman summary table
# ------------------------------------------------------------
spearman_summary <- tibble(
  dataset   = c("area", "area", "area", "polyps", "polyps", "polyps"),
  x_var     = c("L_area", "x0_area", "x0_area", "L_poly", "x0_poly", "x0_poly"),
  y_var     = c("k_area", "k_area", "L_area",  "k_poly", "k_poly", "L_poly"),
  rho       = c(
    unname(cor_area_Lk$estimate),
    unname(cor_area_x0k$estimate),
    unname(cor_area_Lx0$estimate),
    unname(cor_poly_Lk$estimate),
    unname(cor_poly_x0k$estimate),
    unname(cor_poly_Lx0$estimate)
  ),
  p_value   = c(
    cor_area_Lk$p.value,
    cor_area_x0k$p.value,
    cor_area_Lx0$p.value,
    cor_poly_Lk$p.value,
    cor_poly_x0k$p.value,
    cor_poly_Lx0$p.value
  )
)

print(spearman_summary)

write_csv(
  spearman_summary,
  "Gompertz_within_process_spearman_summary.csv"
)



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
  labels = c("A", "B", "C"),
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
  labels = c("D", "E", "F"),
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
  filename = "Gompertz_boxplots_area_top_polyps_bottom.png",
  plot     = figure_combined,
  width    = 330,
  height   = 200,
  units    = "mm",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename    = "Gompertz_boxplots_area_top_polyps_bottom.tif",
  plot        = figure_combined,
  width       = 330,
  height      = 200,
  units       = "mm",
  dpi         = 600,
  compression = "lzw",
  bg          = "white"
)

getwd()