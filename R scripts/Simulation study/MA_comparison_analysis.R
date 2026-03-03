# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : MatchAlign_comparison_analysis.R                           ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 31-01-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# R code to analyze the simulation results of the comparison between the
# MatchAlign and the efficient-RSP versions
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(SimDesign)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Auxiliar function to compute marginal means
marginal_means <- function(data, nombre_factor, cols) {
  tabla_agg <- aggregate(data[, cols], 
                         by = list(Level = data[[nombre_factor]]), 
                         FUN = mean, 
                         na.rm = TRUE)
  tabla_agg$Factor <- nombre_factor
  tabla_agg <- tabla_agg[, c("Factor", "Level", cols)]
  return(tabla_agg)
}

# Auxiliar function to put bold letters in the best metric
format_block <- function(data, cols, type = "min", digits = 3) {
  mat <- as.matrix(data[, cols])
  best_val <- if(type == "min") {
    do.call(pmin, c(as.data.frame(mat), na.rm = TRUE))
  } else {
    do.call(pmax, c(as.data.frame(mat), na.rm = TRUE))
  }
  fmt_vec <- sprintf(paste0("%.", digits, "f"), mat)
  is_best <- abs(mat - best_val) < 1e-8
  fmt_vec[is_best] <- paste0("\\textbf{", fmt_vec[is_best], "}")
  mat_out <- matrix(fmt_vec, nrow = nrow(mat), ncol = ncol(mat))
  return(mat_out)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Agreement between methods
# ─────────────────────────────────────────────────────────────────────────────

# Load simulation results
simres_s2 <- readRDS("Results/Simulation study/MatchAlign comparison/MatchAlign_comp_simres.rds")
all_res <- SimExtract(simres_s2, what = "results")

# Variables to average
vars <- c("abs_fit_base", "abs_fit_MA", "abs_fit_ERSP", 
          "rel_fit_base", "rel_fit_MA", "rel_fit_ERSP", 
          "ratio_fit_MA", "ratio_fit_ERSP", 
          "ESS_MA_bulk", "ESS_MA_tail", 
          "ESS_ERSP_bulk", "ESS_ERSP_tail", 
          "time_MA", "time_ERSP", "time_ratio")

# Compute marginal means
res_N      <- marginal_means(simres_s2, "N",      vars)
res_M      <- marginal_means(simres_s2, "M",      vars)
res_J      <- marginal_means(simres_s2, "J",      vars)
res_Lambda <- marginal_means(simres_s2, "Lambda", vars)

# Change lambda factor level names
res_Lambda$Level <- c("Independent", "Sparse")

# Bind rows and generate final data frame
df <- rbind(res_N, res_M, res_J, res_Lambda)

# Multiply relative measures per 100
df$rel_fit_base <- df$rel_fit_base * 100
df$rel_fit_MA <- df$rel_fit_MA * 100
df$rel_fit_ERSP <- df$rel_fit_ERSP * 100

# Prepare LaTeX table: bold letters in the best values
blk_abs   <- format_block(df, c("abs_fit_base", "abs_fit_MA", "abs_fit_ERSP"), "min", 2)
blk_acc   <- format_block(df, c("rel_fit_base", "rel_fit_MA", "rel_fit_ERSP"), "min", 2)
blk_rat   <- format_block(df, c("ratio_fit_MA", "ratio_fit_ERSP"), "max", 2)
blk_ess_b <- format_block(df, c("ESS_MA_bulk", "ESS_ERSP_bulk"), "max", 2)
blk_ess_t <- format_block(df, c("ESS_MA_tail", "ESS_ERSP_tail"), "max", 2)
blk_time  <- format_block(df, c("time_MA", "time_ERSP"), "min", 2)
col_time_rat <- sprintf("%.2f", df$time_ratio)
delta_time <- sprintf("%.2f", df$time_ERSP - df$time_MA)

# New data frame with using the bold letter information
df_print <- data.frame(
  Level = as.character(df$Level),
  blk_abs,
  blk_acc,
  blk_rat,
  blk_ess_b,
  blk_ess_t,
  blk_time,
  delta_time = delta_time,
  Ratio = col_time_rat
)

# Index for package rows
rle_fact <- rle(as.character(df$Factor))
idx <- cumsum(c(1, rle_fact$lengths))

# Separate results in two different data frames:
simres_2_acc <- df_print[,1:9]
simres_2_eff_time <- df_print[,c(1,10:17)]

# Save results to use in the manuscript
simres_2_tab <- saveRDS(object = list(
  df_acc = simres_2_acc,
  df_ESS_time = simres_2_eff_time,
  rle_fact = rle_fact,
  idx = idx), 
  file = "Results/Simulation study/MatchAlign comparison/simres2_marginal_means.rds")

# ─────────────────────────────────────────────────────────────────────────────

