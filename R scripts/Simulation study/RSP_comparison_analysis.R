# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : RSP_comparison_analysis.R                                  ║
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
# exact-RSP and the efficient-RSP versions
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
# SECTION 3: Agreement between methods
# ─────────────────────────────────────────────────────────────────────────────

# Load simulation results
simres_s1 <- readRDS("Results/Simulation study/RSP comparison/rsp_comparison_simresults.rds")
all_res <- SimExtract(simres_s1, what = "results")

# Correlation between variables
cor(all_res$abs_fit_RSP, all_res$abs_fit_ERSP)

# Bivariate scatterplot
biv_plot <- ggplot(all_res, aes(x = abs_fit_RSP, y = abs_fit_ERSP)) +
  geom_abline(
    slope = 1, 
    intercept = 0, 
    linetype = "dashed", 
    color = "gray50",
    size = 0.6) +
  geom_point(
    shape = 16,          
    fill =  "#2c3e50",  
    size = 2.5,         
    alpha = 0.5         
  ) +
  coord_equal() +
  labs(
    x = "Exact RSP Absolute Error",
    y = "Efficient RSP Absolute Error"
  ) +
  theme_classic(base_size = 12, base_family = "serif") + 
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"), 
    axis.line = element_line(linewidth = 0.8),
    plot.margin = margin(20, 20, 20, 20)
  )

# Save bivariate scatterplot
png(filename = "Figures/RSP_agreement.png", width = 2000, height = 2000, res = 300)
biv_plot
dev.off()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Differences in time-spent aligning draws
# ─────────────────────────────────────────────────────────────────────────────

# Prepare wide data frame
tab_wide <- simres_s1 %>%
  mutate(Lambda = recode(
    Lambda, 
    "dense"      = "Independent", 
    "structured" = "Sparse"
  )) %>%
  filter(J %in% c(50, 100)) %>%
  group_by(M, Lambda, J) %>%
  summarise(
    RSP_min  = mean(time_RSP,  na.rm = TRUE) / 60,
    ERSP_min = mean(time_ERSP, na.rm = TRUE) / 60,
    ratio    = mean(time_RSP / time_ERSP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = J,
    values_from = c(RSP_min, ERSP_min, ratio),
    names_glue  = "{.value}_J{J}"
  ) %>%
  arrange(Lambda, M) %>% 
  select(
    Lambda, M,
    RSP_min_J50,  ERSP_min_J50,  ratio_J50,
    RSP_min_J100, ERSP_min_J100, ratio_J100
  )

# Final LaTeX table
kb <- tab_wide %>%
  select(-Lambda) %>% 
  kable(
    format = "latex",
    booktabs = TRUE,
    linesep = "",   
    escape = FALSE,
    col.names = c(
      "$M$",
      "RSP", "ERSP", "Ratio",
      "RSP", "ERSP", "Ratio"
    ),
    align = c("c","r","r","r","r","r","r"),
    caption = "Computational performance comparison: Execution time (minutes) and speedup ratios comparing Exact and Efficient RSP methods."
  ) %>%
  add_header_above(c(" " = 1, "J = 50" = 3, "J = 100" = 3)) %>%
    pack_rows(
    index = table(as.character(tab_wide$Lambda)), 
    latex_gap_space = "1.5em" 
  ) %>%
  column_spec(4, latex_column_spec = "r@{\\hspace{1cm}}") %>%
  kable_styling(latex_options = c("hold_position"))

# Save the data frame in wide format
saveRDS(object = tab_wide, file = "Results/Simulation study/RSP comparison/simres1_execution_time.rds")

# ─────────────────────────────────────────────────────────────────────────────

