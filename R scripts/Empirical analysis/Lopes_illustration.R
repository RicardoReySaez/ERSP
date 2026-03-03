# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : Lopes_illustration.R                                       ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 04-02-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────
# Detailed description of the script's purpose and functionality
#
#
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(infinitefactor)
library(BayesEFA)
library(tidyverse)
library(posterior)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Load and fit the MSGP Bayesian factor models
# ─────────────────────────────────────────────────────────────────────────────

# Load Poworoznek et al. utils
# Retrieved from here: https://github.com/fedfer/Section_5_MatchAlign/tree/main
source("R scripts/R functions/helper_fcts.R")

# Function to compute ESS metrics using posterior
stan_ESS_est <- function(Lambda_list) {
  colMeans(
    posterior::summarise_draws(
      posterior::as_draws_matrix(
        do.call(
          rbind,
          lapply(
            Lambda_list, function(x) c(x)
          )
        )
      ),
      posterior::ess_bulk,
      posterior::ess_tail
    )[, -1]
  )
}

# Load Lopes y West (2004) dataset (it's standardized by-default)
data <- ier
colnames(data) <- paste0("V", 1:6)

# Fit the same Multiplicative Shrinkage Gamma Priors Model
set.seed(1)
out_LW <- linearMGSP(data,
  verbose = T, nrun = 10000, burn = 5000, thin = 5,
  kinit = ncol(data) - 1, adapt = F, output = "factSamples"
)

# Save posterior draws
Lambda <- out_LW$lambdaSamps
Eta <- out_LW$etaSamps

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Apply MatchAlign, RSP and ERSP alignments
# ─────────────────────────────────────────────────────────────────────────────

# MatchAlign rotation
time_MatchAlign <- system.time({
  out_MatchAlign <- jointRot_NOnorm(Lambda, Eta)
})

# RSP rotation
time_RSP <- system.time({
  out_RSP <- factor.switching::rsp_exact(to_fact_switching(Lambda))
})

# ERSP rotation
time_ERSP <- system.time({
  out_ERSP <- BayesEFA::rsp_align(
    lambda_draws = to_fact_switching(Lambda),
    n_items = ncol(data),
    n_factors = 5, 
    n_chains = 1,
    format = "row_major"
  )
})

# Change RSP and ERSP draws to infinitefactor format
lambda_RSP <- to_infinitefact(
  lambdaMat = out_RSP$lambda_reordered_mcmc,
  k = 5,
  p = ncol(data)
)
lambda_ERSP <- to_infinitefact(
  lambdaMat = out_ERSP$lambda_reordered_mcmc,
  k = 5,
  p = ncol(data)
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Accuracy metric, ESS and estimation time
# ─────────────────────────────────────────────────────────────────────────────

# Identified explained covariance matrix posterior mean
LLT <- lmean(lapply(Lambda, tcrossprod))
norm_LLT <- norm(LLT, type = "F")

# MatchAlign metrics
MatchAlign_cov <- tcrossprod(lmean(out_MatchAlign$lambda))
MatchAlign_norm <- norm(MatchAlign_cov - LLT, type = "F")
MatchAlign_ESS <- stan_ESS_est(out_MatchAlign$lambda)

# RSP metrics
RSP_cov <- tcrossprod(lmean(lambda_RSP))
RSP_norm <- norm(RSP_cov - LLT, type = "F")
RSP_ESS <- stan_ESS_est(lambda_RSP)

# ERSP metrics
ERSP_cov <- tcrossprod(lmean(lambda_ERSP))
ERSP_norm <- norm(ERSP_cov - LLT, type = "F")
ERSP_ESS <- stan_ESS_est(lambda_ERSP)

# Create a results data frame
LW_results <- data.frame(
  Method = c("MatchAlign", "RSP", "ERSP"),
  Acc = c(MatchAlign_norm, RSP_norm, ERSP_norm),
  rel_Acc = 100 * c(MatchAlign_norm, RSP_norm, ERSP_norm) / norm_LLT,
  ESS_bulk = 100 * c(MatchAlign_ESS[1], RSP_ESS[1], ERSP_ESS[1]) / length(Lambda),
  ESS_tail = 100 * c(MatchAlign_ESS[2], RSP_ESS[2], ERSP_ESS[2]) / length(Lambda),
  Time = c(time_MatchAlign["elapsed"], time_RSP["elapsed"], time_ERSP["elapsed"])
)

# Save NHANES results
saveRDS(LW_results, file = "Results/Empirical analysis/LW_results_table.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Aligned factor loading plot
# ─────────────────────────────────────────────────────────────────────────────

# Auxiliar function to prepare data frames for the plot
prep_df <- function(mat, name) {
  as.data.frame(as.matrix(mat)) %>%
    mutate(Item = rownames(.)) %>%
    pivot_longer(-Item, names_to = "Factor", values_to = name)
}

# Generate posterior means and 99% CIs matrices
df_mean <- prep_df(apply(simplify2array(lambda_ERSP), 1:2, mean), "Mean")
df_lower <- prep_df(apply(simplify2array(lambda_ERSP), 1:2, quantile, probs = .005), "Lower")
df_upper <- prep_df(apply(simplify2array(lambda_ERSP), 1:2, quantile, probs = .995), "Upper")

# Generate the plot data frame
df_final <- df_mean %>%
  left_join(df_lower, by = c("Item", "Factor")) %>%
  left_join(df_upper, by = c("Item", "Factor")) %>%
  mutate(
    Factor = gsub("V", "F", Factor),
    Item = paste0("Item ", Item),
    zero_out_99 = ifelse(Lower > 0 | Upper < 0, "bold", "plain"),
    Label = paste0(
      sprintf("%.2f", Mean), "\n(",
      sprintf("%.2f", Lower), ", ", sprintf("%.2f", Upper), ")"
    )
  )

# Order items based on their loadings
item_order <- df_final %>%
  group_by(Item) %>%
  summarise(
    Max_load = max(abs(Mean)),
    Main_factor = Factor[which.max(abs(Mean))]
  ) %>%
  arrange(Main_factor, desc(Max_load)) %>%
  pull(Item)
df_final$Item <- factor(df_final$Item, levels = rev(item_order))

# Now, rename items using their original names and reorder factors
df_final <- df_final |>
  mutate(
    Item = recode(
      Item,
      "Item 1" = "USD",
      "Item 2" = "CAD",
      "Item 3" = "JPY",
      "Item 4" = "FRF",
      "Item 5" = "ITL",
      "Item 6" = "DEM"
    ),
    Item = factor(Item, levels = c("USD", "CAD", "JPY", "FRF", "ITL", "DEM"))
  )

# Final plot
loads_plot <- ggplot(df_final, aes(x = Factor, y = Item, fill = Mean)) +
  geom_tile(color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#ef8a62", mid = "white", high = "#67a9cf",
    midpoint = 0, limits = c(-1, 1),
    name = "Loading"
  ) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = Label, fontface = zero_out_99),
    color = "black",
    size = 4.5,
    family = "serif",
    lineheight = 0.9
  ) +
  theme_bw(base_size = 16, base_family = "serif") +
  labs(
    y = NULL,
    x = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 14),
  ) +
  coord_fixed()

# Save plot
ggsave(
  filename = "Figures/factor_loading_matrix.png",
  plot = loads_plot,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ─────────────────────────────────────────────────────────────────────────────
