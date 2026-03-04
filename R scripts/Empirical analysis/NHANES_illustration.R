# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : NHANES_illustration.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 03-02-2026                                                 ║
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
library(BayesianEFA)
library(tidyverse)
library(coda)
library(posterior)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Load and fit both MSGP Bayesian factor models
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
    )[,-1]
  )
}

# Load NHANES data
df <- readRDS("Data/df_chem.rds")

# Fit the same Multiplicative Shrinkage Gamma Priors Model with 25 dimensions
set.seed(1) # Same as poworoznek et al.
out_25 <- linearMGSP_NA(log(df) |> as.matrix() |> scale(), verbose = T,
                        nrun = 10000, burn = 5000, thin = 5, 
                        kinit = 25, adapt = F, output = "factSamples")

# Fit the same Multiplicative Shrinkage Gamma Priors Model with 50 dimensions
out_50 <- linearMGSP_NA(log(df) |> as.matrix() |> scale(), verbose = T,
                        nrun = 10000, burn = 5000, thin = 5, 
                        kinit = 50, adapt = F, output = "factSamples")

# Save both models to avoid refitting
NHANES_fits <- list(out_25 = out_25, out_50 = out_50)
saveRDS(NHANES_fits, file = "Results/Empirical analysis/NHANES_fits.rds")

# If you want, load from here to avoid model refitting
NHANES_fits <- readRDS("Results/Empirical analysis/NHANES_fits.rds")
out_25 <- NHANES_fits$out_25
out_50 <- NHANES_fits$out_50

# Extract posterior draws
lambdaSamps_25 <- out_25$lambdaSamps
etaSamps_25    <- out_25$etaSamps
lambdaSamps_50 <- out_50$lambdaSamps
etaSamps_50    <- out_50$etaSamps

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: MatchAlign estimation and metrics
# ─────────────────────────────────────────────────────────────────────────────

# MatchAlign rotation for M = 25
time_MatchAlign_25 <- system.time({
  out_MatchAlign_25<- jointRot_NOnorm(lambdaSamps_25, etaSamps_25)
})

# MatchAlign rotation for M = 50
time_MatchAlign_50 <- system.time({
  out_MatchAlign_50<- jointRot_NOnorm(lambdaSamps_50, etaSamps_50)
})

# Save aligned posterior draws
Lambda_MatchAlign_25 <- out_MatchAlign_25$lambda
Lambda_MatchAlign_50 <- out_MatchAlign_50$lambda

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Efficient-RSP estimation and metrics
# ─────────────────────────────────────────────────────────────────────────────

# Prepare draws for Efficient RSP
lambdaSamps_ERSP_25 <- to_fact_switching(lambdaSamps_25)
lambdaSamps_ERSP_50 <- to_fact_switching(lambdaSamps_50)

# Efficient-RSP for M = 25
time_ERSP_25 <- system.time({
  out_ERSP_25 <- BayesianEFA::rsp_align(
    lambda_draws = lambdaSamps_ERSP_25, 
    n_items = ncol(df), 
    n_factors = 25, 
    n_chains = 1,
    format = "row_major")
})

# Efficient-RSP for M = 50
time_ERSP_50 <- system.time({
  out_ERSP_50 <- BayesianEFA::rsp_align(
    lambda_draws = lambdaSamps_ERSP_50, 
    n_items = ncol(df), 
    n_factors = 50, 
    n_chains = 1,
    format = "row_major")
})

# Save aligned posterior draws in a usefull format
Lambda_ERSP_25 <- to_infinitefact(
  lambdaMat = out_ERSP_25$Lambda_hat_mcmc, 
  k = 25, 
  p = ncol(df))
Lambda_ERSP_50 <- to_infinitefact(
  lambdaMat = out_ERSP_50$Lambda_hat_mcmc, 
  k = 50, 
  p = ncol(df))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Accuracy metric, ESS and estimation time
# ─────────────────────────────────────────────────────────────────────────────

# Identified explained covariance matrix posterior mean
LLT_25 <- lmean(lapply(lambdaSamps_25, tcrossprod))
LLT_50 <- lmean(lapply(lambdaSamps_50, tcrossprod))

# Accuracy metrics
norm_base_25 <- norm(LLT_25, type = "F")
norm_base_50 <- norm(LLT_50, type = "F")

# Retuls for M = 25
MatchAlign_cov_25  <- tcrossprod(lmean(Lambda_MatchAlign_25))
MatchAlign_norm_25 <- norm(MatchAlign_cov_25 - LLT_25, type = "F")
MatchAlign_ESS_25  <- stan_ESS_est(Lambda_MatchAlign_25)
ERSP_cov_25  <- tcrossprod(lmean(Lambda_ERSP_25))
ERSP_norm_25 <- norm(ERSP_cov_25 - LLT_25, type = "F")
ERSP_ESS_25  <- stan_ESS_est(Lambda_ERSP_25)

# Retuls for M = 50
MatchAlign_cov_50  <- tcrossprod(lmean(Lambda_MatchAlign_50))
MatchAlign_norm_50 <- norm(MatchAlign_cov_50 - LLT_50, type = "F")
MatchAlign_ESS_50  <- stan_ESS_est(Lambda_MatchAlign_50)
ERSP_cov_50  <- tcrossprod(lmean(Lambda_ERSP_50))
ERSP_norm_50 <- norm(ERSP_cov_50 - LLT_50, type = "F")
ERSP_ESS_50  <- stan_ESS_est(Lambda_ERSP_50)


# Create a results data frame
NHANES_results <- data.frame(
  Method = c("MatchAlign", "ERSP"),
  # M = 25
  Acc_25 = c(MatchAlign_norm_25, ERSP_norm_25),
  rel_Acc_25 = 100 * c(MatchAlign_norm_25, ERSP_norm_25)/norm_base_25,
  ESS_bulk_25 = 100 * c(MatchAlign_ESS_25[1], ERSP_ESS_25[1]) / length(lambdaSamps_25),
  ESS_tail_25 = 100 * c(MatchAlign_ESS_25[2], ERSP_ESS_25[2]) / length(lambdaSamps_25),
  Time_25 = c(time_MatchAlign_25["elapsed"], time_ERSP_25["elapsed"]),
  # M = 50
  Acc_50 = c(MatchAlign_norm_50, ERSP_norm_50),
  rel_Acc_50 = 100 * c(MatchAlign_norm_50, ERSP_norm_50)/norm_base_50,
  ESS_bulk_50 = 100 * c(MatchAlign_ESS_50[1], ERSP_ESS_50[1]) / length(lambdaSamps_50),
  ESS_tail_50 = 100 * c(MatchAlign_ESS_50[2], ERSP_ESS_50[2]) / length(lambdaSamps_50),
  Time_50 = c(time_MatchAlign_50["elapsed"], time_ERSP_50["elapsed"])
)

# Save NHANES results
saveRDS(NHANES_results, file = "Results/Empirical analysis/NHANES_results_table.rds")

# ─────────────────────────────────────────────────────────────────────────────
