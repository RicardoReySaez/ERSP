# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : RSP_comparison_sim.R                                       ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 28-01-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Two simulation studies to compare RSP implementations and MatchAlign Vs. 
# the efficient implementation of RSP algorithm. 

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(SimDesign)
library(infinitefactor)
library(factor.switching)
library(BayesEFA)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Design conditions
Design <- createDesign(
  N = 500, 
  M = seq(2, 7, 1), 
  J = c(50, 100), 
  Lambda = c("dense", "structured"))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: data-generating function
# ─────────────────────────────────────────────────────────────────────────────

# N <- 500; M <- 2; J <- 50; Lambda <- "dense"
# condition <- list(N=N,M=M,J=J,Lambda=Lambda)
Generate <- function(condition, fixed_objects) {
  Attach(condition)
  # Latent parameters
  if(Lambda == "dense") {
    L <- matrix(rnorm(J*M), nrow = J, ncol = M)
  } else {
    L <- matrix(0, nrow = J, ncol = M)
    groups <- split(1:J, cut(1:J, M, labels = FALSE))
    for (m in 1:M) {
      L[groups[[m]], m] <- rnorm(length(groups[[m]]))
    }
  }
  P <- diag(1/rgamma(n = J, shape = .5, rate = .5))
  S <- tcrossprod(L) + P
  
  # Simulate data
  dat <- mvnfast::rmvn(n = N, mu = rep(0, J), sigma = S)
  dat
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: data-analysis function
# ─────────────────────────────────────────────────────────────────────────────

Analyse <- function(condition, dat, fixed_objects) {
  # Estimate bayesian DL model
  fit <- infinitefactor::linearDL(X = dat, nrun = 11000, burn = 1000, 
                                  thin = 10, k = condition$M, 
                                  verbose = FALSE)
  
  # Put draws in row-major order (for factor.switching)
  posterior_lambda <- matrix(NA, ncol = ncol(dat) * condition$M, nrow = 1000)
  colnames(posterior_lambda) <- paste0("LambdaV", rep(1:ncol(dat), each = condition$M), 
                                       "_", 1:condition$M)
  # Row-major format
  for(s in 1:1000) {
    posterior_lambda[s,] <- c(t(fit$lambdaSamps[[s]]))
  }
  
  # Original Rotation-Sign-Permutation algorithm (Papastamoulis & Ntzoufras, 2022)
  time_RSP <- system.time({
    orig_RSP <- factor.switching::rsp_exact(
      lambda_mcmc = posterior_lambda)
  })
  
  # Efficient Rotation-Sign-Permutation algorithm (Rey-Sáez y Revuelta, 2026)
  time_ERSP <- system.time({
    Eff_RSP <- BayesEFA::rsp_align(
      lambda_draws = posterior_lambda, 
      n_items = ncol(dat), 
      n_factors = condition$M, 
      n_chains = 1, 
      format = "row_major")
  })
  
  # Model-implied explained-coraviance matrix
  lambda_postmean <- infinitefactor::lmean(
    lapply(fit$lambdaSamps, tcrossprod)
    )
  
  # Aligned-implied explained-coraviance matrices
  lmean_ERSP <- tcrossprod(Eff_RSP$Lambda_star)
  lmean_RSP  <- tcrossprod(orig_RSP$lambda_hat)
  
  # Absolute Accuracy Metric
  abs_fit_RSP  <- norm(lmean_RSP  - lambda_postmean, type = "F")
  abs_fit_ERSP <- norm(lmean_ERSP  - lambda_postmean, type = "F")
  
  # Relative Accuracy Metric
  rel_fit_RSP  <- abs_fit_RSP / norm(lambda_postmean, type = "F")
  rel_fit_ERSP <- abs_fit_ERSP / norm(lambda_postmean, type = "F")
  
  # Ratio of times
  time_ratio <- time_RSP["elapsed"] / time_ERSP["elapsed"]
  
  # Return results
  res <- nc(
    abs_fit_RSP, 
    abs_fit_ERSP, 
    rel_fit_RSP, 
    rel_fit_ERSP,
    time_RSP = time_RSP["elapsed"], 
    time_ERSP = time_ERSP["elapsed"],
    time_ratio)
  
  # Return results
  return(res)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: summarise function
# ─────────────────────────────────────────────────────────────────────────────

Summarise <- function(condition, results, fixed_objects) colMeans(results)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: simulation
# ─────────────────────────────────────────────────────────────────────────────

# Global seeds per condition
# (inside Generate we control specific seed per replica)
global_seeds <- genSeeds(design = Design, iseed = 2026)

# Run simulation study
res <- runSimulation(
  design         = Design,
  replications   = 20,
  generate       = Generate, 
  analyse        = Analyse,
  summarise      = Summarise,
  parallel       = TRUE,
  store_results  = TRUE,
  save           = TRUE,
  max_errors     = 1000L,
  ncores         = 10, 
  save_details   = list(tmpfilename = "tmp_rsp_comparison_simres.rds",
                        out_rootdir = "Results/Simulation study/RSP comparison"), 
  control        = list(store_Random.seeds = TRUE),
  packages       = c("mvnfast", "infinitefactor", "factor.switching", "BayesEFA"),
  filename       = "rsp_comparison_simresults.rds"
)

# Store again results as rdata
save(res, file = "Results/Simulation study/RSP comparison/rsp_comparison_simresults.rdata")

# ─────────────────────────────────────────────────────────────────────────────
