# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : MatchAlign_comparison.R                                    ║
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
library(mvnfast)
library(infinitefactor)
library(BayesEFA)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# MatchAlign rotation only for lambda (added normalize = FALSE)
# Retrieved from: https://github.com/fedfer/Section_5_MatchAlign/blob/main/helper_fcts.R
jointRot_Lambda = function (lambda, piv = NULL) {
  vari = lapply(lambda, varimax, normalize = FALSE)
  loads = lapply(vari, `[[`, 1)
  norms = sapply(loads, norm, "2")
  if(is.null(piv)){
    piv = loads[order(norms)][[round(length(lambda)/2)]]
  }
  matches = lapply(loads, infinitefactor:::msfOUT, piv)
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  return(list(lambda = lamout))
}

# Design conditions
Design <- createDesign(
  N = c(100, 250, 500), 
  M = seq(5, 30, 5), 
  J = c(50, 100, 200), 
  Lambda = c("dense", "structured"))

# Pass JointRot_Lambda as fixed object
fixed_objects <- list(jointRot_Lambda = jointRot_Lambda)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: data-generating function
# ─────────────────────────────────────────────────────────────────────────────

# N <- 100; M <- 20; J <- 200; Lambda <- "dense"
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
      L[groups[[m]], m] <- rnorm(length(groups[[m]]), mean = 0, sd = 1)
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
  fit <- infinitefactor::linearDL(
    X = dat, nrun = 11000, burn = 1000, 
    thin = 10, k = condition$M, 
    verbose = TRUE)
  
  # Put draws in column-major order (for BayesEFA)
  posterior_lambda <- matrix(NA, ncol = ncol(dat) * condition$M, nrow = 1000)
  
  # column-major format
  for(s in 1:1000) {
    posterior_lambda[s,] <- c(fit$lambdaSamps[[s]])
  }
  
  # MatchAlign Algorithm (Poworoznek et al., 2025)
  time_MA <- system.time({
    MatchAlign <- fixed_objects$jointRot_Lambda(lambda = fit$lambdaSamps)
  })
  
  # Efficient Rotation-Sign-Permutation algorithm (Rey-Sáez y Revuelta, 2026)
  time_ERSP <- system.time({
    Eff_RSP <- BayesEFA::rsp_align(
      lambda_draws = posterior_lambda, 
      n_items = ncol(dat), 
      n_factors = condition$M, 
      n_chains = 1,
      format = "column_major")
  })
  
  # Model-implied explained-coraviance matrix
  lambda_postmean <- infinitefactor::lmean(
    lapply(fit$lambdaSamps, tcrossprod)
  )

  # Aligned-implied explained-coraviance matrices
  lmean_base <- tcrossprod(infinitefactor::lmean(fit$lambdaSamps))
  lmean_MA   <- tcrossprod(infinitefactor::lmean(MatchAlign$lambda))
  lmean_ERSP <- tcrossprod(Eff_RSP$Lambda_star)
  
  # Absolute Accuracy Metric
  abs_fit_base <- norm(lmean_base - lambda_postmean, type = "F")
  abs_fit_MA   <- norm(lmean_MA   - lambda_postmean, type = "F")
  abs_fit_ERSP <- norm(lmean_ERSP - lambda_postmean, type = "F")
  
  # Relative Accuracy Metric
  rel_fit_base <- abs_fit_base / norm(lambda_postmean, type = "F")
  rel_fit_MA   <- abs_fit_MA   / norm(lambda_postmean, type = "F")
  rel_fit_ERSP <- abs_fit_ERSP / norm(lambda_postmean, type = "F")
  
  # Relative accuracy w.r.t. unaligned solution
  ratio_fit_MA   <- abs_fit_base / abs_fit_MA
  ratio_fit_ERSP <- abs_fit_base / abs_fit_ERSP
  
  # Ratio of times
  time_ratio <- time_MA["elapsed"] / time_ERSP["elapsed"]
  
  # Effective Sample Size: MatchAlign
  post_MA <- posterior::summarise_draws(
    posterior::as_draws_matrix(
      do.call(rbind, lapply(MatchAlign$lambda, function(x) c(x)))
    ), 
    posterior::default_convergence_measures()
  )
  
  # Effective Sample Size: Efficient RSP
  post_ERSP <- posterior::summarise_draws(
    posterior::as_draws_matrix(Eff_RSP$Lambda_hat_mcmc),
    posterior::default_convergence_measures())
  
  # Percentage of the total number of samples
  ESS_MA_bulk <- colMeans(post_MA[,3])/10
  ESS_MA_tail <- colMeans(post_MA[,4])/10
  ESS_ERSP_bulk <- colMeans(post_ERSP[,3])/10
  ESS_ERSP_tail <- colMeans(post_ERSP[,4])/10
  
  # Return results
  res <- nc(
    abs_fit_base,
    abs_fit_MA, 
    abs_fit_ERSP, 
    rel_fit_base,
    rel_fit_MA, 
    rel_fit_ERSP,
    ratio_fit_MA,
    ratio_fit_ERSP,
    ESS_MA_bulk, 
    ESS_MA_tail,
    ESS_ERSP_bulk, 
    ESS_ERSP_tail,
    time_MA = time_MA["elapsed"], 
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
  replications   = 100,
  generate       = Generate, 
  analyse        = Analyse,
  summarise      = Summarise,
  fixed_objects  = fixed_objects,
  parallel       = TRUE,
  store_results  = TRUE,
  save           = TRUE,
  max_errors     = 1000L,
  ncores         = 10, 
  save_details   = list(tmpfilename = "tmp_MatchAlign_comp_simres.rds",
                        out_rootdir = "Results/Simulation study/MatchAlign comparison"), 
  control        = list(store_Random.seeds = TRUE),
  packages       = c("infinitefactor", "mvnfast", "BayesEFA", "posterior"),
  filename       = "MatchAlign_comp_simres.rds"
)

# Store again results as rdata
save(res, file = "Results/Simulation study/MatchAlign comparison/MatchAlign_comp_simres.rdata")

# ─────────────────────────────────────────────────────────────────────────────
