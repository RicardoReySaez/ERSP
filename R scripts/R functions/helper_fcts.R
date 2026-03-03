to_fact_switching <- function(lambdaSamps){
  
  S <- length(lambdaSamps)
  p <- nrow(lambdaSamps[[1]])
  q <- ncol(lambdaSamps[[1]])
  
  lambdaSamps_fact_switch <- matrix(0, nrow = S,
                                    ncol = p*q)
  for(s in 1:S){
    lambdaSamps_fact_switch[s,] <- t(lambdaSamps[[s]]) %>% as.vector()
  }
  
  count <- 1
  names_col <- numeric(p*q)
  for(j in 1:p){
    for(k in 1:q){
      names_col[count] <- paste("LambdaV", j, "_", k, 
                                sep = "")
      count <- count + 1
    }
  }
  
  colnames(lambdaSamps_fact_switch) <- names_col
  
  return(lambdaSamps_fact_switch)
}

to_Erosheva <- function(lambdaSamps){
  
  S <- length(lambdaSamps)
  p <- nrow(lambdaSamps[[1]])
  q <- ncol(lambdaSamps[[1]])
  
  lambdaSamps_fact_switch <- matrix(0, nrow = S,
                                    ncol = p*q)
  for(s in 1:S){
    lambdaSamps_fact_switch[s,] <- lambdaSamps[[s]] %>% as.vector()
  }
  
  count <- 1
  names_col <- numeric(p*q)
  for(j in 1:p){
    for(k in 1:q){
      names_col[count] <- paste("Lam[", j, ",", k,"]",
                                sep = "")
      count <- count + 1
    }
  }
  
  colnames(lambdaSamps_fact_switch) <- names_col
  
  return(lambdaSamps_fact_switch)
}

to_infinitefact <- function(lambdaMat, k, p){
  
  lambdaSamps <- list()
  for(s in 1:nrow(lambdaMat)){
    
    lambdaSamps[[s]] <- matrix(lambdaMat[s,] , 
                               ncol = k, 
                               nrow = p, byrow = T)
    
  }
  return(lambdaSamps)
  
}

jointRot_NOnorm = function(lambda, eta){
  vari = lapply(lambda, varimax, normalize = F)
  loads = lapply(vari, `[[`, 1)
  rots = lapply(vari, `[[`, 2)
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  
  norms = sapply(loads, norm, "2")
  piv = loads[order(norms)][[round(length(lambda)/2)]]
  
  matches = lapply(loads, msfOUT, piv)
  
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  
  return(list(lambda = lamout, eta = etaout))
}

jointRot_pivot = function (lambda, eta, piv = NULL) 
{
  vari = lapply(lambda, varimax, normalize = F)
  loads = lapply(vari, `[[`, 1)
  rots = lapply(vari, `[[`, 2)
  rotfact = mapply(`%*%`, eta, rots, SIMPLIFY = FALSE)
  norms = sapply(loads, norm, "2")
  if(is.null(piv)){
    piv = loads[order(norms)][[round(length(lambda)/2)]]
  }
  matches = lapply(loads, infinitefactor:::msfOUT, piv)
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  etaout = mapply(aplr, rotfact, matches, SIMPLIFY = FALSE)
  return(list(lambda = lamout, eta = etaout))
}

plotmat_tile = function (mat, color = "green", size = 1, title = NULL, args = NULL) 
{
  mat = apply(mat, 2, rev)
  longmat = melt(mat)
  Var1 = Var2 = value = NULL
  p = ggplot(longmat, aes(x = Var2, y = Var1)) + geom_tile(aes(fill = value), 
                                                           colour = "grey20",
                                                           size = size)
  if (color == "green") 
    p = p + scale_fill_gradient2(low = "#3d52bf", high = "#33961b", 
                                 mid = "white")
  if (color == "red") 
    p = p + scale_fill_gradient2(low = "#191970", high = "#800000", 
                                 mid = "white")
  if (color == "wes") 
    p = p + scale_fill_gradient2(low = "#046C9A", high = "#D69C4E", 
                                 mid = "white")
  p = p + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
                panel.grid.major = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), axis.ticks = element_blank(), 
                axis.text = element_blank(), legend.title = element_text(), 
                plot.title = element_text(hjust = 0.5)) + labs(fill = " ")
  if (!is.null(title)) 
    p = p + ggtitle(title)
  if (!is.null(args)) 
    p = p + args
  p
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

jointRot_Lambda = function (lambda, piv = NULL) 
{
  vari = lapply(lambda, varimax)
  loads = lapply(vari, `[[`, 1)
  norms = sapply(loads, norm, "2")
  if(is.null(piv)){
    piv = loads[order(norms)][[round(length(lambda)/2)]]
  }
  matches = lapply(loads, infinitefactor:::msfOUT, piv)
  lamout = mapply(aplr, loads, matches, SIMPLIFY = FALSE)
  return(list(lambda = lamout))
}


linearMGSP_NA <- function (X, nrun, burn, thin = 1, prop = 1, epsilon = 0.001, 
                           kinit = NULL, adapt = TRUE, output = c("covMean", "covSamples", 
                                                                  "factSamples", "sigSamples", "numFactors"), verbose = TRUE, 
                           dump = FALSE, filename = "samps.Rds", buffer = 10000, augment = NULL) 
{
  if (nrun <= burn) 
    stop("nrun must be larger than burn")
  if (!is.matrix(X)) 
    stop("X must be a matrix")
  if (!is.null(augment)) 
    if (!is.expression(augment)) 
      stop("augment must be an expression (see expression())")
  cm = any(output %in% "covMean")
  cs = any(output %in% "covSamples")
  fs = any(output %in% "factSamples")
  ss = any(output %in% "sigSamples")
  nf = any(output %in% "numFactors")
  p = ncol(X)
  n = nrow(X)
  as = 1
  bs = 0.3
  df = 3
  ad1 = 2.1
  bd1 = 1
  ad2 = 3.1
  bd2 = 1
  adf = 1
  bdf = 1
  b0 = 1
  b1 = 5e-04
  if (is.null(kinit)) 
    kinit = floor(log(p) * 3)
  sp = floor((nrun - burn)/thin)
  
  
  num = 0
  k = kinit
  ps = rgamma(p, as, bs)
  Lambda = matrix(1, nrow = p, ncol = k)
  psijh = matrix(rgamma(p * k, df/2, df/2), nrow = p, ncol = k)
  delta = c(rgamma(1, ad1, bd1), rgamma(k - 1, ad2, bd2))
  tauh = cumprod(delta)
  Plam = t(t(psijh) * (tauh))
  if (cm) 
    COVMEAN = matrix(0, nrow = p, ncol = p)
  if (cs) 
    OMEGA = array(dim = c(p, p, sp))
  if (fs) {
    LAMBDA = list()
    ETA = list()
  }
  if (ss) 
    SIGMA = array(dim = c(p, sp))
  if (nf) 
    K = rep(NA, sp)
  ind = 1
  at = ceiling(nrun/100)
  if (verbose) {
    pb = txtProgressBar(style = 3)
  }
  
  
  # --- Update missing data --- #
  X_na <- is.na(X)
  eta = matrix(rnorm(n*k), ncol = k, nrow = n)
  X_pred = eta%*%t(Lambda) + mvtnorm::rmvnorm(n, sigma = diag(1/ps))
  X[X_na] = X_pred[X_na]
  
  VY = apply(X, 2, var)
  scaleMat = sqrt((VY) %*% t(VY))
  X = scale(X)
  
  for (i in 1:nrun) {
    
    # --- Update missing data --- #
    X_pred = eta%*%t(Lambda) + mvtnorm::rmvnorm(n, sigma = diag(1/as.vector(ps)))
    X[X_na] = X_pred[X_na]
    
    eta = eta_lin(Lambda, ps, k, n, X)
    Lambda = lam_lin(eta, Plam, ps, k, p, X)
    psijh = psi_mg(Lambda, tauh, ps, k, p, df)
    delta = del_mg(Lambda, psijh, tauh, delta, k, p, ad1, 
                   bd1, ad2, bd2)
    tauh = cumprod(delta)
    ps = sig_lin(Lambda, eta, k, p, n, X, as, bs)
    Plam = plm_mg(psijh, tauh)
    if (!is.null(augment)) 
      eval(augment)
    if (adapt) {
      prob = 1/exp(b0 + b1 * i)
      uu = runif(1)
      lind = colSums(abs(Lambda) < epsilon)/p
      vec = lind >= prop
      num = sum(vec)
      if (uu < prob) {
        if ((i > 20) & (num == 0) & all(lind < 0.995)) {
          k = k + 1
          Lambda = cbind(Lambda, rep(0, p))
          eta = cbind(eta, rnorm(n))
          psijh = cbind(psijh, rgamma(p, df/2, df/2))
          delta[k] = rgamma(1, ad2, bd2)
          tauh = cumprod(delta)
          Plam = t(t(psijh) * tauh)
        }
        else {
          if (num > 0) {
            k = max(k - num, 1)
            Lambda = Lambda[, !vec, drop = F]
            psijh = psijh[, !vec, drop = F]
            eta = eta[, !vec, drop = F]
            delta = delta[!vec]
            tauh = cumprod(delta)
            Plam = t(t(psijh) * tauh)
          }
        }
      }
    }
    if ((i%%thin == 0) & (i > burn)) {
      if (cm | cs) 
        Omega = (tcrossprod(Lambda) + diag(1/c(ps))) * 
          scaleMat
      if (cm) 
        COVMEAN = COVMEAN + Omega/sp
      if (cs) 
        OMEGA[, , ind] = Omega
      if (fs) {
        LAMBDA[[ind]] = Lambda
        ETA[[ind]] = eta
      }
      if (ss) 
        SIGMA[, ind] = 1/ps
      if (nf) 
        K[ind] = k
      ind = ind + 1
    }
    if (verbose & (i%%at == 0)) 
      setTxtProgressBar(pb, i/nrun)
    if (dump & (i%%buffer == 0)) {
      out = list()
      if (cm) 
        out = c(out, list(covMean = COVMEAN))
      if (cs) 
        out = c(out, list(omegaSamps = OMEGA))
      if (fs) 
        out = c(out, list(lambdaSamps = LAMBDA, etaSamps = ETA))
      if (ss) 
        out = c(out, list(sigmaSamps = SIGMA))
      if (nf) 
        out = c(out, list(numFacts = K))
      saveRDS(out, filename, compress = FALSE)
    }
  }
  if (verbose) 
    close(pb)
  out = list()
  if (cm) 
    out = c(out, list(covMean = COVMEAN))
  if (cs) 
    out = c(out, list(omegaSamps = OMEGA))
  if (fs) 
    out = c(out, list(lambdaSamps = LAMBDA, etaSamps = ETA))
  if (ss) 
    out = c(out, list(sigmaSamps = SIGMA))
  if (nf) 
    out = c(out, list(numFacts = K))
  return(out)
}

