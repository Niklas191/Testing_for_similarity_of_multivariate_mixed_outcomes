# === This file is related to =================================================================================================================================
# Testing for similarity of multivariate mixed outcomes using generalised joint regression models with application to efficacy-toxicity responses
# Niklas Hagemann, Giampiero Marra, Frank Bretz and Kathrin Moellenhoff

# === File description ========================================================================================================================================
# This file contains reimplementations of functions from the R package BinNor excluding all print and cat commands.

# === Functions ===============================================================================================================================================

validation.bin_quite <- function(no.bin, prop.vec.bin = NULL){
  if ((no.bin < 0) | (floor(no.bin) != no.bin)) {
    stop("Number of binary variables \nmust be a non-negative integer\n")
  }
  else if (!is.null(prop.vec.bin)) {
    if (no.bin == 0) {
      stop("Proportion vector is specified while no.bin=0")
    }
    else if ((min(prop.vec.bin) <= 0) | (max(prop.vec.bin) >= 
                                         1)) {
      stop("Proportions for binary variables must be between 0 and 1!\n")
    }
    else if (length(prop.vec.bin) != no.bin) {
      stop("Proportion vector is misspecified, dimension is wrong!\n")
    }
  }
  else if (is.null(prop.vec.bin)) {
    if (no.bin > 0) {
      stop("Proportion vector is not specified while no.bin > 0")
    }
  }
}

validation.nor_quiet <- function(no.nor, mean.vec.nor = NULL, var.nor = NULL){
  if ((no.nor < 0) | (floor(no.nor) != no.nor)) {
    stop("Number of normal variables \nmust be an integer whose value or 0 !\n")
  }
  if (!is.null(mean.vec.nor) & no.nor == 0) {
    stop("Mean vector for the normal part is specified while no.nor=0!\n")
  }
  if (!is.null(var.nor) & no.nor == 0) {
    stop("Vector of variances for the normal part is specified while no.nor=0!\n")
  }
  if (is.null(mean.vec.nor) & no.nor > 0) {
    stop("Mean vector for the normal part is not specified while no.nor>0!\n")
  }
  if (is.null(var.nor) & no.nor > 0) {
    stop("Vector of variances for the normal part is not specified while no.nor>0!\n")
  }
  if (!is.null(mean.vec.nor) & !is.null(var.nor) & no.nor > 
      0) {
    if (length(mean.vec.nor) != no.nor) {
      stop("Mean vector for the normal part is misspecified, \ndimension is wrong!\n")
    }
    if (length(var.nor) != no.nor) {
      stop("Vector of variances for the normal part is misspecified, \ndimension is wrong!\n")
    }
    if (min(var.nor) <= 0) {
      stop("Variances must be positive!\n")
    }
  }
}

validation.range_quiet <- function(no.bin, no.nor, prop.vec.bin = NULL, corr.mat){
  d = no.bin + no.nor
  sigma = corr.mat
  p = prop.vec.bin
  q = 1 - p
  L_sigma = diag(d)
  U_sigma = diag(d)
  if (no.bin > 0) {
    for (i in 1:no.bin) {
      for (j in 1:no.bin) {
        if (i != j) 
          L_sigma[i, j] = L_sigma[j, i] = max(-sqrt((p[i] * 
                                                       p[j])/(q[i] * q[j])), -sqrt((q[i] * q[j])/(p[i] * 
                                                                                                    p[j])))
        if (i != j) 
          U_sigma[i, j] = U_sigma[j, i] = min(sqrt((p[i] * 
                                                      q[j])/(q[i] * p[j])), sqrt((q[i] * p[j])/(p[i] * 
                                                                                                  q[j])))
      }
    }
  }
  if (no.bin > 0 & no.nor > 0) {
    for (i in (no.bin + 1):d) {
      for (j in 1:no.bin) {
        L_sigma[i, j] = L_sigma[j, i] = -dnorm(qnorm(p[j]))/sqrt(p[j] * 
                                                                   q[j])
        U_sigma[i, j] = U_sigma[j, i] = dnorm(qnorm(p[j]))/sqrt(p[j] * 
                                                                  q[j])
      }
    }
  }
  if (no.nor > 0) {
    for (i in (no.bin + 1):d) {
      for (j in (no.bin + 1):d) {
        L_sigma[i, j] = L_sigma[j, i] = -1
        U_sigma[i, j] = U_sigma[j, i] = 1
      }
    }
  }
  valid.state = TRUE
  for (i in 1:d) {
    for (j in 1:d) {
      if (j >= i) {
        if (sigma[i, j] < L_sigma[i, j] | sigma[i, j] > 
            U_sigma[i, j]) {
          valid.state = FALSE
        }
      }
    }
  }
  if (valid.state == FALSE) 
    stop("All correlations must be in feasible range!")
}

validation.corr_quiet <- function(no.bin, no.nor, prop.vec.bin = NULL, corr.vec = NULL, corr.mat = NULL){
  d = no.bin + no.nor
  validation.bin_quite(no.bin, prop.vec.bin)
  if (is.null(corr.vec) & is.null(corr.mat)) {
    stop("You must specify full correlation matrix OR vector of elements below the diagonal")
  }
  if (!is.null(corr.vec) & !is.null(corr.mat)) {
    corr.mat.from.corr.vec = lower.tri.to.corr.mat(corr.vec, 
                                                   d)
    if (sum(dim(corr.mat.from.corr.vec) == dim(corr.mat)) != 
        2) {
      stop("corr.vec and corr.mat are non-conformable")
    }
    if (sum(corr.mat.from.corr.vec == corr.mat) != (d * d)) {
      stop("Correlation matrix from corr.vec and corr.mat are not the same")
    }
  }
  if (!is.null(corr.vec)) {
    if (length(corr.vec) != (d * (d - 1)/2)) {
      stop("Vector of correlations is misspecified, dimension is wrong!\n")
    }
    if ((min(corr.vec) <= -1) | (max(corr.vec) >= 1)) {
      stop("Correlations must be between -1 and 1!\n")
    }
    corr.mat.from.corr.vec = lower.tri.to.corr.mat(corr.vec, 
                                                   d)
    if (is.positive.definite(corr.mat.from.corr.vec) == FALSE) {
      stop("Specified correlation matrix (from corr.vec) is not positive definite! \n")
    }
    validation.range_quiet(no.bin, no.nor, prop.vec.bin, corr.mat.from.corr.vec)
  }
  if (!is.null(corr.mat)) {
    if (dim(corr.mat)[1] != d | dim(corr.mat)[2] != d) {
      stop("Correlation matrix dimension is wrong!\n")
    }
    if (is.positive.definite(corr.mat) == FALSE) {
      stop("Specified correlation matrix is not positive definite! \n")
    }
    if (isSymmetric(corr.mat) == FALSE) {
      stop("Specified correlation matrix is not symmetric! \n")
    }
    validation.range_quiet(no.bin, no.nor, prop.vec.bin, corr.mat)
  }
}

compute.sigma.star_quiet <- function(no.bin, no.nor, prop.vec.bin = NULL, corr.vec = NULL, corr.mat = NULL){
  d = no.bin + no.nor
  validation.corr_quiet(no.bin, no.nor, prop.vec.bin, corr.vec = corr.vec, corr.mat = corr.mat)
  if (is.null(corr.mat)) {
    corr.mat = lower.tri.to.corr.mat(corr.vec, d)
  }
  sigma = corr.mat
  p = prop.vec.bin
  q = 1 - p
  if (no.bin != 0) {
    sigmaBB = diag(no.bin)
    for (i in 1:no.bin) {
      for (j in 1:no.bin) {
        if (i != j) 
          sigmaBB[i, j] = phi2tetra(sigma[i, j], c(p[i], 
                                                   p[j]))
      }
    }
  }
  if (no.bin > 0 & no.nor > 0) {
    sigmaBN = sigma
    for (i in (no.bin + 1):d) {
      for (j in 1:no.bin) {
        sigmaBN[i, j] = sigmaBN[i, j]/(dnorm(qnorm(p[j]))/sqrt(p[j] * 
                                                                 q[j]))
      }
    }
    sigmaBN = sigmaBN[(no.bin + 1):d, 1:no.bin]
    sigma_star = sigma
    sigma_star[1:no.bin, 1:no.bin] = sigmaBB
    sigma_star[(no.bin + 1):d, 1:no.bin] = sigmaBN
    sigma_star[1:no.bin, (no.bin + 1):d] = t(sigmaBN)
  }
  if (no.bin > 0 & no.nor == 0) {
    sigma_star = sigmaBB
  }
  if (no.bin == 0 & no.nor > 0) {
    sigma_star = sigma
  }
  PD = TRUE
  temp = NULL
  eigenv = eigen(sigma_star)$value
  if (is.positive.definite(sigma_star) == FALSE) {
    temp = sigma_star
    sigma_star = as.matrix(nearPD(sigma_star, corr = TRUE, 
                                  keepDiag = TRUE)$mat)
    sigma_star = (sigma_star + t(sigma_star))/2
    PD = FALSE
  }
  return(list(sigma_star = sigma_star, nonPD = temp, PD = PD, 
              eigenv = eigenv))
}

jointly.generate.binary.normal_quiet <- function(no.rows, no.bin, no.nor, prop.vec.bin = NULL, mean.vec.nor = NULL, var.nor = NULL, sigma_star = NULL, corr.vec = NULL, corr.mat = NULL, continue.with.warning = TRUE) {
  if ((no.rows < 1) | (floor(no.rows) != no.rows)) {
    stop("Number of rows must be an integer whose value is at least 1!\n")
  }
  d = no.bin + no.nor
  validation.bin_quite(no.bin, prop.vec.bin)
  validation.nor_quiet(no.nor, mean.vec.nor, var.nor)
  if (is.null(sigma_star)) {
    validation.corr_quiet(no.bin, no.nor, prop.vec.bin, corr.vec = corr.vec, 
                    corr.mat = corr.mat)
    sig_star = compute.sigma.star_quiet(no.bin, no.nor, prop.vec.bin, 
                                  corr.vec, corr.mat)
    sigma_star = sig_star$sigma_star
    if (sig_star$PD == FALSE & continue.with.warning == FALSE) {
      stop("User has chosen to stop as the final correlation matrix is not positive definite")
    }
  }
  else {
    if (is.positive.definite(sigma_star) == FALSE) {
      if (continue.with.warning == TRUE) {
        sigma_star = as.matrix(nearPD(sigma_star, corr = TRUE, 
                                      keepDiag = TRUE)$mat)
        sigma_star = (sigma_star + t(sigma_star))/2
      }
      else {
        stop("The final correlation matrix is not positive definite")
      }
    }
  }
  data = rmvnorm(no.rows, mean = rep(0, d), sigma = sigma_star)
  p = prop.vec.bin
  q = 1 - p
  if (no.bin > 0) {
    for (i in 1:no.rows) {
      for (j in 1:no.bin) {
        if (data[i, j] <= qnorm(1 - p[j])) 
          data[i, j] = 0
        else data[i, j] = 1
      }
    }
  }
  if (no.nor > 0) {
    temp = 1
    for (j in (no.bin + 1):d) {
      data[, j] = mean.vec.nor[temp] + (data[, j] * sqrt(var.nor[temp]))
      temp = temp + 1
    }
  }
  return(data)
}