# === This file is related to===================================================================================================================================
# Testing for similarity of multivariate mixed outcomes using generalised joint regression models with application to efficacy-toxicity responses
# Niklas Hagemann, Giampiero Marra, Frank Bretz and Kathrin Moellenhoff

# === File description =========================================================================================================================================
# This file contains the functions for the simulation study and the case study.

# === Bivariate binary =========================================================================================================================================

marg_prob <- function(beta0,beta1,x){
  return(1/(1+exp(-(beta0+beta1*x))))
} 

parameter_range_condition <- function(v){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_theta  = v[5]
  m2_beta01 = v[6]
  m2_beta11 = v[7]
  m2_beta02 = v[8]
  m2_beta12 = v[9]
  m2_theta  = v[10]
  
  r1 <- m1_theta + 1 
  r2 <- (m1_theta - 1) * (-1)
  r3 <- m2_theta + 1
  r4 <- (m2_theta - 1) * (-1)
  r <- c(r1, r2, r3, r4)
  # if(all(c(r1, r2, r3, r4) > 0)){
  #   return(c(r1, r2, r3, r4))
  # } else {
  #   r <- c(r1, r2, r3, r4)
  #   r[which.min(r)] <- -999
  #   return(r)
  # }
  return(min(r))
}
  
softmax_equal_epsilon <- function(v, model1 = NA, model2 = NA, search_grid, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_theta  = v[5]
  m2_beta01 = v[6]
  m2_beta11 = v[7]
  m2_beta02 = v[8]
  m2_beta12 = v[9]
  m2_theta  = v[10]
  
  m1_mp1 = marg_prob(m1_beta01,m1_beta11, x)
  m1_mp2 = marg_prob(m1_beta02,m1_beta12, x)
  m2_mp1 = marg_prob(m2_beta01,m2_beta11, x)
  m2_mp2 = marg_prob(m2_beta02,m2_beta12, x)
  
  max_abs_diff1 = max(abs(m1_mp1 - m2_mp1))
  max_abs_diff2 = max(abs(m1_mp2 - m2_mp2))
  
  return((max(max_abs_diff1, max_abs_diff2)) - epsilon)
}

model1_likelihood=function(v){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_theta_star = v[5]
  
  logLik1 <- llBin(fo = model1, 
                   par.eq1 = c(m1_beta01, m1_beta11), 
                   par.eq2 = c(m1_beta02, m1_beta12), 
                   ass.par.star = m1_theta_star, 
                   ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  return(-(logLik1))
}

model2_likelihood=function(v){
  m2_beta01 = v[1]
  m2_beta11 = v[2]
  m2_beta02 = v[3]
  m2_beta12 = v[4]
  m2_theta_star  = v[5]
  
  logLik2 <- llBin(fo = model2, 
                   par.eq1 = c(m2_beta01, m2_beta11), 
                   par.eq2 = c(m2_beta02, m2_beta12), 
                   ass.par.star = m2_theta_star,
                   ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  return(-(logLik2))
}

joint_likelihood=function(v, model1, model2, search_grid = NA, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_theta_star = v[5]
  m2_beta01 = v[6]
  m2_beta11 = v[7]
  m2_beta02 = v[8]
  m2_beta12 = v[9]
  m2_theta_star = v[10]
  
  logLik1 <- llBin(fo = model1, 
                   par.eq1 = c(m1_beta01, m1_beta11), 
                   par.eq2 = c(m1_beta02, m1_beta12), ass.par.star = m1_theta_star, 
                   ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  logLik2 <- llBin(fo = model2, 
                   par.eq1 = c(m2_beta01, m2_beta11), 
                   par.eq2 = c(m2_beta02, m2_beta12), ass.par.star = m2_theta_star, 
                   ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  return(-(logLik1 + logLik2))
}

llBin <- function(fo, par.eq1, par.eq2, ass.par = NULL, ass.par.star = NULL, ass.par.star_bound){
  
  if(dim(fo$X1)[2] != length(par.eq1)) stop("Check dimension of par.eq1.")
  if(dim(fo$X2)[2] != length(par.eq2)) stop("Check dimension of par.eq2.")
  
  if(all(is.null(ass.par), is.null(ass.par.star))){
    stop("One of ass.par and ass.par.star must be specified.")
  }
  if(! any(is.null(ass.par), is.null(ass.par.star))){
    stop("ass.par and ass.par.star cannot be specified at the same time.")
  }
  
  if(! is.null(ass.par.star)){
    if(! is.numeric(ass.par.star_bound)){
      stop("ass.par.star_bound must be numeric if using ass.par.star.")
    }
    ass.par.star <- ifelse(abs(ass.par.star) > ass.par.star_bound, sign(ass.par.star)*ass.par.star_bound, ass.par.star)
    ass.par <- tanh(ass.par.star)
  }
  
  eta1 <- fo$X1%*%par.eq1
  eta2 <- fo$X2%*%par.eq2
  
  pd1 <- probm(eta1, fo$margins[1], only.pr = FALSE, min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  pd2 <- probm(eta2, fo$margins[2], only.pr = FALSE, min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  
  p1 <- pd1$pr
  p2 <- pd2$pr
  
  Cop <- fo$VC$BivD; if(Cop != "N") stop("Copulae other than Gaussian not allowed.") 
  nC  <- fo$VC$nC 
  
  p11 <- mm(BiCDF(p1, p2, nC, ass.par, fo$VC$dof), min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  
  p10 <- mm(p1 - p11, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  p01 <- mm(p2 - p11, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  p00 <- mm(1 - p11 - p10 - p01, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  
  
  # safety check ##############################
  sall.p <- rowSums(cbind(p11, p10, p01, p00))
  p11 <- p11/sall.p 
  p10 <- p10/sall.p
  p01 <- p01/sall.p
  p00 <- p00/sall.p
  #############################################
  
  l.par <- fo$respvec$y1.y2*log(p11) + fo$respvec$y1.cy2*log(p10) + fo$respvec$cy1.y2*log(p01) + fo$respvec$cy1.cy2*log(p00) 
  
  list(ll = sum(l.par), eta1 = eta1, eta2 = eta2, p1 = p1, p2 = p2)
  
}

gen_data <- function(x, beta01, beta11, beta02, beta12, theta, n, col_names = NULL){
  for(i in 1:length(x)){
    # print (i)
    
    do <- x[i]  
    d1 <- marg_prob(x = do, beta0 = beta01, beta1 = beta11)
    d2 <- marg_prob(x = do, beta0 = beta02, beta1 = beta12)
    
    margp <- c(d1,d2) #marginal probabilities
    
    # if(abs(theta) > 0.99) {print(theta)}
    
    cr <- matrix(c(1,theta,theta,1),nrow=2,ncol=2) #correlation matrix
    
    sim_data <- generate.binary(nObs=n,prop.vec.bin=margp,corr.mat=cr)
    
    sim_data <- cbind(sim_data, rep(do, n))
    
    if (i == 1) {
      res <- sim_data
    } else {
      res <- rbind(res, sim_data)
    }
  }
  
  if (! is.null(col_names)){
    colnames(res) <- col_names
  }
  
  return(res)
}

# === Bivariate continuous =====================================================================================================================================

gen_data_contin <- function(x, beta01, beta11, beta21, beta02, beta12, beta22, variance1, variance2, theta, n, col_names = NULL){
  for(i in 1:length(x)){
    # print (i)
    
    do <- x[i]  
    d1 <- beta01 + beta11 * do + beta21 * do^2
    d2 <- beta02 + beta12 * do + beta22 * do^2
    
    ev <- c(d1,d2) #Expected values 
    
    covariance <- sqrt(variance1) * sqrt(variance2) * theta
    
    S <- matrix(c(variance1,covariance,covariance,variance2),nrow=2,ncol=2) #correlation matrix
    
    sim_data <- mvtnorm::rmvnorm(n = n, mean = ev, sigma = S)
    
    sim_data <- cbind(sim_data, rep(do, n))
    
    if (i == 1) {
      res <- sim_data
    } else {
      res <- rbind(res, sim_data)
    }
  }
  
  if (! is.null(col_names)){
    colnames(res) <- col_names
  }
  
  return(res)
}

llContCont <- function(fo, par.eq1, par.eq2, par.eq3, par.eq4, ass.par = NULL, ass.par.star = NULL, ass.par.star_bound){

  if(dim(fo$X1)[2] != length(par.eq1)) stop("Check dimension of par.eq1.")
  if(dim(fo$X2)[2] != length(par.eq2)) stop("Check dimension of par.eq2.")
  
  if( is.null(fo$X3) && length(par.eq3) != 1) stop("You can only provide one sigma1 parameter.")
  if(!is.null(fo$X3) && dim(fo$X3)[2] != length(par.eq3) ) stop("Check dimension of par.eq3.")
  
  if( is.null(fo$X4) && length(par.eq4) != 1) stop("You can only provide one sigma2 parameter.")
  if(!is.null(fo$X4) && dim(fo$X4)[2] != length(par.eq3) ) stop("Check dimension of par.eq4.")  
  
  if(! is.null(ass.par.star)){
    if(! is.numeric(ass.par.star_bound)){
      stop("ass.par.star_bound must be numeric if using ass.par.star.")
    }
    ass.par.star <- ifelse(abs(ass.par.star) > ass.par.star_bound, sign(ass.par.star)*ass.par.star_bound, ass.par.star)
    ass.par <- tanh(ass.par.star)
  }
  
  eta1 <- fo$X1%*%par.eq1
  eta2 <- eta.tr(fo$X2%*%par.eq2, fo$margins[2])
  
  if(is.null(fo$X3))  sigma.st1 <- etas1 <- par.eq3 
  if(!is.null(fo$X3)) sigma.st1 <- etas1 <- fo$X3%*%par.eq3  
  
  if(is.null(fo$X4))  sigma.st2 <- etas2 <- par.eq4 
  if(!is.null(fo$X4)) sigma.st2 <- etas2 <- fo$X4%*%par.eq4    
  
  sstr1  <- esp.tr(sigma.st1, fo$margins[1])  
  sigma1 <- sstr1$vrb 
  
  sstr2  <- esp.tr(sigma.st2, fo$margins[2])  
  sigma2 <- sstr2$vrb   
  
  dHs  <- distrHsAT(fo$y1, eta1, sigma1, nu = 1, margin2 = fo$margins[1], min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  pdf1 <- dHs$pdf2
  p1   <- dHs$p2   
  
  dHs  <- distrHsAT(fo$y2, eta2, sigma2, nu = 1, margin2 = fo$margins[2], min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  pdf2 <- dHs$pdf2
  p2   <- dHs$p2 
  
  c.copula2.be1be2 <- copgHsAT(p1, p2, ass.par, fo$VC$BivD, Ln = TRUE, min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)$c.copula2.be1be2
  
  l.par <- log(pdf1) + log(pdf2) + log(c.copula2.be1be2) 
  
  list(ll = sum(l.par) )
}

joint_likelihood_contin = function(v, model1, model2, search_grid = NA, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_sigma1.star = v[5]
  m1_sigma2.star = v[6]
  m1_theta_star = v[7]
  
  m2_beta01 = v[8]
  m2_beta11 = v[9]
  m2_beta21 = v[10]
  m2_beta02 = v[11]
  m2_beta12 = v[12]
  m2_beta22 = v[13]
  m2_sigma1.star = v[14]
  m2_sigma2.star = v[15]
  m2_theta_star = v[16]

  logLik1 <- llContCont(fo = model1,
                        par.eq1 = c(m1_beta01, m1_beta11), 
                        par.eq2 = c(m1_beta02, m1_beta12), 
                        par.eq3 = m1_sigma1.star, 
                        par.eq4 = m1_sigma2.star,
                        ass.par.star = m1_theta_star, 
                        ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  logLik2 <- llContCont(fo = model2,
                        par.eq1 = c(m2_beta01, m2_beta11, m2_beta21),
                        par.eq2 = c(m2_beta02, m2_beta12, m2_beta22), 
                        par.eq3 = m2_sigma1.star, 
                        par.eq4 = m2_sigma2.star,
                        ass.par.star = m2_theta_star, 
                        ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  return(-(logLik1 + logLik2))
}

softmax_contin <- function(v, model1 = NA, model2 = NA, search_grid, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_sigma1.star = v[5]
  m1_sigma2.star = v[6]
  m1_theta_star = v[7]
  
  m2_beta01 = v[8]
  m2_beta11 = v[9]
  m2_beta21 = v[10]
  m2_beta02 = v[11]
  m2_beta12 = v[12]
  m2_beta22 = v[13]
  m2_sigma1.star = v[14]
  m2_sigma2.star = v[15]
  m2_theta_star = v[16]
  
  m1_mp1 <- m1_beta01 + m1_beta11 * search_grid                                 
  m1_mp2 <- m1_beta02 + m1_beta12 * search_grid 
  m2_mp1 <- m2_beta01 + m2_beta11 * search_grid + m2_beta21 * search_grid^2
  m2_mp2 <- m2_beta02 + m2_beta12 * search_grid + m2_beta22 * search_grid^2
  
  max_abs_diff1 = max(abs(m1_mp1 - m2_mp1))
  max_abs_diff2 = max(abs(m1_mp2 - m2_mp2))
  
  return((max(max_abs_diff1, max_abs_diff2)) - epsilon)
}

# === Bivariate mixed ==========================================================================================================================================

gen_data_mixed <- function(x, beta01, beta11, beta02, beta12, beta22, variance, theta, n, col_names = NULL){
  for(i in 1:length(x)){
    
    do <- x[i]  
    d1 <- marg_prob(x = do, beta0 = beta01, beta1 = beta11)
    d2 <- beta02 + beta12 * do + beta22 * do^2
     
    cmat <- matrix(c(1, theta, theta, 1), nrow = 2)
    
    Sstar = compute.sigma.star_quiet(no.bin = 1, no.nor = 1, prop.vec.bin = d1, corr.mat = cmat)
    
    sim_data <- jointly.generate.binary.normal_quiet(no.rows = n,
                                               no.bin = 1,
                                               no.nor = 1,
                                               prop.vec.bin = d1,
                                               mean.vec.nor = d2,
                                               var.nor = variance,
                                               sigma_star = Sstar$sigma_star, 
                                               continue.with.warning=TRUE)
    
    sim_data <- cbind(sim_data, rep(do, n))
    
    if (i == 1) {
      res <- sim_data
    } else {
      res <- rbind(res, sim_data)
    }
  }
  
  if (! is.null(col_names)){
    colnames(res) <- col_names
  }
  
  return(res)
}

llBinCont <- function(fo, par.eq1, par.eq2, par.eq3, ass.par = NULL, ass.par.star = NULL, ass.par.star_bound = NULL){
  
  if(dim(fo$X1)[2] != length(par.eq1)) stop("Check dimension of par.eq1.")
  if(dim(fo$X2)[2] != length(par.eq2)) stop("Check dimension of par.eq2.")
  
  if( is.null(fo$X3) && length(par.eq3) != 1) stop("You can only provide one sigma parameter.")
  if(!is.null(fo$X3) && dim(fo$X3)[2] != length(par.eq3) ) stop("Check dimension of par.eq3.")
  
  if(! is.null(ass.par.star)){
    if(! is.numeric(ass.par.star_bound)){
      stop("ass.par.star_bound must be numeric if using ass.par.star.")
    }
    ass.par.star <- ifelse(abs(ass.par.star) > ass.par.star_bound, sign(ass.par.star)*ass.par.star_bound, ass.par.star)
    ass.par <- tanh(ass.par.star)
  }
  
  eta1 <- fo$X1%*%par.eq1
  eta2 <- eta.tr(fo$X2%*%par.eq2, fo$margins[2])
  
  if(is.null(fo$X3))  sigma.st <- etas <- par.eq3 
  if(!is.null(fo$X3)) sigma.st <- etas <- fo$X3%*%par.eq3  
  
  sstr1 <- esp.tr(sigma.st, fo$margins[2])  
  sigma <- sstr1$vrb 
  
  dHs  <- distrHsAT(fo$y2, eta2, sigma, nu = 1, margin2 = fo$margins[2], min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  pdf2 <- dHs$pdf2
  p2   <- dHs$p2 
  
  pd1 <- probm(eta1, fo$margins[1], bc = TRUE, min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr) 
  p1  <- 1 - pd1$pr 
  
  Cop <- fo$VC$BivD 
  
  dH1 <- copgHs(p1, p2, eta1 = NULL, eta2 = NULL, ass.par, ass.par, Cop, fo$dof, min.dn = fo$VC$min.dn, min.pr = fo$VC$min.pr, max.pr = fo$VC$max.pr)
  h   <- dH1$c.copula.be2 
  
  l.par <- fo$respvec$cy*log(h) + fo$respvec$y1*log(1 - h) + log(pdf2)
  
  list(ll = sum(l.par) )
  
}

joint_likelihood_mixed = function(v, model1, model2, search_grid = NA, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_sigma.star = v[5]
  m1_theta_star = v[6]
  
  m2_beta01 = v[7]
  m2_beta11 = v[8]
  m2_beta02 = v[9]
  m2_beta12 = v[10]
  m2_beta22 = v[11]
  m2_sigma.star = v[12]
  m2_theta_star = v[13]
  
  logLik1 <- llBinCont(fo = model1,
                       par.eq1 = c(m1_beta01, m1_beta11), 
                       par.eq2 = c(m1_beta02, m1_beta12), 
                       par.eq3 = m1_sigma.star, 
                       ass.par.star = m1_theta_star, 
                       ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  logLik2 <- llBinCont(fo = model2,
                       par.eq1 = c(m2_beta01, m2_beta11),
                       par.eq2 = c(m2_beta02, m2_beta12, m2_beta22), 
                       par.eq3 = m2_sigma.star,
                       ass.par.star = m2_theta_star, 
                       ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  return(-(logLik1 + logLik2))
}

softmax_mixed <- function(v, model1 = NA, model2 = NA, search_grid, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_sigma.star = v[5]
  m1_theta_star = v[6]
  
  m2_beta01 = v[7]
  m2_beta11 = v[8]
  m2_beta02 = v[9]
  m2_beta12 = v[10]
  m2_beta22 = v[11]
  m2_sigma.star = v[12]
  m2_theta_star = v[13]

  m1_mp1 <- marg_prob(m1_beta01, m1_beta11, x)                              
  m1_mp2 <- m1_beta02 + m1_beta12 * search_grid
  m2_mp1 <- marg_prob(m2_beta01, m2_beta11, x)
  m2_mp2 <- m2_beta02 + m2_beta12 * search_grid + m2_beta22 * search_grid^2
  
  max_abs_diff1 = max(abs(m1_mp1 - m2_mp1))
  max_abs_diff2 = max(abs(m1_mp2 - m2_mp2))
  
  return((max(max_abs_diff1, max_abs_diff2)) - epsilon)
}

# === Case study ===============================================================================================================================================

joint_likelihood_casestudy = function(v, model1, model2, search_grid = NA, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_beta22 = v[5]
  m1_sigma.star = v[6]
  m1_theta_star = v[7]
  
  m2_beta01 = v[8]
  m2_beta11 = v[9]
  m2_beta02 = v[10]
  m2_beta12 = v[11]
  m2_beta22 = v[12]
  m2_sigma.star = v[13]
  m2_theta_star = v[14]
  
  logLik1 <- llBinCont(fo = model1,
                       par.eq1 = c(m1_beta01, m1_beta11), 
                       par.eq2 = c(m1_beta02, m1_beta12, m1_beta22), 
                       par.eq3 = m1_sigma.star, 
                       ass.par.star = m1_theta_star, 
                       ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  logLik2 <- llBinCont(fo = model2,
                       par.eq1 = c(m2_beta01, m2_beta11),
                       par.eq2 = c(m2_beta02, m2_beta12, m2_beta22), 
                       par.eq3 = m2_sigma.star,
                       ass.par.star = m2_theta_star, 
                       ass.par.star_bound = round(atanh(0.9999), 5))$ll
  
  return(-(logLik1 + logLik2))
}

softmax_casestudy <- function(v, model1 = NA, model2 = NA, search_grid, epsilon){
  m1_beta01 = v[1]
  m1_beta11 = v[2]
  m1_beta02 = v[3]
  m1_beta12 = v[4]
  m1_beta22 = v[5]
  m1_sigma.star = v[6]
  m1_theta_star = v[7]
  
  m2_beta01 = v[8]
  m2_beta11 = v[9]
  m2_beta02 = v[10]
  m2_beta12 = v[11]
  m2_beta22 = v[12]
  m2_sigma.star = v[13]
  m2_theta_star = v[14]

  m1_mp1 <- marg_prob(m1_beta01, m1_beta11, x)                              
  m1_mp2 <- m1_beta02 + m1_beta12 * search_grid + m1_beta22 * search_grid^2
  m2_mp1 <- marg_prob(m2_beta01, m2_beta11, x)
  m2_mp2 <- m2_beta02 + m2_beta12 * search_grid + m2_beta22 * search_grid^2
  
  max_abs_diff1 = max(abs(m1_mp1 - m2_mp1))
  max_abs_diff2 = max(abs(m1_mp2 - m2_mp2))
  
  return((max(max_abs_diff1, max_abs_diff2)) - epsilon)
}
