# === This file is related to ==================================================
# Testing for similarity of multivariate mixed outcomes using generalised joint regression models with application to efficacy-toxicity responses
# Niklas Hagemann, Giampiero Marra, Frank Bretz and Kathrin Moellenhoff

# === File description =========================================================
# This file contains the code for the analysis of the case study. 
# For technical reasons the binary response needs to be the first one. Therefore, using the paper notation the parameter vector is (beta_02, beta_12, beta_01, beta_11, beta_21, sigma, rho).

# === Setup (SET ALL THESE MANUALLY) ===========================================
set.seed(12345678) 

# epsilon = 0.15
B <- 300

# === Packages and functions ===================================================

library(alabama) # package for constrained optimization
library(BinNor) # package for generation of mixed data
library(GJRM) # package for Generalised Joint Regression Models
library(dplyr)

source("EffTox_Functions.R")
load(file = "case_study_data.Rdata")

set.seed(12345678) 


# === Simulation ===============================================================

formulae1 <- list(tox ~ dose, eff ~ dose + I(dose^2))
formulae2 <- formulae1

doses <- c(0,0.05,0.2,0.5,1)
doses2 <- c(0,0.1,0.3,0.6,1)

x <- search_grid <- seq(min(doses), max(doses), 0.005)

df1 <- case_study_1
df2 <- case_study_2

model1 <- gjrm(formulae1, data = df1, margins = c("logit","N"), Model = "B", BivD = "N")
model2 <- gjrm(formulae2, data = df2, margins = c("logit","N"), Model = "B", BivD = "N")

param1 <- c(model1$coefficients[1:6], model1$theta)
param2 <- c(model2$coefficients[1:6], model2$theta)

param <- c(param1, param2)


# Max. abs. distance of unconstrained models ===================================

model1_E_mp <- param1[3] + param1[4] * search_grid + param1[5] * search_grid^2 
model2_E_mp <- param2[3] + param2[4] * search_grid + param2[5] * search_grid^2 

m_E <- max(abs(model1_E_mp - model2_E_mp))

model1_T_mp <- marg_prob(beta0 = param1[1], beta1 = param1[2], x = search_grid)
model2_T_mp <- marg_prob(beta0 = param2[1], beta1 = param2[2], x = search_grid)

m_T <- max(abs(model1_T_mp - model2_T_mp))

# Constrained optimization =====================================================

max_dist <- max(m_E, m_T) # maximum of maxima

if (max_dist >= epsilon) {
  minimum <- param
} else {
  
  param[7] <- atanh(param[7])
  param[14] <- atanh(param[14])
  
  minimum_obj <- alabama::auglag(par = param,
                                 fn = joint_likelihood_casestudy,
                                 heq = softmax_casestudy,
                                 control.outer=list(method="nlminb", trace=FALSE),
                                 model1 = model1, 
                                 model2 = model2, 
                                 search_grid = search_grid, 
                                 epsilon = epsilon)
  
  
  minimum  <- as.vector(minimum_obj$par) 
  
  minimum[7] <- tanh(minimum[7])
  minimum[14] <- tanh(minimum[14])
}

cor1 <- df1 %>% dplyr::group_by(dose) %>% summarise(corr = cor(eff, tox)) %>% .$corr %>% mean(., na.rm = TRUE)
cor2 <- df2 %>% dplyr::group_by(dose) %>% summarise(corr = cor(eff, tox)) %>% .$corr %>% mean(., na.rm = TRUE)

# Boostrap =====================================================================

boot_res <- rep(NA, B)
boot_res_E <- rep(NA, B)
boot_res_T <- rep(NA, B)

k <- 1

while(k <= B) {
  
  break_loop_2 <- FALSE
  break_loop_3 <- FALSE
  break_loop_4 <- FALSE
  
  tryCatch(
    {
      step_data1 <- as.data.frame(gen_data_mixed(
        x = doses, 
        n = 30, 
        beta01 = minimum[1], 
        beta11 = minimum[ 2], 
        beta02 = minimum[ 3], 
        beta12 = minimum[ 4], 
        beta22 = minimum[ 5], 
        variance = exp(minimum[ 6])^2, 
        theta = cor1,
        col_names = c("tox", "eff", "doses")
      ))
      
      step_data2 <- as.data.frame(gen_data_mixed(
        x = doses2, 
        n = 30, 
        beta01 = minimum[8], 
        beta11 = minimum[9], 
        beta02 = minimum[10], 
        beta12 = minimum[11], 
        beta22 = minimum[12], 
        variance = exp(minimum[13])^2, 
        theta = cor2,
        col_names = c("tox", "eff", "doses")
      ))
    },
    error = function(e) {
      break_loop_2 <<- TRUE
      cat(paste("Error in gen_data (step_data), i =", i, "\n"))
      cat("ERROR :", conditionMessage(e), "\n")
    }
  )
  
  if (break_loop_2) {
    break
  }
  
  s_step <- min(sd(step_data1[, 1]), sd(step_data1[, 2]), sd(step_data2[, 1]), sd(step_data2[, 2]))
  
  maxCor_step <- max(abs(c(cor(step_data1[, 1], step_data1[, 2]), cor(step_data2[, 1], step_data2[, 2]))))
  
  break_loop_3 <- FALSE
  
  if (s_step == 0 | maxCor_step == 1) {
    
    while_timer <- 0
    
    while (s_step == 0 | maxCor_step == 1) {
      while_timer <- while_timer + 1
      
      if(while_timer >= 1000){
        break_loop_3 <- TRUE
      }
      
      tryCatch(
        {
          step_data1 <- as.data.frame(gen_data_mixed(
            x = doses, 
            n = N, 
            beta01 = minimum[1], 
            beta11 = minimum[ 2], 
            beta02 = minimum[ 3], 
            beta12 = minimum[ 4], 
            beta22 = 0, 
            variance = exp(minimum[ 5])^2, 
            theta = cor1,
            col_names = c("tox", "eff", "doses")
          ))
          
          step_data2 <- as.data.frame(gen_data_mixed(
            x = doses, 
            n = N, 
            beta01 = minimum[7], 
            beta11 = minimum[8], 
            beta02 = minimum[9], 
            beta12 = minimum[10], 
            beta22 = minimum[11], 
            variance = exp(minimum[12])^2, 
            theta = cor2,
            col_names = c("tox", "eff", "doses")
          ))
        },
        error = function(e) {
          break_loop_3 <<- TRUE
          cat(paste("Error in gen_data (step_data), i =", i, "\n"))
          cat("ERROR :", conditionMessage(e), "\n")
        }
      )
      
      if (break_loop_3) {
        break
      }
      
      s_step <- min(sd(step_data1[, 1]), sd(step_data1[, 2]), sd(step_data2[, 1]), sd(step_data2[, 2]))
      maxCor_step <- max(abs(c(cor(step_data1[, 1], step_data1[, 2]), cor(step_data2[, 1], step_data2[, 2]))))
    }
  }
  
  if (break_loop_3) {
    break
  }
  
  names(step_data1) <- c("tox", "eff", "dose")
  names(step_data2) <- c("tox", "eff", "dose")
  
  break_loop_4 <- FALSE
  
  tryCatch(
    {
      step_model1 <- gjrm(formulae1, data = step_data1, margins = c("logit","N"), Model = "B", BivD = "N")
      step_model2 <- gjrm(formulae2, data = step_data2, margins = c("logit","N"), Model = "B", BivD = "N")
    },
    error = function(e) {
      break_loop_4 <<- TRUE
      cat(paste("Error in gjrm (step_data), i =", "\n"))
      cat("ERROR :", conditionMessage(e), "\n")
    }
  )
  
  if (break_loop_4) {
    break
  }
  
  step_param1 <- c(step_model1$coefficients[1:6], step_model1$theta)
  step_param2 <- c(step_model2$coefficients[1:6], step_model2$theta)
  
  step_model1_E_mp <- step_param1[3] + step_param1[4] * search_grid + step_param1[5] * search_grid^2   
  step_model2_E_mp <- step_param2[3] + step_param2[4] * search_grid + step_param2[5] * search_grid^2 
  
  step_m_E <- max(abs(step_model1_E_mp - step_model2_E_mp))
  
  step_model1_T_mp <- marg_prob(beta0 = step_param1[1], beta1 = step_param1[2], x = x)
  step_model2_T_mp <- marg_prob(beta0 = step_param2[1], beta1 = step_param2[2], x = x)
  
  step_m_T <- max(abs(step_model1_T_mp - step_model2_T_mp))
  
  boot_res[k] <- max(step_m_E, step_m_T)
  boot_res_E[k] <- step_m_E
  boot_res_T[k] <- step_m_T
  
  k <- k+1
}

if(any(break_loop_2, break_loop_3, break_loop_4)){
  print("NEXT")
  next
}

ecdf(boot_res)(max_dist)
quantile(boot_res, 0.05, na.rm = TRUE)



