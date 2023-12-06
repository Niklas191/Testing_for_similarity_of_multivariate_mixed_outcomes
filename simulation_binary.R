# === This file is related to =================================================================================================================================
# Testing for similarity of multivariate mixed outcomes using generalised joint regression models with application to efficacy-toxicity responses
# Niklas Hagemann, Giampiero Marra, Frank Bretz and Kathrin Moellenhoff

# === File description ========================================================================================================================================
# This file contains the code for the simulation of the bivarite binary case. 

# === Setup (SET ALL THESE MANUALLY) ==========================================================================================================================
# set.seed()

# epsilon =
# m1_beta01 =
# m1_beta11 =
# m1_beta02 =
# m1_beta12 =
# m2_beta01 =
# m2_beta11 =
# m2_beta02 =
# m2_beta12 =
# Cor =

# N = # Number of observations per dose level per group, i.e. n_g^{(l)}

simu <- 1000
B <- 300

# === Packages and functions ==================================================================================================================================
library(alabama)
library(MultiOrd) # package for generation of binary data
library(GJRM) # package for Generalised Joint Regression Models
library(dplyr)

source("EffTox_Functions.R")

# === Simulation ==============================================================================================================================================
formulae <- list(eff ~ doses, tox ~ doses)

doses <- c(0,0.1,0.2,0.5,1,1.5,2)
x <- search_grid <- seq(min(doses), max(doses), 0.005)

corr <- Cor
  
v <- rep(NA, 10)
    v[1] = m1_beta01
    v[2] = m1_beta11
    v[3] = m1_beta02
    v[4] = m1_beta12
    v[5] = Cor
    v[6] = m2_beta01
    v[7] = m2_beta11
    v[8] = m2_beta02
    v[9] = m2_beta12
    v[10] = Cor
  
  Max_dist <- rep(NA, simu)
  Crit_val <- rep(NA, simu)
  p_Value <- rep(NA, simu)
  NAs <- rep(0, simu)
  softmax_control <- rep(NA, simu)
  maxdist_larger_epsilon <- rep(NA, simu)
  
  i <- 1
  
  while(i <= simu) {
    if (i %in% c(1, 2, 3, seq(0, 1000, 10))) {
      print(i)
    }
    
    skip_to_next <- FALSE
    
    tryCatch(
      {
        df1 <- as.data.frame(gen_data(x = doses, n = N, beta01 = v[1], beta11 = v[2], beta02 = v[3], beta12 = v[4], theta = v[5]))
        df2 <- as.data.frame(gen_data(x = doses, n = N, beta01 = v[6], beta11 = v[7], beta02 = v[8], beta12 = v[9], theta = v[10]))
      },
      error = function(e) {
        skip_to_next <<- TRUE
        cat(paste("Error in gen_data, i =", i, "\n"))
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
    
    if (skip_to_next) {
      print("NEXT")
      next
    }
    
    s <- min(sd(df1[, 1]), sd(df1[, 2]), sd(df2[, 1]), sd(df2[, 2]))
    maxCor <- max(abs(c(cor(df1[, 1], df1[, 2]), cor(df2[, 1], df2[, 2]))))
    
    break_loop_1 <- FALSE
    
    if (s == 0 | maxCor == 1) {
      
      while (s == 0 | maxCor == 1) {

        tryCatch(
          {
            df1 <- as.data.frame(gen_data(x = doses, n = N, beta01 = v[1], beta11 = v[2], beta02 = v[3], beta12 = v[4], theta = v[5]))
            df2 <- as.data.frame(gen_data(x = doses, n = N, beta01 = v[6], beta11 = v[7], beta02 = v[8], beta12 = v[9], theta = v[10]))
          },
          error = function(e) {
            break_loop_1 <<- TRUE
            cat(paste("Error in gen_data, i =", i, "\n"))
            cat("ERROR :", conditionMessage(e), "\n")
          }
        )
        
        if (break_loop_1) {
          break
        }
        
        s <- min(sd(df1[, 1]), sd(df1[, 2]), sd(df2[, 1]), sd(df2[, 2]))
        maxCor <- max(abs(c(cor(df1[, 1], df1[, 2]), cor(df2[, 1], df2[, 2]))))
      }
    }
    
    if (break_loop_1) {
      print("NEXT")
      next
    }
    
    names(df1) <- c("eff", "tox", "doses")
    names(df2) <- c("eff", "tox", "doses")
    
    tryCatch(
      {
        model1 <- gjrm(formulae, data = df1, margins = c("logit", "logit"), Model = "B", BivD = "N")
        model2 <- gjrm(formulae, data = df2, margins = c("logit", "logit"), Model = "B", BivD = "N")
      },
      error = function(e) {
        skip_to_next <<- TRUE
        cat(paste("Error in gjrm, i =", i, "\n"))
        cat("ERROR :", conditionMessage(e), "\n")
      }
    )
    
    if (skip_to_next) {
      print("NEXT")
      next
    }
    
    param1 <- c(model1$coefficients[1:4], model1$theta)
    param2 <- c(model2$coefficients[1:4], model2$theta)
    
    param <- c(param1, param2)
    
    if (param[5] > 0.999) {
      param[5] <- 0.999
    }
    if (param[5] < -0.999) {
      param[5] <- -0.999
    }
    if (param[10] > 0.999) {
      param[10] <- 0.999
    }
    if (param[10] < -0.999) {
      param[10] <- -0.999
    }
    
    # Max. abs. distance of unconstrained models ==============================================================================================================
    
    model1_E_mp <- marg_prob(beta0 = param1[1], beta1 = param1[2], x = search_grid)
    model2_E_mp <- marg_prob(beta0 = param2[1], beta1 = param2[2], x = search_grid)
    
    m_E <- max(abs(model1_E_mp - model2_E_mp))
    
    model1_T_mp <- marg_prob(beta0 = param1[3], beta1 = param1[4], x = search_grid)
    model2_T_mp <- marg_prob(beta0 = param2[3], beta1 = param2[4], x = search_grid)
    
    m_T <- max(abs(model1_T_mp - model2_T_mp))
    
    # Constrained optimization ================================================================================================================================
    
    max_dist <- max(m_E, m_T) # maximum of maxima
    
    if (max_dist >= epsilon) {
      minimum <- param
      maxdist_larger_epsilon[i] <- TRUE
    } else {
      maxdist_larger_epsilon[i] <- FALSE
  
    tryCatch(
      {
    
        param[5] <- atanh(param[5])
        param[10] <- atanh(param[10])
            
    minimum_obj <- alabama::auglag(par = param,
                      fn = joint_likelihood,
                      heq = softmax_equal_epsilon,
                      control.outer=list(method="nlminb", trace=FALSE),
                      model1 = model1, 
                      model2 = model2, 
                      search_grid = search_grid, 
                      epsilon = epsilon)
        
    },
      error = function(e) {
        skip_to_next <<- TRUE
        cat(paste("Error in auglag, i =", i, "\n"))
        cat("ERROR :", conditionMessage(e), "\n")
      })

    if (skip_to_next) {
      print("NEXT")
      next
    }

    minimum  <- as.vector(minimum_obj$par) 
    
    minimum[5] <- tanh(minimum[5])
    minimum[10] <- tanh(minimum[10])
    
    if(abs(softmax_equal_epsilon(minimum, search_grid = x, epsilon = epsilon)) > 0.01){
      print(paste("PROBLEM: softmax_equal_epsilon(minimum) > 0.01", softmax_equal_epsilon(minimum, search_grid = search_grid, epsilon = epsilon)))
      
    } else if (abs(softmax_equal_epsilon(minimum, search_grid = x, epsilon = epsilon)) > 0.001){
      print(paste("PROBLEM: softmax_equal_epsilon(minimum) > 0.001", softmax_equal_epsilon(minimum, search_grid = search_grid, epsilon = epsilon)))
    }
    }
    
    softmax_control[i] <- abs(softmax_equal_epsilon(minimum, search_grid = x, epsilon = epsilon))
    
    if (minimum[5] > 0.999) {
      minimum[5] <- 0.999
    }
    if (minimum[5] < -0.999) {
      minimum[5] <- -0.999
    }
    if (minimum[10] > 0.999) {
      minimum[10] <- 0.999
    }
    if (minimum[10] < -0.999) {
      minimum[10] <- -0.999
    }

    cor1 <- df1 %>% dplyr::group_by(doses) %>% summarise(corr = cor(eff, tox)) %>% .$corr %>% mean(., na.rm = TRUE)
    cor2 <- df2 %>% dplyr::group_by(doses) %>% summarise(corr = cor(eff, tox)) %>% .$corr %>% mean(., na.rm = TRUE)

    # Boostrap ================================================================================================================================================
    
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
          step_data1 <- as.data.frame(gen_data(
            x = doses,
            beta01 = minimum[1],
            beta11 = minimum[2],
            beta02 = minimum[3],
            beta12 = minimum[4],
            theta = cor1,
            n = N,
            col_names = c("eff", "tox", "doses")
          ))
          
          step_data2 <- as.data.frame(gen_data(
            x = doses,
            beta01 = minimum[6],
            beta11 = minimum[7],
            beta02 = minimum[8],
            beta12 = minimum[9],
            theta = cor2,
            n = N,
            col_names = c("eff", "tox", "doses")
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
      
      s <- min(sd(df1[, 1]), sd(df1[, 2]), sd(df2[, 1]), sd(df2[, 2]))
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
              step_data1 <- as.data.frame(gen_data(
                x = doses,
                beta01 = minimum[1],
                beta11 = minimum[2],
                beta02 = minimum[3],
                beta12 = minimum[4],
                theta = cor1,
                n = N,
                col_names = c("eff", "tox", "doses")
              ))
              
              step_data2 <- as.data.frame(gen_data(
                x = doses,
                beta01 = minimum[6],
                beta11 = minimum[7],
                beta02 = minimum[8],
                beta12 = minimum[9],
                theta = cor2,
                n = N,
                col_names = c("eff", "tox", "doses")
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
      
      names(step_data1) <- c("eff", "tox", "doses")
      names(step_data2) <- c("eff", "tox", "doses")
      
      break_loop_4 <- FALSE
      
      tryCatch(
        {
          step_model1 <- gjrm(formulae, data = step_data1, margins = c("logit", "logit"), Model = "B", BivD = "N")
          step_model2 <- gjrm(formulae, data = step_data2, margins = c("logit", "logit"), Model = "B", BivD = "N")
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
      
      step_param1 <- c(step_model1$coefficients[1:4], step_model1$theta)
      step_param2 <- c(step_model2$coefficients[1:4], step_model2$theta)
      
      step_model1_E_mp <- marg_prob(beta0 = step_param1[1], beta1 = step_param1[2], x = x)
      step_model2_E_mp <- marg_prob(beta0 = step_param2[1], beta1 = step_param2[2], x = x)
      
      step_m_E <- max(abs(step_model1_E_mp - step_model2_E_mp))
      
      step_model1_T_mp <- marg_prob(beta0 = step_param1[3], beta1 = step_param1[4], x = x)
      step_model2_T_mp <- marg_prob(beta0 = step_param2[3], beta1 = step_param2[4], x = x)
      
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
    
    Max_dist[i] <- max_dist
    Crit_val[i] <- quantile(boot_res, 0.05, na.rm = TRUE)
    p_Value[i] <- ecdf(boot_res)(max_dist)
    NAs[i] <- sum(is.na((boot_res)), na.rm = TRUE)
    
    i <- i + 1
  }
  
  saveList <- list(Max_dist, Crit_val, p_Value, NAs, corr, N, param_ind, epsilon, seeds[combi], maxdist_larger_epsilon, softmax_control)
