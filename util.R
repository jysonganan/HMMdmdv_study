library(HMMdmdv)

check_convergence <- function(input_data, n, n1, n2, niter = 1000){
  data_df <- remove_case3(dat = input_data, mean_thresholdPV = 0.1, var_thresholdPV = 0.05, n1 = n1, n2 = n2)
  init_para_est <- init_est(dat_df = data_df, mean_thresholdPV = 0.05, var_thresholdPV = 0.1, n1 = n1, n2 = n2, n = n)
  indep_para_est <- runEM(dat_df = data_df, n1 = n1, n2 = n2, init_para_est= init_para_est, niter = niter)
  emissions <- emission_probs(indep_para_est[5], indep_para_est[6], indep_para_est[7], indep_para_est[8], example_data, n1 = n1, n2 = n2)
  HMM_resList <- runHMM(emissions)
  
  #result <- posterior_inference(HMM_resList[[1]], HMM_resList[[2]], train = FALSE)
}


HMM_resList <- runHMM_iters(emissions)

## differential mean z=1,3 vs z = 2,4
posterior_inference_FDR <- function(gamma, fdr_threshold = 0.1){
  pos_state_prob <- gamma
  pos_prob_null <- apply(pos_state_prob, 1, function(x){return(x[1]+x[3])})
  k = min(pos_prob_null)
  while (k <= 1){
    decision <- sum(pos_prob_null[pos_prob_null <= k])/sum(pos_prob_null <= k)
    if (decision <= fdr_threshold){
      k <- k + 0.01
    }
    else{
      break
    }
  }
  rejection_bool <- (pos_prob_null <= k)
  return(list(rejection_bool, k)) ## true: rejection list
}

post_beta <- posterior_inference_FDR(HMM_resList[[6]], fdr_threshold = 0.1)
## evaluation
TP <- sum((Z == 2) & (Z == 4) & (post_beta[[1]]))
FN <- sum((Z == 2) & (Z == 4) & (!post_beta[[1]]))
Sensitivity <- TP/(TP+FN) #TPR = TP/(TP+FN)





############################################################
#### iterations, also output gamma = posterior probabilities
## rewrite runHMM--return(list(alpha, beta, iter, initPI, Tmat))
runHMM_iters <- function(emissions){
  post <- emissions
  initPI <- rep(0.25, 4)
  Tmat <- matrix(0.25, 4, 4)
  iter <- 0
  oldlogProb = 0
  logProb = -Inf
  c <- rep(NA, n)
  alpha <- matrix(NA, n, 4)
  beta <- matrix(NA, n, 4)
  gamma <- matrix(NA, n, 4)
  bigamma <- array(NA, c(4, 4, n))
  
  initPI_iters <- list()
  Tmat_iters <- list()
  logProb_iters <- rep(NA, niter)
  
  while ((iter == 0 | logProb > oldlogProb) & iter <= niter){
    
    oldlogProb <- logProb
    c[1] <- 0
    for (i in 1:4){
      alpha[1,i] <- initPI[i] * post[1,i]
      c[1] <- c[1] + alpha[1,i]
    }
    
    #c[1] <- 1/c[1]
    c[1] <- 1/c[1]
    for (i in 1:4){
      alpha[1,i] <- c[1] * alpha[1,i]
    }
    
    for (t in 2:n){
      c[t] <- 0
      for (i in 1:4){
        alpha[t,i] <- 0
        for (j in 1:4){
          alpha[t,i] <- alpha[t,i] + alpha[t-1,j] * Tmat[j,i]
        }
        alpha[t,i] <- alpha[t,i] * post[t,i]
        c[t] <- c[t] + alpha[t,i]
      }
      
      #c[t] <- 1/c[t]
      c[t] <- 1/c[t]
      for (i in 1:4){
        alpha[t,i] = c[t] * alpha[t,i]
      }
    }
    
    for (i in 1:4){
      beta[n,i] <- c[n]
    }
    
    for (t in (n-1):1){
      for (i in 1:4){
        beta[t,i] <- 0
        for (j in 1:4){
          beta[t,i] <- beta[t,i] + Tmat[i,j] * post[t+1,j] * beta[t+1,j]
        }
        beta[t,i] <- c[t] * beta[t,i]
      }
    }
    
    for (t in 1:(n-1)){
      denom <- 0
      for (i in 1:4){
        for (j in 1:4){
          denom <- denom + alpha[t,i] * Tmat[i,j] * post[t+1,j] * beta[t+1,j]
        }
      }
      for (i in 1:4){
        gamma[t,i] <- 0
        for (j in 1:4){
          bigamma[i,j,t] <- (alpha[t,i] * Tmat[i,j] * post[t+1,j] * beta[t+1,j])/denom
          gamma[t,i] <- gamma[t,i] + bigamma[i,j,t]
        }
      }
    }
    
    denom <- 0
    for (i in 1:4){
      denom <- denom + alpha[n,i]
    }
    for (i in 1:4){
      gamma[n,i] <- alpha[n,i]/denom
    }
    
    for (i in 1:4){
      initPI[i] <- gamma[1,i]
    }
    for (i in 1:4){
      for (j in 1:4){
        numer <- 0
        denom <- 0
        for (t in 1:(n-1)){
          numer <- numer + bigamma[i,j,t]
          denom <- denom + gamma[t,i]
        }
        Tmat[i,j] <- numer/denom
      }
    }
    
    logProb <- 0
    for (t in 1:n){
      logProb <- logProb + log(c[t])
    }
    logProb = -logProb
    iter <- iter + 1
    
    ####
    initPI_iters[[iter]] <- initPI
    Tmat_iters[[iter]] <- Tmat
    logProb_iters[iter] <- logProb
    
  }
  print('HMM initial state probability estimate:')
  print(initPI)
  print('HMM transition probability estimate:')
  print(Tmat)
  print('HMM iterations:')
  print(iter)
  return(list(alpha, beta, iter, initPI, Tmat, gamma, initPI_iters, Tmat_iters, logProb_iters))
}

