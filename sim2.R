source("/gpfs/scratch/jusong/util.R")
library(HMMdmdv)

#### simulation (without model setting)

library(locfdr)
library(NHMMfdr)

set.seed(3)
n = 5000
n1 = 100
n2 = 100
fdr = 0.1
nsim = 100


pi1 = 0.4
p = 0.2 #0.5, 0.8

mu0 = 2
sigma0 = 0.35
sigma = 0.1

library(MCMCpack)
truePara_generation_1 <- function(n_state, pi1, p){
  initProb = c(pi1, (1-pi1)/3, (1-pi1)/3, (1-pi1)/3)
  transProb = matrix(NA, n_state, n_state)
  for (i in 1:n_state){transProb[i,i] <- p}
  
  for (i in 1:4){
    otherProb = rdirichlet(1, c(1,1,1))
    if(i == 1){
      transProb[i,2:4] = otherProb[1,] * (1 - transProb[i,i])
    }
    if (i == 2){
      transProb[i,c(1,3,4)] = otherProb[1,] * (1 - transProb[i,i])
    }
    if (i == 3){
      transProb[i,c(1,2,4)] = otherProb[1,] * (1 - transProb[i,i])
    }
    if (i == 4){
      transProb[i,1:3] = otherProb[1,] * (1 - transProb[i,i])
    }
  }
  
  return(list(initProb, transProb))
}



data_generation_1 <- function(d, mu0, sigma0, sigma, n, n1, n2, initProb, transProb){
  X<-matrix(NA, n, n1)
  Y<-matrix(NA, n, n2)
  Z<-rep(NA,n)
  
  Z[1]<-sample(1:4, 1, prob = initProb)
  
  for (i in 2:n){
    Z[i]<-sample(1:4, 1, prob = transProb[Z[i-1],])
  }


  for (i in 1:n){
    if(Z[i]==1)
    { 
      mu <- rnorm(1, mu0, sigma0)
      X[i,]<-rnorm(n1,mu,sqrt(sigma))
      Y[i,]<-rnorm(n2,mu,sqrt(sigma))
    }
    else if(Z[i]==2){
      mu <- rnorm(1, mu0, sigma0)
      X[i,]<-rnorm(n1,mu,sqrt(sigma))
      Y[i,]<-rnorm(n2,mu + d,sqrt(sigma))
    }
    else if(Z[i]==3){
      mu <- rnorm(1, mu0, sigma0)
      X[i,]<-rnorm(n1,mu,sqrt(sigma))
      Y[i,]<-rnorm(n2,mu,sqrt(10*sigma))
    }
    else{
      mu <- rnorm(1, mu0, sigma0)
      X[i,]<-rnorm(n1,mu,sqrt(sigma1))
      Y[i,]<-rnorm(n2,mu+d,sqrt(10*sigma2))
    }
  }

  A<-cbind(X,Y)
  print('data generated!')
  return(list(A, Z))
}





significantList_length <- matrix(NA, nsim, 5)
significantList_length <- as.data.frame(significantList_length)
colnames(significantList_length) <- c("HMMdmdv","BH-FDR", "IndepStat", "HMMStat", "True non-null")

FDR_tab <- matrix(NA, nsim, 4)
FDR_tab <- as.data.frame(FDR_tab)
colnames(FDR_tab) <-  c("HMMdmdv","BH-FDR", "IndepStat", "HMMStat")

FNDR_tab <- matrix(NA, nsim, 4)
FNDR_tab <- as.data.frame(FNDR_tab)
colnames(FNDR_tab) <-  c("HMMdmdv","BH-FDR", "IndepStat", "HMMStat")

Sensitivity_tab <- matrix(NA, nsim, 4)
Sensitivity_tab <- as.data.frame(Sensitivity_tab)
colnames(Sensitivity_tab) <-  c("HMMdmdv","BH-FDR", "IndepStat", "HMMStat")

Specificity_tab<- matrix(NA, nsim, 4)
Specificity_tab <- as.data.frame(Specificity_tab)
colnames(Specificity_tab) <-  c("HMMdmdv","BH-FDR", "IndepStat", "HMMStat")

tran_prob_sumsqErr = nu0_sqErr = var0_sqErr = k0_sqErr = mu0_sqErr = cluster_concordance_acc <- rep(NA, nsim)






for (i in 1:nsim){
  true_para <- truePara_generation_1(4, pi1, p)

  sim_alldata <- data_generation_1(d, mu0, sigma0, sigma, n, n1, n2, true_para[[1]], true_para[[2]])
  sim_data <- sim_alldata[[1]]
  sim_state <- sim_alldata[[2]]
  
  
  
  
  
  
  
  data_df <- remove_case3(dat = sim_data, mean_thresholdPV = 0.1, var_thresholdPV = 0.05, n1 = n1, n2 = n2)
  init_para_est <- init_est(dat_df = data_df, mean_thresholdPV = 0.05, var_thresholdPV = 0.1, n1 = n1, n2 = n2, n = n)
  
  niter<-1000
  indep_para_est <- runEM(dat_df = data_df, n1 = n1, n2 = n2, init_para_est= init_para_est, niter = niter)
  
  emissions <- emission_probs(indep_para_est[5], indep_para_est[6], indep_para_est[7], indep_para_est[8], sim_data, n1 = n1, n2 = n2)
  
  HMM_resList <- runHMM_iters(emissions)
  
  ## MSE parameter estimates
  nu0_sqErr[i] <- (indep_para_est[5] - true_para[[1]])^2
  var0_sqErr[i] <- (indep_para_est[6] - true_para[[2]])^2
  k0_sqErr[i] <- (indep_para_est[7] - true_para[[3]])^2
  mu0_sqErr[i] <- (indep_para_est[8] - true_para[[4]])^2
  tran_prob_sumsqErr[i] <- sum((HMM_resList[[5]] - true_para[[6]])^2)
  
  ## 1. accuracy concordance (clustering) (compare with DM+DV test threholds for 4 states)
  result <- posterior_inference(HMM_resList[[1]], HMM_resList[[2]], Z = sim_state, train = TRUE)
  cluster_concordance_acc[i] <- result[[2]]
  
  ## 2. differential mean FDR control sensitivity.. (comparison)
  post_beta <- posterior_inference_FDR(HMM_resList[[6]], fdr_threshold = fdr)
  
  ### BH-FDR
  ttest_pvalue <- apply(sim_data, 1, function(x){return(t.test(x[1:n1], x[(n1+1):(n1+n2)])$p.value)})
  ttest_pvalue_adj <- p.adjust(ttest_pvalue, "BH")
  # true: rejection list
  BH_rejection_bool <- (ttest_pvalue_adj <= fdr)
  ttest_zvalue <- qnorm(ttest_pvalue)
  
  ### locfdr
  #locfdr_fdr <- locfdr(ttest_zvalue)$fdr
  #locfdr_rejection_bool <- (locfdr_fdr <= fdr) 
  
  ## IndepStat
  ## HMMstat sun method
  IndepStat_model <- fdr.nhmm(ttest_zvalue, modeltype = "Indep", Z = NULL, dist = NULL)
  IndepStat_LIS <- LIS.adjust(IndepStat_model$LIS, fdr = fdr, adjust = TRUE)
  #IndepStat_LIS$States #1: reject, significant
  HMMStat_model <- fdr.nhmm(ttest_zvalue, modeltype = "HMM", Z = NULL, dist = NULL)
  HMMStat_LIS <- LIS.adjust(HMMStat_model$LIS, fdr = fdr, adjust = TRUE)
  
  
  
  print(paste("true non-null list length: ", sum(sim_state == 2|sim_state == 4)))
  print(paste("significant list length of HMMdmdv: ", sum(post_beta[[1]])))
  print(paste("significant list length of BH-FDR: ", sum(BH_rejection_bool)))
  print(paste("significant list length of IndepStat: ", sum(IndepStat_LIS$States)))
  print(paste("significant list length of HMMStat: ", sum(HMMStat_LIS$States)))
  
  significantList_length[i,] <- c(sum(post_beta[[1]]), sum(BH_rejection_bool), sum(IndepStat_LIS$States), 
                                  sum(HMMStat_LIS$States), sum(sim_state == 2|sim_state == 4))
  
  
  FP_HMMdmdv <- sum((post_beta[[1]] & (sim_state == 1|sim_state == 3)))
  FP_BH <- sum((BH_rejection_bool & (sim_state == 1|sim_state == 3)))
  FP_IndepStat <- sum((IndepStat_LIS$States & (sim_state == 1|sim_state == 3)))
  FP_HMMStat <- sum((HMMStat_LIS$States & (sim_state == 1|sim_state == 3)))
  
  TP_HMMdmdv <- sum((post_beta[[1]] & (sim_state == 2|sim_state == 4)))
  TP_BH <- sum((BH_rejection_bool & (sim_state == 2|sim_state == 4)))
  TP_IndepStat <- sum((IndepStat_LIS$States & (sim_state == 2|sim_state == 4)))
  TP_HMMStat <- sum((HMMStat_LIS$States & (sim_state == 2|sim_state == 4)))
  
  FN_HMMdmdv <- sum((!post_beta[[1]] & (sim_state == 2|sim_state == 4)))
  FN_BH <- sum((!BH_rejection_bool & (sim_state == 2|sim_state == 4)))
  FN_IndepStat <- sum((!IndepStat_LIS$States & (sim_state == 2|sim_state == 4)))
  FN_HMMStat <- sum((!HMMStat_LIS$States & (sim_state == 2|sim_state == 4)))
  
  TN_HMMdmdv <- sum((!post_beta[[1]] & (sim_state == 1|sim_state == 3)))
  TN_BH <- sum((!BH_rejection_bool & (sim_state == 1|sim_state == 3)))
  TN_IndepStat <- sum((!IndepStat_LIS$States & (sim_state == 1|sim_state == 3)))
  TN_HMMStat <- sum((!HMMStat_LIS$States & (sim_state == 1|sim_state == 3)))
  
  ## realized FDR 
  FDR_HMMdmdv <- FP_HMMdmdv/(FP_HMMdmdv + TP_HMMdmdv)
  FDR_BH <- FP_BH/(FP_BH + TP_BH)
  FDR_IndepStat <- FP_IndepStat/(FP_IndepStat + TP_IndepStat)
  FDR_HMMStat <- FP_HMMStat/(FP_HMMStat + TP_HMMStat)
  
  ## FNDR
  FNDR_HMMdmdv <- FN_HMMdmdv/(FN_HMMdmdv + TN_HMMdmdv)
  FNDR_BH <- FN_BH/(FN_BH + TN_BH)
  FNDR_IndepStat <- FN_IndepStat/(FN_IndepStat + TN_IndepStat)
  FNDR_HMMStat <- FN_HMMStat/(FN_HMMStat + TN_HMMStat)
  
  
  ## sensitivity
  Sensitivity_HMMdmdv <- TP_HMMdmdv/(FN_HMMdmdv + TP_HMMdmdv)
  Sensitivity_BH <- TP_BH/(FN_BH + TP_BH)
  Sensitivity_IndepStat <- TP_IndepStat/(FN_IndepStat + TP_IndepStat)
  Sensitivity_HMMStat <- TP_HMMStat/(FN_HMMStat + TP_HMMStat)
  
  # specificity
  Specificity_HMMdmdv <- TN_HMMdmdv/(TN_HMMdmdv + FP_HMMdmdv)
  Specificity_BH <- TN_BH/(TN_BH + FP_BH)
  Specificity_IndepStat <- TN_IndepStat/(TN_IndepStat + FP_IndepStat)
  Specificity_HMMStat <- TN_HMMStat/(TN_HMMStat + FP_HMMStat)
  
  
  FDR_tab[i,] <- c(FDR_HMMdmdv, FDR_BH, FDR_IndepStat, FDR_HMMStat)
  FNDR_tab[i,] <- c(FNDR_HMMdmdv, FNDR_BH, FNDR_IndepStat, FNDR_HMMStat)
  Sensitivity_tab[i,] <- c(Sensitivity_HMMdmdv,Sensitivity_BH,Sensitivity_IndepStat,Sensitivity_HMMStat)
  Specificity_tab[i,] <- c(Specificity_HMMdmdv,Specificity_BH,Specificity_IndepStat,Specificity_HMMStat)
}


save("significantList_length", "FDR_tab", "FNDR_tab", "Sensitivity_tab", "Specificity_tab", "cluster_concordance_acc",
     "nu0_sqErr", "var0_sqErr", "k0_sqErr", "mu0_sqErr", "tran_prob_sumsqErr", file = "/gpfs/scratch/jusong/HMMsim_5.RData")

