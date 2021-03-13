source("util.R")


library(HMMdmdv)

n<-1000
n1 <- 50
n2 <- 50

data_df <- remove_case3(dat = sim_data, mean_thresholdPV = 0.1, var_thresholdPV = 0.05, n1 = n1, n2 = n2)
init_para_est <- init_est(dat_df = data_df, mean_thresholdPV = 0.05, var_thresholdPV = 0.1, n1 = n1, n2 = n2, n = n)

niter<-1000
indep_para_est <- runEM(dat_df = data_df, n1 = n1, n2 = n2, init_para_est= init_para_est, niter = niter)

emissions <- emission_probs(indep_para_est[5], indep_para_est[6], indep_para_est[7], indep_para_est[8], example_data, n1 = n1, n2 = n2)
head(emissions, 10)


HMM_resList <- runHMM_iters(emissions)
post_beta <- posterior_inference_FDR(HMM_resList[[6]], fdr_threshold = 0.1)


## evaluation
TP <- sum((Z == 2) & (Z == 4) & (post_beta[[1]]))
FN <- sum((Z == 2) & (Z == 4) & (!post_beta[[1]]))
Sensitivity <- TP/(TP+FN) #TPR = TP/(TP+FN)



## max posterior prob; clustering
result <- posterior_inference(HMM_resList[[1]], HMM_resList[[2]], train = FALSE)




