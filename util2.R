## suppose only case 1,2,4:indepEM + HMM

data_generation_2 <- function(d, mu, sigma_var, n, n1, n2, initProb, transProb, var_ratio = 10){
  X<-matrix(NA, n, n1)
  Y<-matrix(NA, n, n2)
  Z<-rep(NA,n)

  Z[1]<-sample(1:3, 1, prob = initProb)

  for (i in 2:n){
    Z[i]<-sample(1:3, 1, prob = transProb[Z[i-1],])
  }


  for (i in 1:n){
    if(Z[i]==1)
    {
      X[i,]<-rnorm(n1,mu[i],sqrt(sigma_var[i]))
      Y[i,]<-rnorm(n2,mu[i],sqrt(sigma_var[i]))
    }
    else if(Z[i]==2){
      X[i,]<-rnorm(n1,mu[i],sqrt(sigma_var[i]))
      Y[i,]<-rnorm(n2,mu[i] + d,sqrt(sigma_var[i]))
    }
    else{
      X[i,]<-rnorm(n1,mu[i],sqrt(sigma_var[i]))
      Y[i,]<-rnorm(n2,mu[i] + d,sqrt(var_ratio*sigma_var[i]))
    }
  }

  A<-cbind(X,Y)
  print('data generated!')
  return(list(A, Z))
}




data_generation <- function (nu0, var0, k0, mu0, n, n1, n2, initProb, transProb)
{
  X <- matrix(NA, n, n1)
  Y <- matrix(NA, n, n2)
  Z <- rep(NA, n)
  Z[1] <- sample(1:3, 1, prob = initProb)
  for (i in 2:n) {
    Z[i] <- sample(1:3, 1, prob = transProb[Z[i - 1], ])
  }
  for (i in 1:n) {
    if (Z[i] == 1) {
      sigma <- pscl::rigamma(1, nu0/2, nu0 * var0/2)
      mu <- rnorm(1, mu0, sqrt(sigma/k0))
      X[i, ] <- rnorm(n1, mu, sqrt(sigma))
      Y[i, ] <- rnorm(n2, mu, sqrt(sigma))
    }
    else if (Z[i] == 2) {
      sigma <- pscl::rigamma(1, nu0/2, nu0 * var0/2)
      mu1 <- rnorm(1, mu0, sqrt(sigma/k0))
      mu2 <- rnorm(1, mu0, sqrt(sigma/k0))
      X[i, ] <- rnorm(n1, mu1, sqrt(sigma))
      Y[i, ] <- rnorm(n2, mu2, sqrt(sigma))
    }
    else {
      sigma1 <- pscl::rigamma(1, nu0/2, nu0 * var0/2)
      sigma2 <- pscl::rigamma(1, nu0/2, nu0 * var0/2)
      mu1 <- rnorm(1, mu0, sqrt(sigma1/k0))
      mu2 <- rnorm(1, mu0, sqrt(sigma2/k0))
      X[i, ] <- rnorm(n1, mu1, sqrt(sigma1))
      Y[i, ] <- rnorm(n2, mu2, sqrt(sigma2))
    }
  }
  A <- cbind(X, Y)
  print("data generated!")
  return(list(A, Z))
}






f1<-function(x,y,nu,var,k,mu){
  l1<-length(x)
  l2<-length(y)
  nun<-nu+l1+l2
  kn<-k+l1+l2
  xbar<-mean(x)
  ybar<-mean(y)
  xsq<-var(x)*(l1-1)
  ysq<-var(y)*(l2-1)
  nuvarn<-xsq+ysq+l1*xbar^2+l2*ybar^2+nu*var+k*mu^2-(l1*xbar+l2*ybar+k*mu)^2/(l1+l2+k)
  nom<-lgamma(nun/2) + log(sqrt(k)) + log(nu*var/2)*(nu/2)
  denom1<-log(2*pi)*(l1+l2)/2
  denom2<-lgamma(nu/2) + log(sqrt(kn)) + log(nuvarn/2) * (nun/2)
  return(nom-denom1-denom2)
}




f2<-function(x,y,nu,var,k,mu){
  l1<-length(x)
  l2<-length(y)
  xbar<-mean(x)
  ybar<-mean(y)
  xsq<-var(x)*(l1-1)
  ysq<-var(y)*(l2-1)
  nun<-nu+l1+l2
  w1<-(l1*xbar+k*mu)^2/(l1+k)
  w2<-(l2*ybar+k*mu)^2/(l2+k)
  nuvarn<-xsq+ysq+nu*var+l1*xbar^2+l2*ybar^2+2*k*mu^2-w1-w2
  nom<-log(k)+lgamma(nun/2) + log(nu*var/2)*(nu/2)
  denom2<-log(sqrt(l1+k))+log(sqrt(l2+k))+lgamma(nu/2)+ log(nuvarn/2)*(nun/2)
  denom1<-log(2*pi)*(l1+l2)/2
  return(nom-denom1-denom2)
}




f3<-function(x,y,nu,var,k,mu){
  l1<-length(x)
  l2<-length(y)
  xbar<-mean(x)
  ybar<-mean(y)
  xsq<-var(x)*(l1-1)
  ysq<-var(y)*(l2-1)
  nun1<-nu+l1
  nun2<-nu+l2
  nuvarn1<-nu*var+xsq+l1*k*(xbar-mu)^2/(l1+k)
  nuvarn2<-nu*var+ysq+l2*k*(ybar-mu)^2/(l2+k)
  nom<-log(k)+log(nu*var/2)*nu+lgamma(nun1/2)+lgamma(nun2/2)
  denom1<-log(2*pi)*(l1+l2)/2
  denom2<-log(sqrt(l1+k))+log(sqrt(l2+k))+lgamma(nu/2)*2+log(nuvarn1/2)*(nun1/2)+log(nuvarn2/2)*(nun2/2)
  return(nom-denom1-denom2)
}




emission_probs <- function(nu, var, k, mu, A, n1, n2){

  data_df <- data.frame(A)

  emission_analytical<-apply(data_df,1,function(x){
    a<-f1(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
    b<-f2(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
    d<-f3(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
    return(c(a,b,d))
  })
  emission_analytical <- t(emission_analytical)
  return(emission_analytical)
}






init_est <- function(dat_df, mean_thresholdPV, var_thresholdPV, n1, n2, n){
  # input dat_df (case 3 removed)
  # threshold to subset case1 data
  pv <- apply(dat_df, 1, function(x) t.test(x[1:n1], x[(n1+1):(n1+n2)])$p.value)
  qv <- apply(dat_df, 1, function(x){
    l <- list(x[1:n1], x[(n1+1):(n1+n2)])
    bartlett.test(l)$p.value
  })
  dat_df$pv <- pv
  dat_df$qv <- qv
  case1 <- subset(dat_df, pv >= mean_thresholdPV & qv >= var_thresholdPV)
  case1 <- case1[, -c(n1+n2+1, n1+n2+2)]

  length1 <- dim(case1)[1]

  # use the average values in case 1 dataset for initial hyperparameter estimate
  mu_1 <- apply(case1, 1, mean)
  mu_1 <- as.numeric(mu_1)
  mu_0 <- mean(mu_1)

  mu_1_var <- var(mu_1)
  sigma_1 <- apply(case1, 1, var)
  sigma_1 <- as.numeric(sigma_1)
  sigma_1 <- sigma_1*(n1+n2-1)/(n1+n2-2)
  sigma_1_mean <- mean(sigma_1)
  sigma_1_var <- var(sigma_1)

  k_0 <- sigma_1_mean/mu_1_var
  nu_0 <- 4 + 2*sigma_1_mean^2/sigma_1_var
  var_0 <- (sigma_1_mean^3 + sigma_1_mean*sigma_1_var)/(sigma_1_mean^2 + 2*sigma_1_var)

  case2 <- subset(dat_df, pv <= mean_thresholdPV & qv >= var_thresholdPV)
  case2 <- case2[, -c(n1+n2+1, n1+n2+2)]
  length2 <- dim(case2)[1]

  #p3_0 <- (n - dim(dat_df)[1])/n
  p1_0 <- length1/dim(dat_df)[1]
  p2_0 <- length2/dim(dat_df)[1]
  p4_0 <- 1 - p1_0 - p2_0
  #p1_0 <- (1 - p3_0) * p1_0
  #p2_0 <- (1 - p3_0) * p2_0
  #p4_0 <- (1 - p3_0) * p4_0
  print('initial parameter estimate:')
  #print(c(p1_0, p2_0, p3_0, p4_0, nu_0, var_0, k_0, mu_0))
  #return(c(p1_0, p2_0, p3_0, p4_0, nu_0, var_0, k_0, mu_0))
  print(c(p1_0, p2_0, p4_0, nu_0, var_0, k_0, mu_0))
  return(c(p1_0, p2_0, p4_0, nu_0, var_0, k_0, mu_0))
}




runEM_L_BFGS_B_1 <- function(dat_df, EM_threshold = 0.0001, n1, n2, init_para_est, niter){

  init_para_est_1 <- rep(NA, 7)
  init_para_est_1[1:7] <- init_para_est[1:7]
  #total <- init_para_est[1] + init_para_est[2] + init_para_est[4]
  #prop1 <- init_para_est[1]/total
  #prop2 <- init_para_est[2]/total
  #prop4 <- init_para_est[4]/total
  #init_para_est_1[1] <- prop1
  #init_para_est_1[2] <- prop2
 # init_para_est_1[3] <- prop4

  theta <- matrix(NA, nrow = niter, ncol = 7)
  theta[1,]<-c(0.5,0.25,0.25,6,4,2,0.2)
  theta[2,] <- init_para_est_1

  i <- 2
  while(max(abs(theta[i,]-theta[i-1,]))>0.0001){
    p1<-theta[i,1]
    p2<-theta[i,2]
    p3<-theta[i,3]
    nu<-theta[i,4]
    var<-theta[i,5]
    k<-theta[i,6]
    mu<-theta[i,7]

    epc<-apply(dat_df,1,function(x){
      d1<-f1(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
      d2<-f2(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
      d3<-f3(x[1:n1],x[(n1+1):(n1+n2)],nu,var,k,mu)
      common_term<-log(p1)+d1+log(1+exp(log(p2)-log(p1)+d2-d1)+exp(log(p3)-log(p1)+d3-d1))
      fai1<-log(p1)+d1-common_term
      fai2<-log(p2)+d2-common_term
      fai3<-log(p3)+d3-common_term
      fai1<-exp(fai1)
      fai2<-exp(fai2)
      fai3<-exp(fai3)
      return(c(fai1,fai2,fai3))
    })
    epct<-t(epc)

    keep <- complete.cases(epct)
    epct <- epct[keep,]
    dat_df <- dat_df[keep,]
    theta[i+1,1]<-sum(epct[,1])/dim(dat_df)[1]
    theta[i+1,2]<-sum(epct[,2])/dim(dat_df)[1]
    theta[i+1,3]<-sum(epct[,3])/dim(dat_df)[1]


    loglik<-function(para){
      B<-apply(dat_df,1,function(x){
        a1<-f1(x[1:n1],x[(n1+1):(n1+n2)],para[1],para[2],para[3],para[4])
        b1<-f2(x[1:n1],x[(n1+1):(n1+n2)],para[1],para[2],para[3],para[4])
        c1<-f3(x[1:n1],x[(n1+1):(n1+n2)],para[1],para[2],para[3],para[4])
        return(c(a1,b1,c1))
      })
      return(-sum(diag(epct%*%B)))}

    theta[i+1,4:7]<-optim(c(nu,var,k,mu),loglik, method = "L-BFGS-B")$par
    i <- i + 1
  }
  return(c(i,theta[i,]))
}



runEM_L_BFGS_B <- function(dat_df, EM_threshold = 0.0001, n1, n2, init_para_est, niter){

  EM_para_est <- runEM_L_BFGS_B_1(dat_df = dat_df, EM_threshold = EM_threshold, n1 = n1, n2 = n2, init_para_est = init_para_est, niter = niter)
  i <- EM_para_est[1]
  theta_est <- EM_para_est[-1]

  #PI <- (1-init_para_est[3])*theta_est[1:3]
  #PI<-c(PI[1],PI[2],init_para_est[3],PI[3])
  PI <- theta_est[1:3]
  nu<-theta_est[4]
  var<-theta_est[5]
  k<-theta_est[6]
  mu<-theta_est[7]

  indep_para_est <- c(PI, nu, var, k, mu)
  print('independent EM parameter estimate:')
  print(indep_para_est)
  return(indep_para_est)
}



















runHMM_iters_ThreeCases <- function(emissions, niter, n){
  post <- emissions
  initPI <- c(0.5, 0.25, 0.25)
  Tmat <- t(matrix(c(0.4, 0.3, 0.3, 0.4, 0.3, 0.3, 0.4, 0.3, 0.3), 3, 3))##
  iter <- 0
  oldlogProb = 0
  logProb = -Inf
  c <- rep(NA, n)
  alpha <- matrix(NA, n, 3)
  beta <- matrix(NA, n, 3)
  gamma <- matrix(NA, n, 3)
  bigamma <- array(NA, c(3, 3, n))

  initPI_iters <- list()
  Tmat_iters <- list()
  logProb_iters <- rep(NA, niter)

  while ((iter == 0 | logProb > oldlogProb) & iter <= niter){

    oldlogProb <- logProb
    c[1] <- 0
    for (i in 1:3){
      alpha[1,i] <- initPI[i] * post[1,i]
      c[1] <- c[1] + alpha[1,i]
    }

    #c[1] <- 1/c[1]
    c[1] <- 1/c[1]
    for (i in 1:3){
      alpha[1,i] <- c[1] * alpha[1,i]
    }

    for (t in 2:n){
      c[t] <- 0
      for (i in 1:3){
        alpha[t,i] <- 0
        for (j in 1:3){
          alpha[t,i] <- alpha[t,i] + alpha[t-1,j] * Tmat[j,i]
        }
        alpha[t,i] <- alpha[t,i] * post[t,i]
        c[t] <- c[t] + alpha[t,i]
      }

      #c[t] <- 1/c[t]
      c[t] <- 1/c[t]
      for (i in 1:3){
        alpha[t,i] = c[t] * alpha[t,i]
      }
    }

    for (i in 1:3){
      beta[n,i] <- c[n]
    }

    for (t in (n-1):1){
      for (i in 1:3){
        beta[t,i] <- 0
        for (j in 1:3){
          beta[t,i] <- beta[t,i] + Tmat[i,j] * post[t+1,j] * beta[t+1,j]
        }
        beta[t,i] <- c[t] * beta[t,i]
      }
    }

    for (t in 1:(n-1)){
      denom <- 0
      for (i in 1:3){
        for (j in 1:3){
          denom <- denom + alpha[t,i] * Tmat[i,j] * post[t+1,j] * beta[t+1,j]
        }
      }
      for (i in 1:3){
        gamma[t,i] <- 0
        for (j in 1:3){
          bigamma[i,j,t] <- (alpha[t,i] * Tmat[i,j] * post[t+1,j] * beta[t+1,j])/denom
          gamma[t,i] <- gamma[t,i] + bigamma[i,j,t]
        }
      }
    }

    denom <- 0
    for (i in 1:3){
      denom <- denom + alpha[n,i]
    }
    for (i in 1:3){
      gamma[n,i] <- alpha[n,i]/denom
    }

    for (i in 1:3){
      initPI[i] <- gamma[1,i]
    }
    for (i in 1:3){
      for (j in 1:3){
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



posterior_inference_FDR <- function(gamma, fdr_threshold = 0.1){
  pos_state_prob <- gamma
  pos_prob_null <- apply(pos_state_prob, 1, function(x){return(x[1])})
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



posterior_inference<- function(alpha, beta, Z = NULL, train = TRUE)
{
  pos_state_prob <- alpha * beta
  pos_state <- apply(pos_state_prob, 1, which.max)
  print(table(pos_state))
  if (train == TRUE) {
    print(table(Z))
    acc = sum(pos_state == Z)/length(Z)
    return(list(pos_state, acc))
  }
  else {
    return(pos_state)
  }
}

