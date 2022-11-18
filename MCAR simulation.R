#################################################################
###   This file reproduces tables of simulations under MCAR   ###
#################################################################
require("geeCRT")
require("nlme")
#function calculating correction factor
correction_factor <- function(CV, m_b, rho_o, rho_c){
  #Input
  #CV: coefficient of variation of sample size m_i
  #m_b: average sample size 
  #rho_o: o-ICC
  #rho_c: c-ICC
  
  #output:
  #CF: correction factor
  CF <- (1 - CV^2 * (m_b*rho_o*(1-rho_o)*(rho_c-rho_o))/((1 + (m_b-2)*rho_o - (m_b - 1)*rho_c*rho_o)*
                                                           (1 + (m_b - 1)*rho_o)^2) )^{-1}
  return(CF)
}

#function calculating number of clusters for unequal cluster size
number_of_clusters_unequal <- function(CV, alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m_b, rho_o, rho_c, delta){
  #Input
  #CV: coefficient of variation
  #alpha: type I error, default value is 0.05
  #zeta: type II error, default power is 0.8
  #sigma_o: conditional standard deviation of outcome
  #sigma_c: marginal standard deviation of covariates
  #sigma_i: standard deviation of intervention, default value is 0.5
  #m_b: mean cluster size
  #rho_o: o-ICC
  #rho_c: c-ICC
  #delta: pre-specified test size
  
  #output:
  #n: number of clusters
  n <- correction_factor(CV=CV, m_b=m_b, rho_o=rho_o, rho_c=rho_c) *
    (qnorm(1-alpha/2) + qnorm(1-zeta))^2 * sigma_o^2 * (1 - rho_o)*(1 + (m_b - 1)*rho_o)/(m_b * delta^2 * sigma_i^2 * sigma_c^2 * (1 + (m_b - 2)*rho_o - (m_b - 1)*rho_c*rho_o))
  return(n)
}

############################Simulation Codes#####################################################################

#Power calculation function
power_cal <- function(alpha = 0.05, sigma_o, rho_o, m, pi, delta, sigma_c, sigma_i = 0.5, rho_c, n, tau, CV){
  
  #Input
  #alpha: type I error, default = 0.05
  #sigma_o: standard deviation of outcome
  #rho_o: o-ICC
  #m: cluster size
  #pi: observing rate
  #delta: pre-specified effect size
  #sigma_c: standard deviation of covariate
  #sigma_i: standard deviation of intervention
  #rho_c: c-ICC
  #n: the number of cluster
  #tau: missing parameter
  
  m_b <- m*pi
  CF <- correction_factor(CV, m_b, rho_o, rho_c)
  ss <- n*m_b*delta^2*sigma_i^2*sigma_c^2*(1 + (m_b - 2)*rho_o - (m_b - 1)*rho_c*rho_o) / ( CF * sigma_o^2*(1-rho_o)*(1+(m_b-1)*rho_o) )
  power <- pnorm(sqrt(ss) - qnorm(1-alpha/2))
  return(power)
}

#data generation function
data_gene <- function(m, pi, tau, rho_c, rho_o, delta, nc, tyc = 0){

  #Input
  #m: cluster size
  #pi: mean observing rate
  #tau: correlation for the missingness
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #nc: number of clusters
  #tyc: tyc of covariates, 0 refers to continuous and 1 refers to binary
  sigma_i = 0.5 #standard deviation of intervention
  alpha = 0.05 #type I error
  zeta = 0.2 #type II error
  sigma_o = 1 #standard deviation of outcome
  #effect size
  beta1 = 0
  beta2 = 0.25
  beta3 = 0.1
  beta4 = delta
  
  
  #generate cluster size accounting for attrition  
  #MCAR
  cs <- c()
  for (i in 1:nc){
    if (tau<1){
    mu <- rep(pi, m)
    sigma <- matrix(tau, nrow = m, ncol = m)
    for (i in 1:m){
      sigma[i,i] <- 1
    }
    cs <- c(sum(simbinCLF(mu,sigma,n =1)),cs)
    } else if (tau==1){
      cs[i] <- rbinom(1,1,pi)*m
    }
  }

  #generate equal intervention allocation 
  ia <- rep(0, nc)
  rs <- sample(1:nc, size = nc/2)
  ia[rs] = 1
  to_ia <- rep(ia, cs)
  
  #generate covariates
  if (tyc == 0) {
    sigma_c = 1
  } 
  if (tyc ==1) {
    sum_q1_q2 <- 1/rho_c - 1 
    q1 <- 0.3*sum_q1_q2
    q2 <- sum_q1_q2 - q1
  }
  
  if (tyc == 0) {
    var_miu <- rho_c*sigma_c^2
    var_tau <- sigma_c^2 - var_miu
    miu <- rep(rnorm(nc,0,sqrt(var_miu)) ,cs)
    ta <- rnorm(sum(cs),0,sqrt(var_tau))
    X <- 1/2 + miu + ta
  } 
  if (tyc == 1){
    pi_i <- rbeta(nc,q1,q2)
    X <- c()
    for (i in 1:nc){
      X <- c(X,rbinom(cs[i],1,pi_i[i]))
    }
  }
  #generate outcome
  var_gamma <- rho_o*sigma_o^2
  var_epsilon <- sigma_o^2 - var_gamma
  gamma <- rep(rnorm(nc,0,sqrt(var_gamma)) ,cs)
  epsilon <- rnorm(sum(cs),0,sqrt(var_epsilon))
  #alternative hypothesis
  Y_alt <- beta1 + beta2*to_ia + beta3*X + beta4*X*to_ia + gamma + epsilon
  #null hypothesis
  Y_null <-  beta1 + beta2*to_ia + beta3*X + gamma + epsilon
  #whole index
  ind <- rep(1:nc, cs)
  
  sim_data <- data.frame(ind, to_ia, X, Y_alt, Y_null)
  
  return(sim_data)
}

#function implementing simulation
simula <- function(m, pi, rho_c, rho_o, delta,  tyc = 0, tau, simt = 1000){
  #Input
  #m: cluster size
  #pi: mean observing rate
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #tyc: tyc of covariates, 0 refers to continuous and 1 refers to binary
  #simt: the number of simulations per scenario, default = 1000
  sigma_i = 0.5 #standard deviation of intervention
  alpha = 0.05 #type I error
  sigma_o = 1 #standard deviation of outcome
  zeta = 0.2 #type II error rate
  
  
  if (tyc == 0){sigma_c <- 1}
  if (tyc == 1){
    sum_q1_q2 <- 1/rho_c - 1 
    q1 <- 0.3*sum_q1_q2
    q2 <- sum_q1_q2 - q1
    sigma_c <- sqrt(q1*q2/sum_q1_q2^2)
  }
  
  #calculate the sample size via our formula
  m_b <- m*pi
  CV <- sqrt((1-pi)*(1+tau*(m-1))/pi/m)
  nc <- ceiling(number_of_clusters_unequal(CV, alpha = alpha, zeta = zeta, sigma_o, sigma_c, sigma_i = sigma_i, m_b, rho_o, rho_c, delta))
  #round up to the nearest even integer
  if (nc %% 2 ==1){nc <- nc + 1}
  #store p-values for null distribution
  pvn <- array(NA,dim=simt)
  #store p-values for alternative distribution
  pva <- array(NA,dim=simt)
  for (i in 1:simt){
    #generate data
    
    #indicator for singular-fitting
    exit <- FALSE
    while (exit==FALSE){
      df <- data_gene(m, pi, tau, rho_c, rho_o, delta, nc, tyc = tyc)
      lmm1 <- try(lme(Y_null ~  to_ia*X, random = ~ 1|ind, data = df))
      lmm2 <- try(lme(Y_alt ~  to_ia*X, random = ~ 1|ind, data = df))
      if(class(lmm1)!="try-error"&class(lmm2)!="try-error"){
        exit <- TRUE
      }
    }
    pvn[i] <- summary(lmm1)$tTable[4,5] 
    pva[i] <- summary(lmm2)$tTable[4,5] 
  }
  
  
  #calculate true power
  tp <- power_cal(alpha = alpha, sigma_o, rho_o, m, pi, delta, sigma_c, sigma_i = sigma_i, rho_c, nc, tau, CV)
  output <- as.data.frame(cbind(nc,round(mean(pvn[]<0.05, na.rm=T),4),round(tp,3),round(mean(pva[]<0.05, na.rm=T),3)))
  names(output) <- c("number.of.cluster","empirical.typeI.error","predicted.power","empirical.power")
  return(output)
  
}

m <- c(20,50,100)
pi <- c(0.7, 0.9)
rho_c <- c(0.1, 0.5)
rho_o <- c(0.01, 0.1)
delta_c <- c(0.1, 0.25)
delta_b <- c(0.25, 0.45)
tau <- c(0.05,0.3,0.6,1)

#function generating simulation results
simre <- function(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 1000,seed=666){
  #Input
  #m: cluster size without attrition
  #pi: mean follow up rate
  #rho_c: covariate ICC
  #rho_o: outcome ICC
  #delta: effect size
  #tyc: type of covariates 0=continuous 1=binary
  #tau: correlation of missing indicator
  #simt: simulation time
  #seed: seed for reproducibility

  set.seed(seed)
  result <- c()
      for (i in 1:length(m)){
         for (k in 1:length(rho_c)){
           for (l in 1:length(rho_o)){
               for (j in 1:length(pi)){
            result <- rbind(result,simula(m[i],pi[j],rho_c[k],rho_o[l],delta=delta, tyc = tyc, tau=tau,simt = simt))
          }
        }
      }
  }
  return(result)
}
#continous covariate delta=0.25 tau=0.05
delta <- delta_c[1]
tau <- tau[1]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.45 tau=0.05
delta <- delta_c[2]
tau <- tau[1]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.25 tau=0.3
delta <- delta_c[1]
tau <- tau[2]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.45 tau=0.3
delta <- delta_c[2]
tau <- tau[2]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.25 tau=0.6
delta <- delta_c[1]
tau <- tau[3]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.45 tau=0.6
delta <- delta_c[2]
tau <- tau[3]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.25 tau=1
delta <- delta_c[1]
tau <- tau[4]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=666)
#continous covariate delta=0.45 tau=1
delta <- delta_c[2]
tau <- tau[4]
simre(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 3000,seed=988)

#Binary covariate delta=0.25 tau=0.05
delta <- delta_b[1]
tau <- tau[1]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.45 tau=0.05
delta <- delta_b[2]
tau <- tau[1]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.25 tau=0.3
delta <- delta_b[1]
tau <- tau[2]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.45 tau=0.3
delta <- delta_b[2]
tau <- tau[2]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.25 tau=0.6
delta <- delta_b[1]
tau <- tau[3]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.45 tau=0.6
delta <- delta_b[2]
tau <- tau[3]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.25 tau=1
delta <- delta_b[1]
tau <- tau[4]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=666)
#Binary covariate delta=0.45 tau=1
delta <- delta_b[2]
tau <- tau[4]
simre(m, pi, rho_c, rho_o, delta, tyc = 1, tau, simt = 3000,seed=988)

