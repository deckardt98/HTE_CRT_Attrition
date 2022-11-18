##########################################################################
###   This file reproduces sample size using simple inflation method   ###
##########################################################################
#function calculating number of clusters for equal cluster size
number_of_clusters_equal <- function(alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m, rho_o, rho_c, delta, tyc=0){
  #Input
  #CV: coefficient of variation
  #alpha: type I error, default value is 0.05
  #zeta: type II error, default power is 0.8
  #sigma_o: conditional standard deviation of outcome
  #sigma_c: marginal standard deviation of covariates
  #sigma_in: standard deviation of intervention, default value is 0.5
  #m: cluster size
  #rho_o: o-ICC
  #rho_c: c-ICC
  #delta: pre-specified test size
  #tyc: type of covariate 0=continuous 1=binary
  
  #output:
  #n: number of clusters
  
  #calculate standard deviation of binary covariate
  if (tyc==1){                  
    sum_q1_q2 <- 1/rho_c - 1 
    q1 <- 0.3*sum_q1_q2
    q2 <- sum_q1_q2 - q1
    sigma_c <- sqrt(q1*q2/sum_q1_q2^2)
  }
  n <- (qnorm(1-alpha/2) + qnorm(1-zeta))^2 * sigma_o^2 * (1 - rho_o)*(1 + (m - 1)*rho_o)/ #n: number of clusters
    (m * delta^2 * sigma_i^2 * sigma_c^2 * (1 + (m - 2)*rho_o - (m - 1)*rho_c*rho_o))
  return(n)
}
#function calculating number of clusters by naive method
nc_naive <- function(pi, alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m, rho_o, rho_c, delta, tyc=0){
  #Input
  #pi: mean follow up rate
  #alpha: type I error, default value is 0.05
  #zeta: type II error, default power is 0.8
  #sigma_o: conditional standard deviation of outcome
  #sigma_c: marginal standard deviation of covariates
  #sigma_i: standard deviation of intervention, default value is 0.5
  #m: cluster size
  #rho_o: o-ICC
  #rho_c: c-ICC
  #delta: pre-specified test size
  #tyc: type of covariate 0=continuous 1=binary
  
  #output:
  #n: number of clusters
  
  #calculate standard deviation of binary covariate
  if (tyc==1){
    sum_q1_q2 <- 1/rho_c - 1 
    q1 <- 0.3*sum_q1_q2
    q2 <- sum_q1_q2 - q1
    sigma_c <- sqrt(q1*q2/sum_q1_q2^2)
  }
  n <- number_of_clusters_equal(alpha=alpha, zeta=zeta, sigma_o=sigma_o,  #n: number of clusters
                                sigma_c=sigma_c, sigma_i=sigma_i, m=m, rho_o=rho_o, rho_c=rho_c, delta=delta,tyc=tyc)/pi
  #round up to the nearest even integer
  n <- ceiling(n)
  if (n%%2==1){n<-n+1}
  return(n)
}

m <- c(20,50,100)            #m: cluster size without attrition
pi <- c(0.7, 0.9)            #pi: average follow-up rate
rho_c <- c(0.1, 0.5)         #rho_c: covariate ICC
rho_o <- c(0.01, 0.1)        #rho_o: outcome ICC
delta_c <- c(0.1, 0.25)      #delta_c: effect size for continuous covariate
tau <- c(0.05,0.3,0.6,1)     #tau: correlation of missing indicator

result_naive_c <- c()        #result_naive_c: data set storing sample size of continuous covariate by simple inflation
for (i in 1:length(delta_c)){
  for (l in 1:length(m)){
    for (j in 1:length(rho_c)){
      for (k in 1:length(rho_o)){
        for (v in 1:length(pi)){
          result_naive_c <- c(result_naive_c,nc_naive(pi[v], alpha = 0.05, zeta = 0.2, sigma_o=1, sigma_c=1, sigma_i = 0.5, m[l], rho_o[k], rho_c[j], delta_c[i],tyc=0))
        }
      }
    }
  }
}
result_naive_c

result_naive_b <- c()        #result_naive_b: data set storing sample size of binary covariate by simple inflation
for (i in 1:length(delta_b)){
  for (l in 1:length(m)){
    for (j in 1:length(rho_c)){
      for (k in 1:length(rho_o)){
        for (v in 1:length(pi)){
          result_naive_b <- c(result_naive_b,nc_naive(pi[v], alpha = 0.05, zeta = 0.2, sigma_o=1, sigma_c=1, sigma_i = 0.5, m[l], rho_o[k], rho_c[j], delta_b[i],tyc=1))
        }
      }
    }
  }
}
result_naive_b
