######################################################
###   This file reproduces table 3 of WFHS trial   ###
######################################################
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
  n <- ceiling(n)
  if (n%%2==1){n<-n+1}
  return(n)
}
#function calculating number of clusters for equal cluster size
number_of_clusters_equal <- function(alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m, rho_o, rho_c, delta){
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
  
  #output:
  #n: number of clusters
  n <- (qnorm(1-alpha/2) + qnorm(1-zeta))^2 * sigma_o^2 * (1 - rho_o)*(1 + (m - 1)*rho_o)/
    (m * delta^2 * sigma_i^2 * sigma_c^2 * (1 + (m - 2)*rho_o - (m - 1)*rho_c*rho_o))
  return(n)
}
#function calculating number of clusters by naive method
nc_naive <- function(pi, alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m, rho_o, rho_c, delta){
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
  
  #output:
  #n: number of clusters
  n <- number_of_clusters_equal(alpha=alpha, zeta=zeta, sigma_o=sigma_o, sigma_c=sigma_c, sigma_i=sigma_i, m=m, rho_o=rho_o, rho_c=rho_c, delta=delta)/pi
  n <- ceiling(n)
  if (n%%2==1){n<-n+1}
  return(n)
}

n_mcar <- c()                     #n_mcar: storing sample size under MCAR
n_simple <- c()                   #n_simple: storing sample size using simple inflation method
pi <- c(0.935,0.87,0.61)          #pi: mean follow up rate 
tau <- c(0.05,0.3,0.6)            #tau: correlation of missing indicator
delta <- c(0.05,0.2,0.3)          #delta: standardized effect size
m <- 29                           #m: average cluster size              
sigma_i <- 0.5                    #sigma_i: standard error of treatment indicator
rho_o <- 0.14                     #rho_o: outcome icc
rho_c <- 0.058                    #rho_c: covariate icc
sigma_o <- sqrt(0.23)             #sigma_o: standard deviation of outcome 
sigma_c <- sqrt(0.4)              #sigma_c: standard deviation of covariate
#calculate sample size under MCAR and simple inflation adjustment
for (i in 1:length(pi)){
  for (k in 1:length(delta)){
     for (j in 1:length(tau)){
      m_b <- m*pi[i]                        #m_b: average cluster size adjusting for attrition
      CV <- sqrt((1-pi[i])*tau[j]/pi[i])    #CV: coefficient of variation of cluster size
      #sample size under MCAR
      n_mcar <- c(n_mcar,number_of_clusters_unequal(CV, alpha = 0.05, zeta = 0.2, sigma_o=sigma_o,
                                                    sigma_c=sigma_c, sigma_i = sigma_i, m_b=m_b, rho_o=rho_o, rho_c=rho_c, delta=delta[k]))
      #sample size using simple inflation adjustment
      n_simple <- c(n_simple,nc_naive(pi[i],alpha = 0.05, zeta = 0.2, sigma_o=sigma_o, 
                                      sigma_c=sigma_c, sigma_i = sigma_i, m=m, rho_o=rho_o, rho_c=rho_c, delta=delta[k]))
    }
  }
}
n_mcar
n_simple
