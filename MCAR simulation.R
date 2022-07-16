#function calculating correction factor
correction_factor <- function(CV, m_b, rho_o, rho_c){
  #Input
  #CV: coefficient of variation of sample size m_i
  #m_b: average sample size 
  #rho_o: o-ICC
  #rho_c: c-ICC
  
  #output:
  #CF: correction factor
  CF <- (1 - CV^2 * (m_b*rho_o*(1-rho_o)*(rho_c-rho_o))/((1 + (m_b-2)*rho_o - (m_b - 1)*rho_c*rho_o)*(1 + (m_b - 1)*rho_o)) )^{-1}
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
  n <- correction_factor(CV, m_b, rho_o, rho_c) *
    (qnorm(1-alpha/2) + qnorm(1-zeta))^2 * sigma_o^2 * (1 - rho_o)*(1 + (m_b - 1)*rho_o)/(m_b * delta^2 * sigma_i^2 * sigma_c^2 * (1 + (m_b - 2)*rho_o - (m_b - 1)*rho_c*rho_o))
  return(n)
}

#function calculating number of clusters for equal cluster size
number_of_clusters_equal <- function(alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m_b, rho_o, rho_c, delta){
  #Input
  #CV: coefficient of variation
  #alpha: type I error, default value is 0.05
  #zeta: type II error, default power is 0.8
  #sigma_o: conditional standard deviation of outcome
  #sigma_c: marginal standard deviation of covariates
  #sigma_in: standard deviation of intervention, default value is 0.5
  #m_b: mean cluster size
  #rho_o: o-ICC
  #rho_c: c-ICC
  #delta: pre-specified test size
  n <- (qnorm(1-alpha/2) + qnorm(1-zeta))^2 * sigma_o^2 * (1 - rho_o)*(1 + (m_b - 1)*rho_o)/(m_b * delta^2 * sigma_i^2 * sigma_c^2 * (1 + (m_b - 2)*rho_o - (m_b - 1)*rho_c*rho_o))
  return(n)
}

#function calculating number of clusters by naive method
nc_naive <- function(pi, alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m, rho_o, rho_c, delta){
  m_b <- m*pi
  return(number_of_clusters_equal(alpha, zeta, sigma_o, sigma_c, sigma_i, m_b, rho_o, rho_c, delta)/pi)
}

#tau <- 0.6
#m <- 29
#pi <- 0.61
#sigma_o <- sqrt(0.23)
#sigma_c <- sqrt(0.4)
#m_b <- m*pi
#rho_o <- 0.14
#rho_c <- 0.058
#delta <- 0.3
#CV <- sqrt((1-pi)*(1+tau*(m-1))/pi/m)

#nc_naive(pi, alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m, rho_o, rho_c, delta)
#number_of_clusters_unequal(CV, alpha = 0.05, zeta = 0.2, sigma_o, sigma_c, sigma_i = 0.5, m_b, rho_o, rho_c, delta)


#function calculating the approximated ratio of correct methods vs naive methods given sigma_intervention = 0.5
ratio_correct_naive <- function(pi, m, rho_o, rho_c){
  #Input
  #pi: the rate of outcomes or covariates being observed
  #m: cluster size without attrition
  #rho_o: o-ICC
  #rho_c: c-ICC
  ratio_discard_CF <- (1 - rho_o * (1 + (m-1)*rho_c)/(1 + (m-1)*rho_o)) / (1 - rho_o*(1+(m*pi-1)*rho_c)/(1 + (m*pi - 1)*rho_o))
  return(ratio_discard_CF)
}

#Assess whether naive approach is conservative with MAR data
#Since we use asymptotic CV, CF would be very close to 1 and the ratio is approximately ratio_discard_CF
#Set m = 20, 50 ,100, pi = 0.9, 0.8, 0.6, rho_o = 0.01 0.05 0.1 0.2, rho_c = 0.05 0.25 0.5 0.75
library(ggplot2)
library(gridExtra)
library(tidyverse)

#function generating heatmap dataset
heatdata <- function(m,pi,tau){
  #Input
  #m: cluster size
  #pi: follow-up rate
  #tau: correlation for missing indicator
  
  rho_o <- rep(seq(0.001,0.991,by=0.01),100)
  rho_c <- rep(seq(0.001,0.991,by=0.01),rep(100,100))
  m_b <- m*pi
  CV <- sqrt((1-pi)*(1+tau*(m-1))/pi/m)
  nc_mcar <- number_of_clusters_unequal(CV, alpha = 0.05, zeta = 0.2, 1, 1, sigma_i = 0.5, m_b, rho_o, rho_c, delta=0.2)
  nc_na <- nc_naive(pi, alpha = 0.05, zeta = 0.2, 1, 1, sigma_i = 0.5, m, rho_o, rho_c, delta=0.2)
  z <- round(nc_mcar/nc_na,4)
  df <- as.data.frame(cbind(rho_o,rho_c,z))
  return(df)
}


#plot heat map
phm <- function(pi,tau,m, label){
  m_b <- m*pi
  CV <- sqrt((1-pi)*(1+tau*(m-1))/pi/m)
  x <- seq(0,0.2,by=0.001)
  y <- seq(0,1,by=0.001)
  z <- matrix(data= NA, nrow=201,ncol=1001)
  for (i in 1:201){
    for (j in 1:1001){
      nc_mcar <- number_of_clusters_unequal(CV, alpha = 0.05, zeta = 0.2, 1, 1, sigma_i = 0.5, m_b, x[i], y[j], delta=0.2)
      nc_na <- nc_naive(pi, alpha = 0.05, zeta = 0.2, 1, 1, sigma_i = 0.5, m, x[i], y[j], delta=0.2)
      z[i,j] <- nc_mcar/nc_na
    }
  }
  
  filled.contour(x=seq(0,0.2,by=0.001),
                 y=seq(0,1,by=0.001),
                 z=z,
                 levels = seq(0.5,1.3,length.out = 25),
                 plot.title =   title(paste(label, ": m = ", m, "\u03c0 = ", pi, "\u03c4 = ", tau), xlab="Outcome ICC",
                                      ylab="Covariate ICC"),
                 plot.axes = {
                   axis(1, at=seq(0,0.2,by=0.05));
                   axis(2, at=seq(0,1,by=0.10));
                   contour(x,y,z,nlevels =10,add=TRUE,labcex = 1,col="black")
                 },
                 key.title = {par(cex.main=1.2);title(main="ratio")},
                 key.axes = axis(4, at=seq(0.5,1.3,by=0.1)))
}

m <- c(20,100)
pi <- c(0.9,0.6)
tau <- c(0.05,0.6,1)
a <- phm(pi[1],tau[1],m[1],"(a)")
b <- phm(pi[2],tau[1],m[1],"(b)")
c <- phm(pi[1],tau[1],m[2],"(c)")
d <- phm(pi[2],tau[1],m[2],"(d)")
e <- phm(pi[1],tau[2],m[1],"(e)")
f <- phm(pi[2],tau[2],m[1],"(f)")
g <- phm(pi[1],tau[2],m[2],"(g)")
h <- phm(pi[2],tau[2],m[2],"(h)")
i <- phm(pi[1],tau[3],m[1],"(a)")
j <- phm(pi[2],tau[3],m[1],"(b)")
k <- phm(pi[1],tau[3],m[2],"(c)")
l <- phm(pi[2],tau[3],m[2],"(d)")


#Simulation Codes:

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
  require(geeCRT)
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
  #require nlme package for estimation
  require(nlme)
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
  output <- as.data.frame(cbind(nc,round(mean(pvn[]<0.05, na.rm=T),3),round(tp,3),round(mean(pva[]<0.05, na.rm=T),3)))
  names(output) <- c("number.of.cluster","empirical.typeI.error","predicted.power","empirical.power")
  return(output)
  
}

m <- c(20,50,100)
pi <- c(0.9, 0.7)
rho_c <- c(0.1, 0.5)
rho_o <- c(0.01, 0.1)
delta_c <- c(0.1, 0.25)
delta_b <- c(0.25, 0.45)
tau <- c(0.05,0.3,0.6,1)

#function generating simulation results
simre <- function(m, pi, rho_c, rho_o, delta, tyc = 0, tau, simt = 1000){
  result <- c()
  for (i in 1:length(m)){
  for (j in 1:length(pi)){
    for (k in 1:length(rho_c)){
      for (l in 1:length(rho_o)){
        for (v in 1:length(delta)){
          for (t in 1:length(tau)){
            result <- rbind(result,simula(m[i],pi[j],rho_c[k],rho_o[l],delta[v], tyc = tyc, tau[t],simt = simt))
          }
        }
      }
    }
  }
  }
  return(result)
}

library(geeCRT)
library(nlme)

#function adding label
addla <- function(resul,m,pi,rho_c,rho_o,delta){
  resul$tau <- rep(tau,48)
  resul$delta <- rep(rep(delta,c(4,4)),24)
  resul$rho_o <- rep(rep(rho_o,c(8,8)),12)
  resul$rho_c <- rep(rep(rho_c,c(16,16)),6)
  resul$pi <- rep(rep(pi,c(32,32)),3)
  resul$m <- rep(m,c(64,64,64))
  colnames(resul)[5:10] <- c("tau","effect size","o-ICC",
                             "c-ICC", "pi", "m")
  return(resul)
}

#continuous
mcarcon <- simre(m,pi,rho_c,rho_o,delta_c,tyc=0,tau,simt=3000)
mcarconla <- addla(mcarcon,m,pi,rho_c,rho_o,delta_c)
write.csv(mcarconla,"/Users/deckard/desktop/Fan Li Project/Project1//Continuous MCAR.csv")

#binary
mcarbi <- simre(m,pi,rho_c,rho_o,delta_b,tyc=1,tau,simt=3000)
mcarbila <- addla(mcarbi,m,pi,rho_c,rho_o,delta_b)
write.csv(mcarbila,"/users/deckard/desktop/Fan Li Project/Project//\Binary MCAR.csv")


result_naive_c <- c()
for (i in 1:length(delta_c)){
  for (j in 1:length(rho_c)){
    for (k in 1:length(rho_o)){
      for (l in 1:length(m)){
        for (v in 1:length(pi)){
          result_naive_c <- c(result_naive_c,nc_naive(pi[v], alpha = 0.05, zeta = 0.2, sigma_o=1, sigma_c=1, sigma_i = 0.5, m[l], rho_o[k], rho_c[j], delta_c[i]))
        }
      }
    }
  }
}

result_naive_c <- cbind(result_naive_c,rep(pi,24),rep(rep(m,c(2,2,2)),8),rep(rep(rho_o,c(6,6)),4),rep(rep(rho_c,c(12,12)),2),rep(delta_c,c(24,24)))
colnames(result_naive_c)[2:6] <- c("pi","m","outcome ICC",
                                   "covariate ICC", "effect size")
write.csv(result_naive_c,"/Users/deckard/desktop/Fan Li Project/Project 1/nc_naive_continuous.csv")
