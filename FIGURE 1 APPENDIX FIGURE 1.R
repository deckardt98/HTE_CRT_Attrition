###############################################################
###   This file reproduces Figure 1 and Appendix Figure 1   ###
###############################################################
library(ggplot2)
library(gridExtra)
library(tidyverse)
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
  return(n)
}



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
  nc_mcar <- number_of_clusters_unequal(CV, alpha = 0.05, zeta = 0.2, sigma_o=1, sigma_c=1, sigma_i = 0.5, m_b, rho_o, rho_c, delta=0.2)
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
                 levels = seq(0.7,1.15,length.out = 25),
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