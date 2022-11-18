########################################################################################
###    This file reproduces tables of simulations under MAR for continuous covariate  ###
########################################################################################
require("nlme")
#####################################sample size calculation based on semi-Monte Carlo algorithms#############################################
#Function searching for logistic regression parameters
searlog <- function(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m, B=5000){
  #Input
  #alpha0, alpha1: logistic regression coefficients
  #m: cluster size
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #nc: number of clusters
  #tau_m: correlation of missing indicators
  #B: simulation times 
  
  sigma_i = 0.5 #standard deviation of intervention
  alpha = 0.05 #type I error
  sigma_o = 1  #standard deviation of outcome
  
  #effect size
  beta4 = delta
  cs <- rep(m, nc)                                           #cs: vector of cluster size without attrition
  
  #generate covariates
  sigma_c = 1                                                #sigma_c: standard deviation of covariate
  marpi <- c()                                               #marpi: marginal follow up rate
  for (j in 1:B){
    var_miu <- rho_c*sigma_c^2                               #var_miu: variance of mu
    var_tau <- sigma_c^2 - var_miu                           #var_tau: variance of tau
    miu <- rep(rnorm(nc,0,sqrt(var_miu)) ,cs)                #miu: vector of random effect mu 
    ta <- rnorm(sum(cs),0,sqrt(var_tau))                     #ta: vector of random effect tau
    X <- 1/2 + miu + ta                                      #X: simulated vector of covariate of the entire study
   
    #generate cluster size accounting for attrition  
    #MAR
    PI <- 3.14159265359                                      #PI: constant 
    sigma_b <- sqrt(PI^2/3/(1-tau_m)-PI^2/3)                 #sigma_b: standard deviation of b
    bi <- rnorm(nc,mean = 0, sd = sigma_b)                   #bi: vector of random effect bi of each cluster
    b <- rep(bi, rep(m,nc))                                  #b:  vector of random effect bi of the entire study
    logitp <- c()                                            #logitp: simulated logit of probability of missing indicator
    logitp <- alpha0 + alpha1*X+ b                           
    p <- 1/(1+exp(-logitp))                                  #p: simulated probability of missing indicator
    
    O <- c()                                                 #O: simulated missing indicator 
    for (i in 1:length(p)){
      O[i] <- rbinom(1,1,p[i])
    }
    
    marpi[j] <- mean(O)
    
  }
  return(mean(marpi))
}

#Function generating predicted power by the number of clusters
searc <- function(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m, B=1000, seed=666){
  set.seed(666)
  #Input
  #alpha0, alpha1: logistic regression coefficients
  #m: cluster size
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #nc: number of clusters
  #B: simulation times for eastimating variance beta_4
  #seed: seed for reproducibility
  
  sigma_i = 0.5        #sigma_i: standard deviation of intervention
  alpha = 0.05         #alpha: type I error
  sigma_o = 1          #sigma_o: standard deviation of outcome
  
  #effect size
  beta4 = delta
  cs <- rep(m, nc)                                           #cs: vector of cluster size without attrition
  #generate covariates
  sigma_c = 1                                                #sigma_c: standard deviation of covariate              
  temp <- matrix(0,4,4)                                      #temp: storing each simulated variance covariance matrix
  
  for (j in 1:B) {
    var_miu <- rho_c*sigma_c^2                               #var_miu: variance of mu               
    var_tau <- sigma_c^2 - var_miu                           #var_tau: variance of tau
    miu <- rep(rnorm(nc,0,sqrt(var_miu)) ,cs)                #miu: vector of random effect mu 
    ta <- rnorm(sum(cs),0,sqrt(var_tau))                     #ta: vector of random effect tau
    X <- 1/2 + miu + ta                                      #X: simulated vector of covariate of the entire study
    
    #generate cluster size accounting for attrition  
    #MAR
    PI <- 3.14159265359                                      #PI: constant 
    sigma_b <- sqrt(PI^2/3/(1-tau_m)-PI^2/3)                 #sigma_b: standard deviation of b
    bi <- rnorm(nc,mean = 0, sd = sigma_b)                   #bi: vector of random effect bi of each cluster
    b <- rep(bi, rep(m,nc))                                  #b:  vector of random effect bi of the entire study
    logitp <- c()                                            #logitp: simulated logit of probability of missing indicator
    logitp <- alpha0 + alpha1*X+ b                           
    p <- 1/(1+exp(-logitp))                                  #p: simulated probability of missing indicator
    
    O <- c()                                                 #O: simulated missing indicator 
    for (i in 1:length(p)){
      O[i] <- rbinom(1,1,p[i])
    }
    
    csa <- c()                                               #csa: cluster size after attrition
    for (i in 1:nc){
      csa[i] <- sum(O[((i-1)*m+1):(i*m)])
    }
    
    #generate equal intervention allocation 
    ia <- rep(0, nc)                                         #ia: vector of intervention allocation indicator for each cluster
    rs <- sample(1:nc, size = nc/2)                          #rs: vector of intervention allocation indicator=1 for each cluster
    ia[rs] = 1
    to_ia <- rep(ia, csa)                                    #to_ia: vector of intervention allocation indicator for entire study
    
    #generate new X
    Xnew <- X*O                                              #Xnew: vector of covariate after attrition
    Xnew <- Xnew[which(O!=0)]                                              
    #whole index
    ind <- rep(1:nc, csa)                                    #ind: vector of patients' cluster number
    
    #empirical var4
    var4 <- c()                                              #var4: vector of simulated variance of beta4
    U <- matrix(0, 4, 4)                                     #U: matrix U
    for (i in 1:nc){
      if (csa[i]!=0){
        #Z_i
        col1 <- rep(1, csa[i])                               #col1: column 1 of Z_i
        col2 <- rep(ia[i]-1/2,csa[i])                        #col2: column 2 of Z_i
        if (i > 1) {
          col3 <- Xnew[(sum(csa[1:(i-1)])+1):sum(csa[1:i])]  #col3: column 3 of Z_i
        }else {
          col3 <- Xnew[1:csa[1]]
        }
        col4 <- col2*col3                                    #col4: column 4 of Z_i
        Zi <- as.matrix(cbind(col1,col2,col3,col4))          #Zi: design matrix for each cluster 
        #R_i^{-1}
        Riinv <- (1/(1-rho_o))*diag(1,csa[i])-(rho_o/((1-rho_o)*(1+(csa[i]-1)*rho_o)))*matrix(1, csa[i], csa[i])
        #Z_i^T*R_i^{-1}*Z_i
        Ui <- t(Zi)%*%Riinv%*%Zi
        U <- U+Ui
      }
    }
    U <- U/(nc*sigma_o^2)
    temp <- U + temp
  }
  final <- solve((temp/B))
  var4es <- final[4,4]
  #predicted power
  power <- pnorm(sqrt(nc*delta^2/var4es)-qnorm(1-alpha/2))   #power: predicted power by proposed algorithm
  return(power)
}

#searching for alpha0
alpha0f <- function(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=1000){
  #Input
  #start: lower bound of alpha0
  #end: upper bound of alpha0
  #by: searching precision
  #alpha1: coefficient alpha1
  #m: cluster size without attrition
  #rho_c: covariate ICC
  #rho_o: outcome ICC
  #delta: effect size
  #nc: proposed number of clusters
  #tau: correlation of missing indicator
  #pi: mean follow up rate
  #B: simulation time
  
  for (i in seq(start, end, by=by)){
    mpi <- searlog(i, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)       #mpi: marginal mean follow up rate by proposed number of clusters
    if (mpi>=pi) {
      alpha0 <- i
      break;
    }
  }
  return(alpha0)
}

#function implementing searching 
semimonte <- function(start, nc, pi, m, rho_c, rho_o, delta, tau, alpha1, B=1000,seed=666){
  #Input
  #start: lower bound of alpha0
  #nc: number of clusters from MCAR
  #pi: marginal follow-up rate
  #m: original cluster size 
  #rho_c: covariate ICC
  #rho_o: outcome ICC
  #delta: pre-specified effect size
  #correlation of missing indicator
  #alpha1: fixed logistic coeffcient
  #B: simulation times
  
  end <- 10             #end: reasonable upper bounds of alpha0
  by <- 0.01            #by: searching precision
  
  #searching for alpha0
  alpha0 <- alpha0f(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=B)  #alpha0: coefficient alpha0
  prepower <- searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B,seed=seed) #prepower: predicted power
  if (prepower>0.8){
    while (sprepower>0.8){
      nc <- nc - 2
      alpha0 <- alpha0f(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=B)
      prepower <- searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B,seed=seed)
      if (prepower<=0.8){
        nc <- nc + 2
        break;
      }
    }
  }else {
    while (prepower<0.8){
      nc <- nc + 2
      alpha0 <- alpha0f(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=B)
      prepower <- searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B,seed=seed)
      if (searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B,seed=seed)>=0.8){
        break;
      }
    }
  }
  prepower <- searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B,seed=seed)
  return(c(nc,alpha0,
           prepower,
           searlog(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)))
}

#####################################Simulation for empirical power and empirical type I error rate#############################################
#data generating function
data_gene <- function(m, rho_c, rho_o, delta, nc, alpha0, alpha1, tau_m){
  #Input
  #m: cluster size
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #nc: number of clusters
  #alpha0: searched logistic coefficients
  #alpha1: searched logistic coefficients
  #tau_m: missing indicator correlations
  
  sigma_o = 1 #standard deviation of outcome
  #effect size
  beta1 = 0
  beta2 = 0.25
  beta3 = 0.1
  beta4 = delta
  
  cs <- rep(m, nc)                                    #cs: total cluster size
  
  #generate covariates
  sigma_c = 1                                         #sigma_c: standard deviation of covariate
  
  var_miu <- rho_c*sigma_c^2                 
  var_tau <- sigma_c^2 - var_miu
  miu <- rep(rnorm(nc,0,sqrt(var_miu)) ,cs)
  ta <- rnorm(sum(cs),0,sqrt(var_tau))
  X <- 1/2 + miu + ta
  #generate cluster size accounting for attrition  
  #MAR
  
  PI <- 3.14159265359                                 #PI: constant pi
  sigma_b <- sqrt(PI^2/3/(1-tau_m)-PI^2/3)
  bi <- rnorm(nc,mean = 0, sd = sigma_b)
  b <- rep(bi, rep(m,nc))
  logitp <- c()
  logitp <- alpha0 + alpha1*X+ b
  p <- 1/(1+exp(-logitp))
  
  O <- c()
  for (i in 1:length(p)){
    O[i] <- rbinom(1,1,p[i])
  }
  
  csa <- c()
  for (i in 1:nc){
    csa[i] <- sum(O[((i-1)*m+1):(i*m)])
  }
  #generate equal intervention allocation 
  ia <- rep(0, nc)
  rs <- sample(1:nc, size = nc/2)
  ia[rs] = 1
  to_ia <- rep(ia, csa)
  #generate new X
  Xnew <- X*O
  Xnew <- Xnew[which(O!=0)]
  #whole index
  ind <- rep(1:nc, csa)
  
  #generate outcome
  var_gamma <- rho_o*sigma_o^2
  var_epsilon <- sigma_o^2 - var_gamma
  gamma <- rep(rnorm(nc,0,sqrt(var_gamma)) ,csa)
  epsilon <- rnorm(sum(csa),0,sqrt(var_epsilon))
  
  Y <- beta1 + beta2*to_ia + beta3*Xnew + beta4*Xnew*to_ia + gamma + epsilon
  Y_null <- beta1 + beta2*to_ia + beta3*Xnew + gamma + epsilon
  
  sim_data <- data.frame(ind, to_ia, Xnew, Y, Y_null)
  
  return(sim_data)
}
#function implementing simulation
simula <- function(m, rho_c, rho_o, delta, tau_m, alpha0, alpha1, nc, simt = 1000, seed=666){
  set.seed(seed)      #for reproducibility
  #Input
  #m: cluster size
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #tau_m: missing indicator correlations
  #alpha0: searched logistic coefficients
  #alpha1: searched logistic coefficients
  #nc: the number of clusters
  #simt: the number of simulations per scenario, default = 1000
  
  sigma_i = 0.5  #standard deviation of intervention
  alpha = 0.05   #type I error
  sigma_o = 1    #standard deviation of outcome
  zeta = 0.2     #type II error rate

  pva <- c()
  pvn <- c()
  for (i in 1:simt){
    #generate data
    exit <- FALSE
    while (exit==FALSE){
      df <- data_gene(m, rho_c, rho_o, delta, nc, alpha0 = alpha0, alpha1 = alpha1, tau_m = tau_m)
      lmm1 <- try(lme(Y ~  to_ia*Xnew, random = ~ 1|ind, data = df), silent = T)
      lmm2 <- try(lme(Y_null ~  to_ia*Xnew, random = ~ 1|ind, data = df), silent = T)
      if(class(lmm1)!="try-error"&class(lmm2)!="try-error"){
        exit <- TRUE
    }
    pva[i] <- summary(lmm1)$tTable[4,5] 
    pvn[i] <- summary(lmm2)$tTable[4,5] 
  }
  
  }
  return(c(round(mean(pvn[]<0.05, na.rm=T),3), round(mean(pva[]<0.05, na.rm=T),3)))
}

m <- c(20,50,100)
pi <- c(0.7, 0.9)
rho_c <- c(0.1, 0.5)
rho_o <- c(0.01, 0.1)
delta_c <- c(0.1, 0.25)
tau <- c(0.05,0.3,0.6)

#function searching coefficients and calculating sample size as well as predicted power 
simre <- function(start, nc, pi, m, rho_c, rho_o, delta, tau, alpha1, B=1000,seed=666){
  #Input
  #start: lower bound of alpha0
  #nc: number of clusters from MCAR
  #pi: marginal follow-up rate
  #m: original cluster size 
  #rho_c: covariate ICC
  #rho_o: outcome ICC
  #delta: pre-specified effect size
  #correlation of missing indicator
  #alpha1: fixed logistic coeffcient
  #B: simulation times
  
  set.seed(seed)
  result <- c()
  for (i in 1:length(m)){
    for (k in 1:length(rho_c)){
      for (l in 1:length(rho_o)){
        for (j in 1:length(pi)){
          for (q in 1:length(tau)){
            for (v in 1:length(tau)){
              new <- c(semimonte(start=start,nc=nc,pi[j],m[i],rho_c[k],rho_o[l],delta=delta[q], tau=tau[v],alpha1=0.5,B = B,seed=seed),
                       m[i],rho_c[k],rho_o[l],delta[q],tau[v])
              result <- rbind(result,new)
            }
        }
      }
    }
    }
  }
  return(result)
}
simre(0.05,nc=10,pi=pi,m=m,rho_c=rho_c,rho_o=rho_o,delta=delta_c,tau=tau,alpha1=0.5,B=1000,seed=666)

#simulate again to obtain empirical power 
final <- c()
for (i in 1:nrow(result)){
  final <- rbind(final,simula(m=result[i,5], rho_c=result[i,6], rho_o=result[i,6], delta=result[i,7], tau_m=result[i,8], alpha0=result[i,2], 0.5, nc=result[i,2], 
         simt = 3000, seed=666))
}
