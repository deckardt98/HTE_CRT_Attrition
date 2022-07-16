#MAR
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
  sigma_o = 1 #standard deviation of outcome
  #effect size
  
  beta4 = delta
  cs <- rep(m, nc)
  
  
  #generate covariates
  sum_q1_q2 <- 1/rho_c - 1 
  q1 <- 0.3*sum_q1_q2
  q2 <- sum_q1_q2 - q1
  
  marpi <- c()
  for (j in 1:B){
    pi_i <- rbeta(nc,q1,q2)
    X <- c()
    for (i in 1:nc){
      X <- c(X,rbinom(cs[i],1,pi_i[i]))
    }
    
    #generate cluster size accounting for attrition  
    #MAR
    
    PI <- 3.14159265359
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
    
    marpi[j] <- mean(O)
    
  }
  return(mean(marpi))
}

#Function generating predicted power by the number of clusters
searc <- function(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m, B=1000){
  #Input
  #alpha0, alpha1: logistic regression coefficients
  #m: cluster size
  #rho_c: c-ICC
  #rho_o: i-ICC
  #delta: pre-specified effect size
  #nc: number of clusters
  #B: simulation times for eastimating variance beta_4
  sigma_i = 0.5 #standard deviation of intervention
  alpha = 0.05 #type I error
  sigma_o = 1 #standard deviation of outcome
  #effect size
  
  beta4 = delta
  cs <- rep(m, nc)
  
  
  #generate covariates
  
  sum_q1_q2 <- 1/rho_c - 1 
  q1 <- 0.3*sum_q1_q2
  q2 <- sum_q1_q2 - q1
  
  temp <- matrix(0,4,4)
  
  for (j in 1:B) {
    
    pi_i <- rbeta(nc,q1,q2)
    X <- c()
    for (i in 1:nc){
      X <- c(X,rbinom(cs[i],1,pi_i[i]))
    }
    #generate cluster size accounting for attrition  
    #MAR
    
    PI <- 3.14159265359
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

    #empirical var4
    var4 <- c()
    U <- matrix(0, 4, 4)
    for (i in 1:nc){
      if (csa[i]!=0){
        #Z_i
        col1 <- rep(1, csa[i])
        col2 <- rep(ia[i]-1/2,csa[i])
        if (i > 1) {
          col3 <- Xnew[(sum(csa[1:(i-1)])+1):sum(csa[1:i])]
        }else {
          col3 <- Xnew[1:csa[1]]
        }
        col4 <- col2*col3
        Zi <- as.matrix(cbind(col1,col2,col3,col4))
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
  
  power <- pnorm(sqrt(nc*delta^2/var4es)-qnorm(1-alpha/2))
  
  return(power)
}



#searching for alpha0
alpha0f <- function(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=1000){
  for (i in seq(start, end, by=by)){
    mpi <- searlog(i, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)
    if (mpi>=pi) {
      alpha0 <- i
      break;
    }
  }
  
  return(alpha0)
}

#function implementing searching 
semimonte <- function(start, nc, pi, m, rho_c, rho_o, delta, tau, alpha1, B=1000){
  #Input
  #nc: number of clusters from MCAR
  #pi: marginal follow-up rate
  #m: original cluster size 
  #rho_c: covariate ICC
  #rho_o: outcome ICC
  #delta: pre-specified effect size
  #correlation of missing indicator
  #alpha1: fixed logistic coeffcient
  #B: simulation times
  
  end <- 10
  by <- 0.01
  
  
  #searching for alpha0
  alpha0 <- alpha0f(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=B)
  
  if (searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)>0.8){
    while (searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)>0.8){
      nc <- nc - 2
      alpha0 <- alpha0f(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=B)
      if (searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)<=0.8){
        nc <- nc + 2
        break;
      }
    }
  }else {
    while (searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)<0.8){
      nc <- nc + 2
      alpha0 <- alpha0f(start, end, by, alpha1, m, rho_c, rho_o, delta, nc, tau, pi, B=B)
      if (searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)>=0.8){
        break;
      }
    }
  }
  
  return(c(nc,alpha0,
           searc(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B),
           searlog(alpha0, alpha1, m, rho_c, rho_o, delta, nc, tau_m=tau, B=B)))
}

#function calculating empirical var4 and estimated var4
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
  
  cs <- rep(m, nc)
  
  #generate covariates
  sum_q1_q2 <- 1/rho_c - 1 
  q1 <- 0.3*sum_q1_q2
  q2 <- sum_q1_q2 - q1
  
  pi_i <- rbeta(nc,q1,q2)
  X <- c()
  for (i in 1:nc){
    X <- c(X,rbinom(cs[i],1,pi_i[i]))
  }
  #generate cluster size accounting for attrition  
  #MAR
  
  PI <- 3.14159265359
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


#b <- read.csv("Binary MCAR.csv")
#b <- b[b$tau!=1,]
#b <- b[,-c(1)]

#2.1
#2.48
#3.4

#0.7
#0.9
#1.3

#nc.mar.b <- c()
#predicted.power.b <- c()
#alpha0.b <- c()

#for (i in 1:24){
#  if (i%%3==1){
#    semob <- semimonte(2.1, b[i,1], b[i,9], b[i,10], b[i,8], b[i,7], b[i,6], b[i,5], 0.5, B=1000)
#    nc.mar.b[i] <- semob[1]
#    alpha0.b[i] <- semob[2]
#    predicted.power.b[i] <- semob[3]
#  } else if (i%%3==2){
#    semob <- semimonte(2.48, b[i,1], b[i,9], b[i,10], b[i,8], b[i,7], b[i,6], b[i,5], 0.5, B=1000)
#    nc.mar.b[i] <- semob[1]
#    alpha0.b[i] <- semob[2]
#    predicted.power.b[i] <- semob[3]
#  } else if (i%%3==0){
#    semob <- semimonte(3.4, b[i,1], b[i,9], b[i,10], b[i,8], b[i,7], b[i,6], b[i,5], 0.5, B=1000)
#    nc.mar.b[i] <- semob[1]
#    alpha0.b[i] <- semob[2]
#    predicted.power.b[i] <- semob[3]
#  }
#}

#marb1 <- cbind(nc.mar.b,predicted.power.b,alpha0.b)


#for (i in 25:48){
#  if (i%%3==1){
#    semob <- semimonte(0.7, b[i,1], b[i,9], b[i,10], b[i,8], b[i,7], b[i,6], b[i,5], 0.5, B=1000)
#    nc.mar.b[i] <- semob[1]
#    alpha0.b[i] <- semob[2]
#    predicted.power.b[i] <- semob[3]
#  } else if (i%%3==2){
#    semob <- semimonte(0.9, b[i,1], b[i,9], b[i,10], b[i,8], b[i,7], b[i,6], b[i,5], 0.5, B=1000)
#    nc.mar.b[i] <- semob[1]
#    alpha0.b[i] <- semob[2]
#   predicted.power.b[i] <- semob[3]
# } else if (i%%3==0){
#    semob <- semimonte(1.3, b[i,1], b[i,9], b[i,10], b[i,8], b[i,7], b[i,6], b[i,5], 0.5, B=1000)
#    nc.mar.b[i] <- semob[1]
#    alpha0.b[i] <- semob[2]
#    predicted.power.b[i] <- semob[3]
#  }
#}

#marb2 <- cbind(nc.mar.b,predicted.power.b,alpha0.b)
#marb2 <- read.csv("marb2.csv")
#marb2 <- marb2[25:48,-c(1)]
#marb3 <- read.csv("marb3.csv")
#marb3 <- marb3[49:72,-c(1)]
#marb4 <- read.csv("marb4.csv")
#marb4 <- marb4[73:96,-c(1)]
#marb5 <- read.csv("marb5.csv")
#marb5 <- marb5[97:120,-c(1)]
#arb6 <- read.csv("marb6.csv")
#marb6 <- marb6[121:144,-c(1)]
#marb <- rbind(marb1,marb2,marb3,marb4,marb5,marb6)
#b <- cbind(b,marb)
#rownames(b) <- c(1:nrow(b))
#b <- b[,-c(2,3,4)]
#write.csv(b, "~/Desktop/Fan Li Project/Project 1/Paper Material/Simulation results/MAR Binary Final.csv")

#function implementing simulation
simula <- function(m, rho_c, rho_o, delta, tau_m, alpha0, alpha1, nc, simt = 1000){
  #require nlme package for estimation
  require(nlme)
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
  sigma_i = 0.5 #standard deviation of intervention
  alpha = 0.05 #type I error
  sigma_o = 1 #standard deviation of outcome
  zeta = 0.2 #type II error rate
  
  
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
    }
    pva[i] <- summary(lmm1)$tTable[4,5] 
    pvn[i] <- summary(lmm2)$tTable[4,5] 
  }
  
  return(c(round(mean(pvn[]<0.05, na.rm=T),3), round(mean(pva[]<0.05, na.rm=T),3)))
  
}

marfi <- read.csv("MAR Binary Final.csv")
marfi <- marfi[,-c(1)]
empirical.power.mar1 <- c()
for (i in 1:24){
  empirical.power.mar1[i] <- simula(marfi[i,7], marfi[i,5], marfi[i,4], 
                                    marfi[i,3], marfi[i,2], marfi[i,10], 
                                    0.5, marfi[i,8], simt = 3000)
}

