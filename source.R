#LIST OF FUNCTIONS

#DATA
datapreparation <- function(N,counts) {
  mono  <- counts[1:4]
  poly  <- rep(NA,6)
  kpoly <- rep(NA,6)
  Npoly <- rep(NA,6)
  p<-5
  for (i in 1:6) {
    scounts <- 0 
    kcounts <- 0 
    Ncounts <- 0 
    for (j in 1:(N-1)){
      scounts <- scounts + counts[p]
      kcounts <- kcounts + (j-1)*counts[p]
      Ncounts <- Ncounts + (N-j-1)*counts[p]
      p<- p+1
    }
    poly[i]  <- scounts
    kpoly[i] <- kcounts
    Npoly[i] <- Ncounts
  }
  return(list(s=sum(counts),m=mono,p=poly,kp=kpoly,Np=Npoly))
}

#SCHRINK
newcounts <- function(counts,M) {
  N      <- (length(counts)-4)/6+1
  
  ncounts <- rep(NA,4+6*(M-1))
  ncounts[1:4] <- counts[1:4]
  
  #matrix
  wmatrix <- matrix(NA,ncol=M-1,nrow=N-1)
  m <- 1:(M-1)
  for (i in 1:(N-1)) {
    wmatrix[i,] <- dbinom(m,M,i/N)/sum(dbinom(m,M,i/N))  
  }
  
  for (i in 1:6) {
    nrange <- ((i-1)*(N-1)+5):(i*(N-1)+4)
    mrange <- ((i-1)*(M-1)+5):(i*(M-1)+4)
    ncounts[mrange] <- colSums(counts[nrange]*wmatrix)
  }
  return(ncounts)
}


#NORMALIZATION CONSTANT
normalizationconstants <- function(N,pi,rho,sigma){
  
  i1 <- c(1,1,1,2,2,3)
  i2 <- c(2,3,4,3,4,4)
  
  nc <- sum(pi*(sigma^(N-1)))
  for (i in 1:6) {
    for (j in 1:(N-1)) {
      nc <- nc + pi[i1[i]]*pi[i2[i]]*rho[i]*(sigma[i1[i]]^(j-1))*(sigma[i2[i]]^(N-j-1))*N/(j*(N-j))
    }
  }
  return(nc)
}


#STATIONARY DISTRIBUTION
stationarydistributions <- function(N,pi,rho,sigma){
  
  sd <- rep(NA,4+6*(N-1))
  sd[1:4] <- pi*(sigma^(N-1))
  
  i1 <- c(1,1,1,2,2,3) 
  i2 <- c(2,3,4,3,4,4)
  
  p <- 5
  for (i in 1:6) {
    for (j in 1:(N-1)){
      sd[p] <- pi[i1[i]]*pi[i2[i]]*rho[i]*(sigma[i1[i]]^(j-1))*(sigma[i2[i]]^(N-j-1))*N/(j*(N-j))  
      p <- p+1
    }
  }
  sd <- sd/normalizationconstants(N,pi,rho,sigma)
  return(sd)
}


#LOG CONDITIONAL POSTERIORS

#log posterior pi
logpi <- function(N,pi,rho,sigma,data,a){
  i1 <- c(1,1,1,2,2,3)
  i2 <- c(2,3,4,3,4,4)
  lpost <- -data$s*log(normalizationconstants(N,pi,rho,sigma)) + sum((data$m+a-1)*log(pi)) + sum(data$p*(log(pi[i1])+log(pi[i2])))
  return(lpost)
}

#log posterior rho
logrho <- function(N,pi,rho,sigma,data,c){
  lpost <- -data$s*log(normalizationconstants(N,pi,rho,sigma)) - sum(c*rho) + sum(data$p*log(rho))
  return(lpost)
}

#log posterior sigma
logsigma <- function(N,pi,rho,sigma,data,b) {
  i1 <- c(1,1,1,2,2,3)
  i2 <- c(2,3,4,3,4,4)
  lpost <- -data$s*log(normalizationconstants(N,pi,rho,sigma)) + sum(data$m*(N-1)*log(sigma) - sigma*b) + sum(data$kp*log(sigma[i1])) + sum(data$Np*log(sigma[i2]))
  return(lpost)
}


#PROPOSALS
library(MCMCpack)

#metropolis-hastings pi
mhpi <- function(N,pi0,rho,sigma,data,a,tun){
  lpost0 <- logpi(N,pi0,rho,sigma,data,a)
  tun1   <- tun
  pi1    <- rdirichlet(1,pi0*tun1)
  lpost1 <- logpi(N,pi1,rho,sigma,data,a)
  ratio  <- min(1,exp(lpost1-lpost0)) 
  if (ratio > runif(1)) {
    pi0    <- pi1
    lpost0 <- lpost1
  }
  return(pi0)
}

#metropolis-hastings rho
mhrho <- function(N,pi,rho0,sigma,data,c,tun){
  lpost0 <- logrho(N,pi,rho0,sigma,data,c)
  epson  <- 1/tun 
  mul    <- runif(6,0,1)-0.5
  rho1   <- rho0*exp(epson*mul)
  lpost1 <- logrho(N,pi,rho1,sigma,data,c)
  ratio  <- min(1,exp(lpost1-lpost0)) 
  if (ratio > runif(1)) {
    rho0   <- rho1
    lpost0 <- lpost1
  }
  return(rho0)
}

#metropolis-hastings sigma
mhsigma <- function(N,pi,rho,sigma0,data,b,tun){
  lpost0 <- logsigma(N,pi,rho,sigma0,data,b)
  epson  <- 1/tun
  mul    <- runif(3,0,1)-0.5
  sigma1 <- c(1,sigma0[2:4]*exp(epson*mul))
  lpost1 <- logsigma(N,pi,rho,sigma1,data,b)
  while (is.na(lpost1)) {
    tun    <- runif(3,0,1)-0.5
    sigma1 <- c(1,sigma0[2:4]*exp(epson*tun))
  }
  ratio  <- min(1,exp(lpost1-lpost0)) 
  if (ratio > runif(1)) {
    sigma0 <- sigma1
    lpost0 <- lpost1
  }
  return(sigma0)
}

#TUNING
mcmctunning <- function(N,data,a,b,c,sim,lag,tun) {
  
  mtun <- matrix(NA,ncol=14,nrow=length(tun))
  ipar <- rep(0,14)
  
  for (t in 1:length(tun)) {
  
    #RANDOM STARTING POINTS
    pi    <- rgamma(4,1,1)
    pi    <- pi/sum(pi)
    rho   <- rexp(6)
    sigma <- c(1,rexp(3))
  
    #MCMC
    usim <- seq(lag,sim,lag)
  
    smcmc <- matrix(NA,ncol=14,nrow=length(usim))
  
    p <- 1
    for (i in 1:sim) {
    
      pi    <- mhpi(N,pi,rho,sigma,data,a,tun[t])
      rho   <- mhrho(N,pi,rho,sigma,data,c,tun[t])
      sigma <- mhsigma(N,pi,rho,sigma,data,b,tun[t])
    
      if (i == usim[p]) {
        print(paste("% Tunning: ",round(i*100/sim,1)," (",t,"/",length(tun),")",sep=""))
        smcmc[p,] <- c(pi,rho,sigma)
        p<-p+1
      }
    }
  
    ipar <- smcmc[length(usim),] + ipar
    for (j in 1:14) {
      mtun[t,j] <- (length(unique(smcmc[,j]))-1)/nrow(smcmc)
    }
  }

  y1 <- mtun[sum(mtun[,1] < 0.234),1]
  y2 <- mtun[length(tun)-sum(mtun[,1] > 0.234)+1,1]
  x1 <- tun[sum(mtun[,1] < 0.234)]
  x2 <- tun[length(tun)-sum(mtun[,1] > 0.234)+1]

  tun1 <- (0.234-y1)*(x2-x1)/(y2-y1)+x1

  y1 <- mtun[sum(mtun[,5] < 0.234),5]
  y2 <- mtun[length(tun)-sum(mtun[,5] > 0.234)+1,5]
  x1 <- tun[sum(mtun[,5] < 0.234)]
  x2 <- tun[length(tun)-sum(mtun[,5] > 0.234)+1]

  tun2 <- (0.234-y1)*(x2-x1)/(y2-y1)+x1

  y1 <- mtun[sum(mtun[,12] < 0.234),12]
  y2 <- mtun[length(tun)-sum(mtun[,12] > 0.234)+1,12]
  x1 <- tun[sum(mtun[,12] < 0.234)]
  x2 <- tun[length(tun)-sum(mtun[,12] > 0.234)+1]

  tun3 <- (0.234-y1)*(x2-x1)/(y2-y1)+x1 
  
  return(list(ipar=ipar/length(tun),vtun=c(tun1,tun2,tun3)))
}

#mcmc mutation + selection
mcmc <- function(N,data,a,b,c,sim,lag,ipar,vtun) {
  
  tun1 <- vtun[1]
  tun2 <- vtun[2]
  tun3 <- vtun[3]
  
  #RANDOM STARTING POINTS
  pi    <- ipar[1:4]
  rho   <- ipar[5:10]
  sigma <- ipar[11:14]
  
  usim <- seq(lag,sim,lag)

  smcmc <- matrix(NA,ncol=14,nrow=length(usim))

  p <- 1
  for (i in 1:sim) {
    
    pi    <- mhpi(N,pi,rho,sigma,data,a,tun1)
    rho   <- mhrho(N,pi,rho,sigma,data,c,tun2)
    sigma <- mhsigma(N,pi,rho,sigma,data,b,tun3)
    
    if (i == usim[p]) {
      print(paste("% MCMC: ", round(i*100/sim,1)," | pi: ",paste(round(pi,3),collapse=" ")," rho:",paste(round(rho,5),collapse=" ")," sigma:",paste(round(sigma,3),collapse=" "),sep=""))
      smcmc[p,] <- c(pi,rho,sigma)
      p<-p+1
    }
  }
  
  return(smcmc)
}

#mcmc mutation
mcmcm <- function(N,data,a,b,c,sim,lag,ipar,vtun) {
  
  tun1 <- vtun[1]
  tun2 <- vtun[2]
  tun3 <- vtun[3]
  
  #RANDOM STARTING POINTS
  pi    <- ipar[1:4]
  rho   <- ipar[5:10]
  sigma <- rep(1,4)
  
  usim <- seq(lag,sim,lag)
  
  smcmc <- matrix(NA,ncol=14,nrow=length(usim))
  
  p <- 1
  for (i in 1:sim) {
    
    pi    <- mhpi(N,pi,rho,sigma,data,a,tun1)
    rho   <- mhrho(N,pi,rho,sigma,data,c,tun2)
    sigma <- rep(1,4)
    
    if (i == usim[p]) {
      print(paste("% MCMC: ", round(i*100/sim,1)," | pi: ",paste(round(pi,3),collapse=" ")," rho:",paste(round(rho,5),collapse=" ")," sigma:",paste(round(sigma,3),collapse=" "),sep=""))
      smcmc[p,] <- c(pi,rho,sigma)
      p<-p+1
    }
  }
  
  return(smcmc)
}


plotmcmc <- function(smcmc,ipar) {
  npar <- c("piA","piC","piG","piT","rhoAC","rhoAG","rhoAT","rhoCG","rhoCT","rhoGT","sigmaA","sigmaC","sigmaG","sigmaT")
  par(mfrow=c(7,2),mar=c(2,2,2,2))
  for (i in 1:14) {
    plot(smcmc[,i],col="grey",xlab="",ylab="",main=npar[i])
    lines(cumsum(smcmc[,i])/(1:length(smcmc[,i])),col="blue")
    abline(h=par1[i],col="red")
  }
}

plotmutationrates <- function(smcmc,N,slog)  {
  mmut <- matrix(NA,ncol=12,nrow=nrow(smcmc))
  i1  <- c(1,1,2,2,3,3,4,4,5,5,6,6)
  i2  <- c(2,1,3,1,4,1,3,2,4,2,4,3)
  mut <- c("muAC","muCA","muAG","muGA","muAT","muTA","muCG","muGC","muCT","muTC","muGT","muTG") 
  for (i in 1:12) {
    mmut[,i] <- N*smcmc[,i1[i]+4]*smcmc[,i2[i]]
  }
  
  boxplot(mmut,ylab="scaled mut rates",xaxt='n')
  axis(1,at=1:12,labels=mut)
  if (slog==T) {
    boxplot(log(mmut),ylab="scaled mut rates",xaxt='n')
    axis(1,at=1:12,labels=mut)
  }
} 
  
plotselectionrates <- function(smcmc,N) {
  sig <- c("sigmaC","sigmaG","sigmaT") 
  boxplot(N*(smcmc[,12:14]-1),ylab="scaled sel coefficients",xaxt='n')
  axis(1,at=1:3,labels=sig)
}
