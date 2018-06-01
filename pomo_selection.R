dir <- getwd()
setwd(dir)

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

#SIMULATED DATA
pi    <- rep(0.25,4)
rho   <- rep(0.001,6)
sigma <- rep(1,4)
par1 <- c(pi,rho,sigma)

N   <- 10
S   <- 10000000
sd  <- stationarydistributions(N,pi,rho,sigma)
counts <- S*sd
data   <- datapreparation(N,counts)

#PRIOR PARAMETERS
a <- rep(1,4)
b <- rep(0.01,4)
c <- rep(0.01,6)
tun <- 10^(-1:ceiling(log(data$s,10)))


#TUNING
mtun <- matrix(NA,ncol=14,nrow=length(tun))

for (t in 1:length(tun)) {
  
  #RANDOM STARTING POINTS
  pi    <- rgamma(4,1,1)
  pi    <- pi/sum(pi)
  rho   <- rexp(6)
  sigma <- c(1,rexp(3))
  
  #MCMC
  sim   <- 10000
  lag   <- 1
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


#RANDOM STARTING POINTS
pi    <- rgamma(4,1,1)
pi    <- pi/sum(pi)
rho   <- rexp(6)
sigma <- c(1,rexp(3))

#MCMC
sim   <- 10000000
lag   <- 1000
usim <- seq(lag,sim,lag)

smcmc <- matrix(NA,ncol=14,nrow=length(usim))
t1 <- Sys.time()
p <- 1
for (i in 1:sim) {
  
  pi    <- mhpi(N,pi,rho,sigma,data,a,tun1)
  rho   <- mhrho(N,pi,rho,sigma,data,c,tun2)
  sigma <- mhsigma(N,pi,rho,sigma,data,b,tun3)
  
  if (i == usim[p]) {
    print(paste("% MCMC: ", round(i*100/sim,1)," | pi: ",round(pi,3)," rho:",round(rho,5)," sigma:",round(rho,3),sep=""))
    smcmc[p,] <- c(pi,rho,sigma)
    p<-p+1
  }
}
t2 <- Sys.time()
write.table(smcmc,"mcmc.txt",quote=F,row.names=F,col.names=F)

print(t1)
print(t2)


pdf("parameters.pdf", 10, 20)
npar <- c("piA","piC","piG","piT","rhoAC","rhoAG","rhoAT","rhoCG","rhoCT","rhoGT","sigmaA","sigmaC","sigmaG","sigmaT")
par(mfrow=c(7,2),mar=c(2,2,2,2))
for (i in 1:14) {
  plot(smcmc[,i],col="grey",xlab="",ylab="",main=npar[i])
  lines(cumsum(smcmc[,i])/(1:length(smcmc[,i])),col="blue")
}
dev.off()


