dir <- getwd()
setwd(dir)

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
sim   <- 10000
lag   <- 1

#TUNING
mcmctunninf <- function(tun,sim,lag) {
  
 mtun <- matrix(NA,ncol=14,nrow=length(tun))

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
  
  return(c(tun1,tun2,tun3))
}


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
    print(paste("% MCMC: ", round(i*100/sim,1)," | pi: ",paste(round(pi,3),collapse=" ")," rho:",paste(round(rho,5),collapse=" ")," sigma:",paste(round(sigma,3),collapse=" "),sep=""))
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
