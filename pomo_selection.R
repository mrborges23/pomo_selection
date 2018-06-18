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

#TUNNING MCMC
tun <- 10^(-1:ceiling(log(data$s,10)))
sim   <- 10000
lag   <- 1

#function here
vtun <- mcmctunning(data,a,b,c,tun,sim,lag)


#MCMC
sim   <- 10000000
lag   <- 1000

mcmc <- function(data,a,b,c,sim,lag,vtun) {

  tun1 <- vtun[1]
  tun2 <- vtun[1]
  tun3 <- vtun[1]
  
  #RANDOM STARTING POINTS
  pi    <- rgamma(4,1,1)
  pi    <- pi/sum(pi)
  rho   <- rexp(6)
  sigma <- c(1,rexp(3))

  smcmc <- matrix(NA,ncol=14,nrow=length(usim))

  usim <- seq(lag,sim,lag)

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

smcmc <- mcmc(data,a,b,c,sim,lag,vtun)
  
write.table(smcmc,"mcmc.txt",quote=F,row.names=F,col.names=F)



pdf("parameters.pdf", 10, 20)
npar <- c("piA","piC","piG","piT","rhoAC","rhoAG","rhoAT","rhoCG","rhoCT","rhoGT","sigmaA","sigmaC","sigmaG","sigmaT")
par(mfrow=c(7,2),mar=c(2,2,2,2))
for (i in 1:14) {
  plot(smcmc[,i],col="grey",xlab="",ylab="",main=npar[i])
  lines(cumsum(smcmc[,i])/(1:length(smcmc[,i])),col="blue")
}
dev.off()
