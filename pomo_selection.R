dir <- getwd()
setwd(dir)
source("source.R")

#SIMULATED DATA
pi    <- rep(0.25,4)
rho   <- rep(0.001,6)
sigma <- rep(1,4)
par1 <- c(pi,rho,sigma)

N   <- 5
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

vtun <- mcmctunning(N,data,a,b,c,tun,sim,lag)

#MCMC
sim   <- 10000000
lag   <- 1000

smcmc <- mcmc(N,data,a,b,c,sim,lag,vtun)
  
#write.table(smcmc,"mcmc.txt",quote=F,row.names=F,col.names=F)

npar <- c("piA","piC","piG","piT","rhoAC","rhoAG","rhoAT","rhoCG","rhoCT","rhoGT","sigmaA","sigmaC","sigmaG","sigmaT")
par(mfrow=c(7,2),mar=c(2,2,2,2))
for (i in 1:14) {
  plot(smcmc[,i],col="grey",xlab="",ylab="",main=npar[i])
  lines(cumsum(smcmc[,i])/(1:length(smcmc[,i])),col="blue")
}

