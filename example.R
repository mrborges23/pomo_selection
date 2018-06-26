#Make sure all the packages are installed
#install.packages("MCMCpack")

#Create a directory with the folder source.R and example.txt
dir <- getwd()
setwd(dir)
source("source.R")

#EXAMPLE 1: Simulating allele counts in neutrality

#Simulating som
pi    <- rep(0.25,4)    #base composition
rho   <- rep(0.0001,6)   #exchangeabilities
sigma <- rep(1,4)       #allelic fitness: neutral
par1 <- c(pi,rho,sigma) 

N   <- 10               #number of individuals
S   <- 100000000        #number of sites
sd  <- stationarydistributions(N,pi,rho,sigma)  #allele frequencies
counts <- S*sd                                  #allele counts
data   <- datapreparation(N,counts)             #calculating some statistics

#PRIOR PARAMETERS
a <- rep(1,4)     #pi
b <- rep(0.01,4)  #sigma
c <- rep(0.01,6)  #rho

#obtained with mcmctunning()
vtun <- c(4150000,19,60)
ipar <- c(0.11,0.13,0.47,0.29, 
          0.7001653,0.2251310,0.1709783,0.8133506,0.5477473,0.5095183,
          1.00,1.00,1.00,1.00)
  
#MCMC
sim   <- 100000  #number of iterations
lag   <- 100     #number of utile iterations: sim/lag

smcmc <- mcmcm(N,data,a,b,c,sim,lag,ipar,vtun)  #MCMC run
  
write.table(smcmc,"mcmc_sim1.txt",quote=F,row.names=F,col.names=F)

plotmcmc(smcmc,ipar)

#burnin 25%
smcmc <- smcmc[round(sim*0.25/lag):sim/lag,]

plotmcmc(smcmc)

#plotting mu
plotmutationrates(smcmc,N,F)


#EXAMPLE 2: Simulating allele counts in selection

#Simulating som
pi    <- rep(0.25,4)      #base composition
rho   <- rep(0.0001,6)    #exchangeabilities
sigma <- c(1,1.01,1.01,1) #allelic fitness: neutral
par1 <- c(pi,rho,sigma) 

N   <- 10               #number of individuals
S   <- 100000000        #number of sites
sd  <- stationarydistributions(N,pi,rho,sigma)  #allele frequencies
counts <- S*sd                                  #allele counts
data   <- datapreparation(N,counts)             #calculating some statistics

#PRIOR PARAMETERS
a <- rep(1,4)     #pi
b <- rep(0.01,4)  #sigma
c <- rep(0.01,6)  #rho

#obtained with mcmctunning()
vtun <- c(4150000,19,800)
ipar <- c(0.11,0.13,0.47,0.29, 
          0.7001653,0.2251310,0.1709783,0.8133506,0.5477473,0.5095183,
          1.00,1.1,0.98,0.96)

#MCMC
sim   <- 10000000  #number of iterations
lag   <- 10000     #number of utile iterations: sim/lag

smcmc <- mcmc(N,data,a,b,c,sim,lag,ipar,vtun)  #MCMC run

write.table(smcmc,"mcmc_sim2.txt",quote=F,row.names=F,col.names=F)

plotmcmc(smcmc)

#burnin 25% chain
smcmc <- smcmc[round(sim*0.25/lag):sim/lag,]

plotmcmc(smcmc)

#Plotting mu and sigma
par(mfrow=c(2,1))
plotmutationrates(smcmc,N,F)
plotselectionrates(smcmc,N)


#EXAMPLE 3: Real data: allele counts from pongo pygmaeus

counts <- scan("counts_pongo_pygmaeus.txt")  #reading counts
N <- (length(counts)-4)/6 +1                 #Nstates=4+6*(N-1)

data   <- datapreparation(N,counts)          #calculating some statistics

#PRIOR PARAMETERS
a <- rep(1,4)     #pi
b <- rep(0.01,4)  #sigma
c <- rep(0.01,6)  #rho

#obtained with
#tun <- 10^(-1:ceiling(log(data$s,10)))
#sim   <- 100000
#lag   <- 100
#mcmctun <- mcmctunning(N,data,a,b,c,sim,lag,tun)
#vtun <- mcmctun$vtun
#ipar <- mcmctun$ipar

vtun <- c(7.576923e+06,2.931211e+00,2.516304e+02)
ipar <- c(0.11,0.13,0.47,0.29, 
          0.7001653,0.2251310,0.1709783,0.8133506,0.5477473,0.5095183,
          1.00,1.1,0.98,0.96)

#MCMC
sim   <- 100000  #number of iterations
lag   <- 100     #number of utile iterations: sim/lag

smcmc <- mcmc(N,data,a,b,c,sim,lag,ipar,vtun)  #MCMC run

write.table(smcmc,"mcmc_real.txt",quote=F,row.names=F,col.names=F)

plotmcmc(smcmc)

#burnin 25%
smcmc <- smcmc[round(sim*0.25/lag):(sim/lag),]
plotmcmc(smcmc)

#Plotting mu and sigma
par(mfrow=c(2,1))
plotmutationrates(smcmc,N,F)
plotselectionrates(smcmc,N)
