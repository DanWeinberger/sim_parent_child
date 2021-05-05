#In this simulation, we generate the colonization status each day for a child and their parent
#We specify the community acuistion rate and the probability of transmitting
#to a HH member if you are colonized
#We track the household colonization status in 4 states: neither colonized; child colonized, parent colonized, both colonized
#We also track individual colonization status as uncolonized, colonized from community acquisition or colonized from HH contact


source('sirFunc.R')
source('LikelihoodFunc.R')
library(dplyr)
library(pbapply)
library(rjags)
library(HDInterval)
set.seed(123)

#Generates simulated data from a Markov model for 365 time points
# At each time point, can be 1: neither colonized; 2: child colonized; 3: adult colonized; 4: both colonized
#Each column is a household

logit <- function(x){
  log(x/(1-x))
}

ilogit <-function(x){
  exp(x)/(1+exp(x))
}

#Generate synthetic cross sectional data and tabulate the number of households
#uninfected, with kid infected, adult infected, or both
#set beta to be small, representing low risk from individual ST
TrueData <-  replicate(1000, sirHH( time=400, 
                                  logit.beta= logit(c(1/1200, 1/7000  ))  , #community infection rate for kids and adults #Make this Very small to represent ST-specfici risk
                                  durR=c(60,60), #waning of protection from subsequent infection
                                  durInf=c(60,21), #1/duration of colonization for kids and adults
                                  burn.days=1
                                          ), 
                       simplify='array')
TrueData80 <- t(TrueData[,150,] )
#table(TrueData80[,1],TrueData80[,2])
Y <- as.vector(table(factor(TrueData80[,1],levels=c('0','1')), factor(TrueData80[,2],levels=c('0','1'))))

#Check if TrueData is at steady state
prev <- apply(TrueData,c(1,2), mean )

matplot(t(prev))
plot(prev[1,],ylim=c(0,0.1))
plot(prev[2,],ylim=c(0,0.01))
plot(prev[1,]*prev[2,],ylim=c(0,0.001)) #Co-colonization


#Set constraints for likelihood to prevent extreme values (causes problems for optimization)
 lower.prob <- logit(0.0001)
 upper.prob <- logit(0.9999)
 parms.l <- rep(lower.prob,length(parms))
 parms.u <- rep(upper.prob,length(parms))

parms <-  c(0,0)
#parms <- c(logit(c(1/150, 1/(365*2))), logit(c(0.1*1/120, 0.5*1/120)) ) #feed in correct parms
ptm <- proc.time()
  mod1 <- optim(parms,LL.hh, lower=parms.l, obs=Y, upper=parms.u, method='L-BFGS-B' )
  #mod1 <- nlm(LL.hh, parms, obs=Y )

proc.time() - ptm

parm.est <- mod1$estimate
ilogit(parm.est)


parms <-  rep(0,2)

ptm <- proc.time()
mod1 <- nlm(LL.hh, parms, obs=Y)
proc.time() - ptm

parm.est <- mod1$estimate
ilogit(parm.est)

#Josh's suggestion
biv_mod <- "model{
  for(i in 1:N){
    COL1[i] ~ dbern (P1[i]) 
    COL2[i] ~ dbern (P2[i])
    
    
    logit(P1[i]) <- alpha[1] +phi[i]
    logit(P2[i]) <- alpha[2] + phi[i]

    phi[i] ~dnorm(0, tau.disp)
   
  }
  # Prior distribution of model parameters
  alpha[1] ~ dnorm(0,0.001)
  alpha[2] ~ dnorm(0,0.001)
  tau.disp ~dgamma(0.001, 0.001)

}"


##############################################################
#Model Fitting
##############################################################
inits1=list(".RNG.seed"=c(123), ".RNG.name"='base::Wichmann-Hill')
inits2=list(".RNG.seed"=c(456), ".RNG.name"='base::Wichmann-Hill')
inits3=list(".RNG.seed"=c(789), ".RNG.name"='base::Wichmann-Hill')


##############################################
#Model Organization
##############################################
model_spec<-textConnection(biv_mod)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('COL1' = TrueData[1,200,],
                                 'COL2' = TrueData[2,200,],
                                 'N'=dim(TrueData)[3]),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('alpha',
          'tau.disp')

##############################################
#Posterior Sampling
##############################################
posterior_samples<-coda.samples(model_jags, 
                                params, 
                                n.iter=10000)
posterior_samples.all<-do.call(rbind,posterior_samples)
#post1.summary<-summary(posterior_samples)
#post_means<-colMeans(posterior_samples.all)

post_means<-apply(posterior_samples.all, 2, median)
sample.labs<-names(post_means)
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
ci<-matrix(sprintf("%.1f",round(ci,1)), ncol=2)
row.names(ci)<-sample.labs
#post_means<-sprintf("%.1f",round(post_means,1))
#names(post_means)<-sample.labs

prev <- ilogit(post_means)
#Can recover acquisition rate as 1/prev*duration 
