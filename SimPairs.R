#In this simulation, we generate the colonization status each day for a child and their parent
#We specify the community acuistion rate and the probability of transmitting
#to a HH member if you are colonized
#We track the household colonization status in 4 states: neither colonized; child colonized, parent colonized, both colonized
#We also track individual colonization status as uncolonized, colonized from community acquisition or colonized from HH contact


source('sirFunc.R')
source('sirFuncGen.R') #Generate observed 2x2 table of co-occurrence
source('LikelihoodFunc.R')
library(dplyr)
library(pbapply)
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
TrueData <-  replicate(2000, sirHHGen( time=180, 
                                  logit.beta=logit(c(1/150, 1/(365*2))/10  ) , #comunity infection rate for kids and adults #Make this Very small to represent ST-specfici risk
                                  logit.lambda= logit(c(0.1*1/120, 0.5*1/120)), ##H infection rate for adult-kid and kid-adults  #HH txn rate per day from adult to child prev*prob.transmit.day
                                  logit.mu=logit(c(1/30,1/30)), #waning of protection from subsequent infection
                                  logit.delta=logit(c(1/60,1/21)), #1/duration for kids and adults
                                  burn.days=100
                                          ), 
                       simplify='array')
TrueData80 <- t(TrueData[,80,] )
#table(TrueData80[,1],TrueData80[,2])
Y <- as.vector(table(factor(TrueData80[,1],levels=c('0','1')), factor(TrueData80[,2],levels=c('0','1'))))

#Set constraints for likelihood to prevent extreme values (causes problems for optimization)
# lower.prob <- logit(0.0001)
# upper.prob <- logit(0.9999)
# parms.l <- rep(lower.prob,length(parms))
# parms.u <- rep(upper.prob,length(parms))

parms <-  c(0,0,0,0)
#parms <- c(logit(c(1/150, 1/(365*2))), logit(c(0.1*1/120, 0.5*1/120)) ) #feed in correct parms
ptm <- proc.time()
  #mod1 <- optim(parms,LL.hh, lower=parms.l, obs=Y, upper=parms.u, method='L-BFGS-B' )
  mod1 <- nlm(LL.hh, parms, obs=Y )

proc.time() - ptm

parm.est <- mod1$estimate
1/ilogit(parm.est)


#TEST holding adult constant
parms <-  logit(rep(0.001,4))
#parms <- c(logit(c(1/150, 1/(365*2))), logit(c(0.1*1/120, 0.5*1/120)) ) #feed in correct parms
ptm <- proc.time()
#mod1 <- optim(parms,LL.hh, lower=parms.l, obs=Y, upper=parms.u, method='L-BFGS-B' )
mod1 <- nlm(LL.hh.test, parms, obs=Y, steptol=1e-3 )

proc.time() - ptm

parm.est <- mod1$estimate
ilogit(parm.est)
