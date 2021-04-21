#In this simulation, we generate the colonization status each day for a child and their parent
#We specify the community acuistion rate and the probability of transmitting
#to a HH member if you are colonized
#We track the household colonization status in 4 states: neither colonized; child colonized, parent colonized, both colonized
#We also track individual colonization status as uncolonized, colonized from community acquisition or colonized from HH contact


source('sirFunc.R')
library(dplyr)
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

#600 people, adult and child
TrueData <-  replicate(600, sirHH( time=180, 
                                  logit.beta=logit(c(1/60, 1/300)) , #comunity infection rate for kids and adults
                                  logit.lambda= logit(c(0.00167,0.0167)), ##H infection rate for adult-kid and kid-adults
                                  logit.mu=logit(c(1/60,1/60)), #waning of protection from subsequent infection
                                  logit.delta=logit(c(1/60,1/30)), #1/duration for ids and adults
                                  burn.days=100
                                          ), 
                       simplify='array')
TrueData50 <- t(TrueData[,50,] )
Y <- as.vector(table(TrueData50[,1], TrueData50[,2]))

LL.hh <- function(parms){
  simDat1 <-  t(replicate(sum(Y), sirHH(time=101, burn.days=100,logit.beta=parms[c(1,2)] ,
                                                            logit.lambda=parms[c(3,4)]  
                                        ), simplify='array'))
  Y.hat <- as.vector(table(simDat1[,1], simDat1[,2]))
  pi <- Y.hat/sum(Y.hat)
  LL <- dmultinom(x=Y,  prob=pi, log = TRUE) #fits obs Y|pi
  return(-LL)
}

#Set constraints to prevent extreme values (causes problems for optimization)
lower.prob <- log(0.00001/(1- 0.00001))
upper.prob <- log(0.99999/(1- 0.99999))
parms <-  c(0,0,0,0 )
parms.l <- rep(lower.prob,length(parms))
parms.u <- rep(upper.prob,length(parms))

ptm <- proc.time()
mod1 <- optim(parms,LL.hh, lower=parms.l, upper=parms.u, method='L-BFGS-B' )
proc.time() - ptm

parms <- mod1$par
ilogit(parms)
