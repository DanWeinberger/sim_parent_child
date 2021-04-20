#In this simulation, we generate the colonization status each day for a child and their parent
#We specify the community acuistion rate and the probability of transmitting
#to a HH member if you are colonized
#We track the household colonization status in 4 states: neither colonized; child colonized, parent colonized, both colonized
#We also track individual colonization status as uncolonized, colonized from community acquisition or colonized from HH contact


source('SimFunc.R')
set.seed(123)
#Many of these parameters baed on Melegaro https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2870123/
DurKid <- 51 # 2 month duration
DurAdult <- 19 #Results in 54% clearing in 14 days

#these 2 values just initialize data but aren't critical to get right
init.prev.kid=0.5 #initial prevalence for the strain of interest
init.prev.adult= init.prev.kid/4

#These are a bit arbitrary
#COMMUNITY acquisition rates
set.acq.rate.kid <- 0.012
set.acq.rate.adult <- set.acq.rate.kid/3 #based on melegaro

set.logit.acq.rate.kid = log(set.acq.rate.kid/(1-set.acq.rate.kid))
set.logit.acq.rate.adult = log(set.acq.rate.adult/(1-set.acq.rate.adult))

#Probability of transmitting to HH members if colonized
prob.transmit.kid <- 0.5 #prob kid transmits to HH member
prob.transmit.adult <- prob.transmit.kid/10 #prob adult transmits to HH member

#Generates simulated data from a Markov model for 365 time points
# At each time point, can be 1: neither colonized; 2: child colonized; 3: adult colonized; 4: both colonized
#Each column is a household

#600 people, adult and child
TrueData <-  replicate(600, gen.pair.data(ntimes=180, burn.days=100), simplify='array')
TrueData50 <- t(TrueData[,50,] )
Y <- as.vector(table(TrueData50[,1], TrueData50[,2]))

LL.hh <- function(parms){
  simDat1 <-  t(replicate(sum(Y), gen.pair.data(ntimes=101, burn.days=100,logit.acq.rate.kid=parms[1], logit.acq.rate.adult=parms[2]), simplify='array'))
  Y.hat <- as.vector(table(simDat1[,1], simDat1[,2]))
  pi <- Y.hat/sum(Y.hat)
  
  LL <- dmultinom(x=Y,  prob=pi, log = TRUE)
  return(-LL)
}

mod1 <- nlm(f=LL.hh, c(alpha=0, beta=0 ) )

