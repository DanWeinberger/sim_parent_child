#4 states: neither colonized; child colonized, parent colonized, both colonized
#In Okie dokie, most common serotypes ~12% prevalence

source('SimFunc.R')
set.seed(123)
DurKid <- 60
DurAdult <- 18 #Results in 54% clearing in 14 days
set.acq.rate.adult <- 0.1/365
set.acq.rate.kid <- 0.6/365

init.prev.kid=0.12
init.prev.adult= init.prev.kid/2

prob.transmit.kid <-0.70 #prob transmit to HH member
prob.transmit.adult <- prob.transmit.kid/10 #prob transmit to HH member

#Note this model isn't quite right--there is no immunity built in. If a kid was just colonized with a strain, they are not likely to immediately be re-colonized

#Generates simulated data from a Markov model for 365 time points
# At each time point, can be 1: neither colonized; 2: child colonized; 3: adult colonized; 4: both coonized
#Each column is a household
sim.data <- replicate(100,gen.pair.data())
kid.colonized <- matrix(sim.data %in% c(2,4), ncol=ncol(sim.data))
adult.colonized <- matrix(sim.data %in% c(3,4), ncol=ncol(sim.data))

kid.prevalence <- apply(kid.colonized,1,mean)
adult.prevalence <- apply(adult.colonized,1,mean)

plot(kid.prevalence)
abline(h=0.12)

plot(adult.prevalence)
abline(h=0.2)
