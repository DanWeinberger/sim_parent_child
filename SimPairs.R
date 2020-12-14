#In this simulation, we generate the colonization status each day for a child and their parent
#We specify the community acuistion rate and the probability of transmitting
#to a HH member if you are colonized
#We track the household colonization status in 4 states: neither colonized; child colonized, parent colonized, both colonized
#We also track individual colonization status as uncolonized, colonized from community acquisition or colonized from HH contact


source('SimFunc.R')
set.seed(123)
DurKid <- 60 # 2 month duration
DurAdult <- 18 #Results in 54% clearing in 14 days

#these 2 values just initialize data but aren't critical to get right
init.prev.kid=0.12 #initial prevalence for the strain of interest
init.prev.adult= init.prev.kid/2

#These are a bit arbitrary
#COMMUNITY acquisition rates
set.acq.rate.adult <- 0.1/365
set.acq.rate.kid <- 0.6/365

#Probability of transmitting to HH members if colonized
prob.transmit.kid <-0.5 #prob kid transmits to HH member
prob.transmit.adult <- prob.transmit.kid/10 #prob adult transmits to HH member

#Generates simulated data from a Markov model for 365 time points
# At each time point, can be 1: neither colonized; 2: child colonized; 3: adult colonized; 4: both coonized
#Each column is a household
sim.data <- replicate(100,gen.pair.data(), simplify=F)

HH.States <- sapply(sim.data,'[[', 'HH.States')
Indiv.States <- sapply(sim.data,'[[','Individual_States', simplify='array')

kid.colonized <- matrix(HH.States %in% c(2,4), ncol=ncol(HH.States))
adult.colonized <- matrix(HH.States %in% c(3,4), ncol=ncol(HH.States))

kid.prevalence <- apply(kid.colonized,1,mean)
adult.prevalence <- apply(adult.colonized,1,mean)

plot(kid.prevalence)
abline(h=0.12)

plot(adult.prevalence)
abline(h=0.2)


#sORT OUT which colonization events are acquired from parent or from 


#1=uncolonized; 2=kid colonized; 3=adult colonized; 4=both colonized
table(HH.States)


#Prevalence, by source of acquisition (3=HH, 1=community)
table(Indiv.States[1,,]) #Kids
table(Indiv.States[2,,]) #Adults

#Prevalence, by source of acquisition, when colonization detected in both (3=HH, 1=community)

table(Indiv.States[1,,]*(HH.States==4 ) )  #kids
table(Indiv.States[2,,]*(HH.States==4 ) )  #Adults


