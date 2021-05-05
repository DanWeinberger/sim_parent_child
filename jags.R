#In this simulation, we generate the colonization status each day for a child and their parent
#We specify the community acuistion rate and the probability of transmitting
#to a HH member if you are colonized
#We track the household colonization status in 4 states: neither colonized; child colonized, parent colonized, both colonized
#We also track individual colonization status as uncolonized, colonized from community acquisition or colonized from HH contact


source('sirFunc.R')
source('LikelihoodFunc.R')
source('jags_mod_hh.R')

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
                                    burn.days=1,
                                    hh.txn=F
), 
simplify='array')

#Check if TrueData is at steady state
prev <- apply(TrueData,c(1,2), mean )

plot(prev[1,])
plot(prev[2,])


#SIMPLE JAGS MODEL, ASSUMING NO HH TXN
#Can recover acquisition rate as 1/prev*duration 
mod1 <- biv.mod.func(Y1=TrueData[1,200,], Y2=TrueData[2,200,])

comm.infect.rate.k <- 1/1200
comm.infect.rate.a <- 1/7000
hh.infect.rate.k <- 1/120
hh.infect.rate.a <- 1/120
DurInf.k <-60
DurInf.a <- 21
TrueData.txn <-  replicate(5000, sirHH( time=400, 
                                    logit.beta= logit(c(comm.infect.rate.k, comm.infect.rate.a ))  , #community infection rate for kids and adults #Make this Very small to represent ST-specfici risk
                                    durR=c(60,60), #waning of protection from subsequent infection
                                    durInf=c(DurInf.k,DurInf.a), #1/duration of colonization for kids and adults
                                    burn.days=1,
                                    logit.lambda=log(hh.infect.rate.k, hh.infect.rate.a),
                                    hh.txn=T
), 
simplify='array')

#Check if TrueData is at steady state
prev <- apply(TrueData.txn,c(1,2), mean )

plot(prev[1,])
plot(prev[2,])
mod2 <- biv.mod.func(Y1=TrueData.txn[1,200,], Y2=TrueData.txn[2,200,])

#expected prevalence of colonization resulting from HH transmission?
#HH txn:
#lambda[1]*DurInf[1]*PREV_adult + lambda[2]*DurInf[2]*PREV_kids
#First part gives prevalence in kids resulting from txn from adults 2nd part of the reverse
1/120*60*(1/7000*21) + (1/120)*21*(1/1200*60)

p1 <-ilogit(post_means[1])
p2 <-ilogit(post_means[2])
p3 <-ilogit(post_means[1] + post_means[2] + post_means[3])

#################################################################
#################################################################
##Try this outside a model##
co.colcat.time <- apply(TrueData.txn,2, function(x){
  colcat <- rep(NA, nrow(x))
  colcat[x[1,]==0 & x[2,]==0] <- 0
  colcat[x[1,]==1 & x[2,]==0] <- 1
  colcat[x[1,]==0 & x[2,]==1] <- 2
  colcat[x[1,]==1 & x[2,]==1] <- 3
  return(colcat)
})

prev.cocol <- apply(co.colcat.time,2,function(x){
    #Take a time slice
    x <- factor(x, levels=c('0','1','2','3'))
    hh.status <- as.vector(table(x))
    #Prevalence person 1
    prev.person1 <- sum(hh.status[c(2,4)])/sum(hh.status)
    
    #Prevalence person 2
    prev.person2<- sum(hh.status[c(3,4)])/sum(hh.status)
    
    #expected co-occurrence if they were independent
    prob.cocol.indep <- prev.person1 * prev.person2
    
    #Observed co-occurrence
    prev.cocol.obs <- sum(hh.status[c(4)])/sum(hh.status)
    
    #Co-occurrence that (maybe) arises from HH exposure
    prev.cocol.hh <- prev.cocol.obs - prob.cocol.indep
    return(prev.cocol.hh)

})


#expected prevalence of colonization resulting from HH transmission?
#HH txn:
#lambda[1]*DurInf[1]*PREV_adult + lambda[2]*DurInf[2]*PREV_kids
#But this formula ignores the fact that the index case might clear their colonization before
#contact does. So we need to adjust for that. On average, if we assume
#risk of transmission is uniform between time=0 and time=DurInf,
#then the average remaining duration AFTER transmission is DurInf/2
#If duration of the index is shorter than duration of the contact
#then duration when they are co-colonized is DurInfIndex/2, otherwise
# if the duration of the index is longer, then duration of colonization is DurInfIndex


if(DurInf.a<DurInf.k){
  dur.cocol.a <- DurInf.a 
  dur.cocol.k <- DurInf.a/2 #cocol of kids resulting from adult
}else{
  dur.cocol.a <- DurInf.k/2 
  dur.cocol.k <- DurInf.k #cocol of kids resulting from adult
}
expected.cocol <-   dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) + #Prev of kidresulting from txn from adult 
                    dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k)     #Prev of adult resulting from txn from kid
    

hist(prev.cocol)
abline(v=expected.cocol)
##################################

dist.dur <-rexp(10000, 1/60)
day.transmit <- sapply(dist.dur, function(x) runif(n=1,min=0, max=x))
remain.days <- dist.dur - day.transmit #how many days after transmit is index colonized? 
