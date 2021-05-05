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

comm.infect.rate.k <- 1/50
comm.infect.rate.a <- 1/720
hh.infect.rate.k <- 1/30 #this is conditional on adult being colonized
hh.infect.rate.a <- 1/60 #conditional on kid being colonized
DurInf.k <-60
DurInf.a <- 21

TrueData.txn <-  replicate(1000, sirHH( time=400, 
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


Y.dat <- TrueData.txn[,200,]

t.Y <- t(Y.dat)
Y.cat <- rep(4, nrow(t.Y))
Y.cat[t.Y[,1]==1 & t.Y[,2]==1] <- 1
Y.cat[t.Y[,1]==1 & t.Y[,2]==0] <- 2
Y.cat[t.Y[,1]==0 & t.Y[,2]==1] <- 3


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
    #prev.cocol.hh <- prev.cocol.obs - prob.cocol.indep
    return(prev.cocol.obs)

})

ObsTab <- apply(co.colcat.time,2,function(x){
  #Take a time slice
  x <- factor(x, levels=c('0','1','2','3'))
  hh.status <- as.vector(table(x))
  
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
# if the duration of the index is longer, then duration of colonization is DurInfContact


if(DurInf.a<DurInf.k){
  dur.cocol.a <- DurInf.a 
  dur.cocol.k <- DurInf.a/2 #cocol of kids resulting from adult
}else{
  dur.cocol.a <- DurInf.k/2 
  dur.cocol.k <- DurInf.k #cocol of kids resulting from adult
}
expected.cocol <-   dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) + #Prev of kidresulting from txn from adult 
                    dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k) +    #Prev of adult resulting from txn from kid
                    (comm.infect.rate.a*DurInf.a)*(comm.infect.rate.k*DurInf.k) #Probability of co-occurrence if indepdenent
  hist(prev.cocol)
abline(v=expected.cocol)
##################################


#Now can estimate the 2x2 table with MLE?
#estimate: comm.infect.rate.a, comm.infect.rate.k,hh.infect.rate.k,hh.infect.rate.a
ml.simple <- function(params, DurInf.a=21, DurInf.k=60, Y ){
  hh.infect.rate.k <- ilogit(params[1])
  hh.infect.rate.a <- ilogit(params[2])
 # comm.infect.rate.k  <- ilogit(params[3])
  #comm.infect.rate.a <- ilogit(params[4])
     comm.infect.rate.a  <- 1/1200
     comm.infect.rate.k <- 1/7000
  
  #Define duration of co-colonization if index is parent or child
  if(DurInf.a<DurInf.k){
    dur.cocol.a <- DurInf.a 
    dur.cocol.k <- DurInf.a/2 #cocol of kids resulting from adult
  }else{
    dur.cocol.a <- DurInf.k/2 
    dur.cocol.k <- DurInf.k #cocol of kids resulting from adult
  }
  
  prob.col1 <-   comm.infect.rate.k*DurInf.k
  prob.col2 <-   comm.infect.rate.a*DurInf.a
  
  prob.cocol3 <-   dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) + #Prev of kidresulting from txn from adult 
    dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k)+     #Prev of adult resulting from txn from kid
    (comm.infect.rate.a*DurInf.a)*(comm.infect.rate.k*DurInf.k) #Probability of co-occurrence if indepdenentp

  prob.uncol0 <- max(0, (1 - (prob.col1 +prob.col2 + prob.cocol3)) )
  
  probs.est <- c(prob.col1,prob.col2, prob.cocol3, prob.uncol0)
  probs.est <- probs.est/sum(probs.est) + 1e-7
  
  LL <- dmultinom(Y, prob=probs.est, size=sum(Y), log=T)
  return(-LL)
}

parm.inits <- rep(log(0.001),4)

comm.infect.rate.k <- 1/1200
comm.infect.rate.a <- 1/7000
hh.infect.rate.k <- 1/120
hh.infect.rate.a <- 1/120

#parm.inits <- logit(c(1/100, 1/100, 1/1000, 1/1000))
parm.inits <- logit(c(1/100,1/100))
#init.lower <- rep(log(0.00001),4)
#init.upper <- rep(log(0.2),4)

mod1 <- nlm(ml.simple, parm.inits, Y=ObsTab[,200])
ilogit(mod1$estimate)

mod1 <- optim(par=parm.inits, fn=ml.simple,  Y=ObsTab[,200])
ilogit(mod1$par)


#To do in a regular regression setup--not sure this still makes sense unless we define as :
prob.col1 <-   comm.infect.rate.k*DurInf.k + dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) 
prob.col2 <-   comm.infect.rate.a*DurInf.a  + dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k)



#prob.cocol3 <-   dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) + #Prev of kidresulting from txn from adult 
 # dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k)+     #Prev of adult resulting from txn from kid
#  (comm.infect.rate.a*DurInf.a)*(comm.infect.rate.k*DurInf.k) #Probability of co-occurrence if indepdenentp

Xcat <- matrix( c(c(1,1), c(1,0), c(0,1), c(0,0)), nrow=4, byrow=T)

#Can w combine this with JAGS model, similar to https://esajournals.onlinelibrary.wiley.com/doi/10.1890/10-0173.1
##use the model to estimate prevalence of co-occurrence, and then back out
#the infection rate if we fix duration?
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12587
#see p118 here:https://academicworks.cuny.edu/cgi/viewcontent.cgi?article=4068&context=gc_etds
biv_mod2 <- "model{
  for (j in 1:n.pairs){
 
  #this formulation with the denom[i] ensures that the probs don't exceed 1;
  #as 
    denom[j] <- (comm.infect.rate.a*DurInf.a + 
                comm.infect.rate.k*DurInf.k + 
                (dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) + #Prev of kidresulting from txn from adult 
                  dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k)+     #Prev of adult resulting from txn from kid
                  (comm.infect.rate.a*DurInf.a)*(comm.infect.rate.k*DurInf.k) ) +
                  1)

    psi[1,j] <-  (dur.cocol.k*hh.infect.rate.k*(comm.infect.rate.a*DurInf.a) + #Prev of kidresulting from txn from adult 
                  dur.cocol.a*hh.infect.rate.a*(comm.infect.rate.k*DurInf.k)+     #Prev of adult resulting from txn from kid
                  (comm.infect.rate.a*DurInf.a)*(comm.infect.rate.k*DurInf.k) )/denom[j]
    psi[2,j] <- comm.infect.rate.k*DurInf.k/denom[j] #psi_10
    psi[3,j] <- comm.infect.rate.a*DurInf.a/denom[j] #psi_01
    psi[4,j] <-  1/denom[j] #psi_00
    

    y.cat[j] ~ dcat( psi[,j] )
    
    
  }
  
  DurInf.a <- 21
  DurInf.k <- 60

    dur.cocol.a <- DurInf.a 
    dur.cocol.k <- DurInf.a/2 #cocol of kids resulting from adult

  comm.infect.rate.k <- exp(a1)/(1+exp(a1)) #ilogit

  comm.infect.rate.a <- 1/7000 #FIX this based on HH without kids

  hh.infect.rate.k <- exp(c1)/(1+exp(c1)) #ilogit
  hh.infect.rate.a <- exp(d1)/(1+exp(d1)) #ilogit
  
  #For kids and adult, constrain the HH txn rate to be > community rate
  a1 ~dnorm(-3, 1/4)
  a1.diff ~ dnorm(0, 1)
  c1 <- a1 + exp(a1.diff)

  b1 <- logit(comm.infect.rate.a)
  b1.diff ~ dnorm(0, 1)
  d1 <- b1 + exp(b1.diff) 

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
model_spec<-textConnection(biv_mod2)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('y.cat' = Y.cat,
                                 'n.pairs' =ncol(Y.dat)),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('psi','comm.infect.rate.k',
          'comm.infect.rate.a','hh.infect.rate.k', 'hh.infect.rate.a','a1','b1','c1','d1')
##############################################
#Posterior Sampling
##############################################
posterior_samples<-coda.samples(model_jags, 
                                params, 
                                n.iter=10000)
posterior_samples.all<-do.call(rbind,posterior_samples)

post_means<-apply(posterior_samples.all, 2, median)
sample.labs<-names(post_means)
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
#ci<-matrix(sprintf("%.1f",round(ci,1)), ncol=2)
row.names(ci)<-sample.labs

#comm inf.a=0.000142 ;comm.k=0.00067; hh.a=0.018; kk.k=0.03
View(as.data.frame(post_means))

plot(logit(posterior_samples.all[,'comm.infect.rate.k']), type='l')
plot(logit(posterior_samples.all[,'comm.infect.rate.a']), type='l')
plot(logit(posterior_samples.all[,'hh.infect.rate.k']), type='l')
#plot(logit(posterior_samples.all[,'hh.infect.rate.a']), type='l')

plot(posterior_samples.all[,'comm.infect.rate.k'],posterior_samples.all[,'hh.infect.rate.k'])

#Community infection rate in kids and HH infection rate in adults not identifiable
plot(posterior_samples.all[,'comm.infect.rate.k'],posterior_samples.all[,'hh.infect.rate.a'])


plot(posterior_samples.all[,'a1'],posterior_samples.all[,'b1'])
plot(posterior_samples.all[,'a1'],posterior_samples.all[,'c1'])
plot(posterior_samples.all[,'a1'],posterior_samples.all[,'d1'])
plot(posterior_samples.all[,'b1'],posterior_samples.all[,'c1'])
plot(posterior_samples.all[,'b1'],posterior_samples.all[,'d1'])
plot(posterior_samples.all[,'d1'],posterior_samples.all[,'c1'])


plot(posterior_samples.all[,'psi[1,1]'],posterior_samples.all[,'psi[2,1]'])
plot(posterior_samples.all[,'psi[1,1]'],posterior_samples.all[,'psi[3,1]'])
plot(posterior_samples.all[,'psi[1,1]'],posterior_samples.all[,'psi[4,1]'])


plot(posterior_samples.all[,'c1']+posterior_samples.all[,'d1'], type='l')
plot(posterior_samples.all[,'a1']+posterior_samples.all[,'b1'], type='l')
plot(posterior_samples.all[,'a1']+posterior_samples.all[,'c1'], type='l')
plot(posterior_samples.all[,'b1']+posterior_samples.all[,'d1'], type='l')
plot(posterior_samples.all[,'a1']+posterior_samples.all[,'b1']+posterior_samples.all[,'c1'], type='l')


a_b_combined_post_samp <- posterior_samples.all[,'a1']+posterior_samples.all[,'b1']
c_d_combined_post_samp <- posterior_samples.all[,'c1']+posterior_samples.all[,'d1']
a_c_combined_post_samp <- posterior_samples.all[,'a1']+posterior_samples.all[,'c1']
b_d_combined_post_samp <- posterior_samples.all[,'b1']+posterior_samples.all[,'d1']

ilogit(hdi(a_b_combined_post_samp, credMass = 0.95))
