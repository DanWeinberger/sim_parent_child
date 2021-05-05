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
    

    x[j] ~ dcat( psi[,j] )
    
    
  }
  
  DurInf.a <- 21
  DurInf.k <- 60
  #if(DurInf.a<DurInf.k){
    dur.cocol.a <- DurInf.a 
    dur.cocol.k <- DurInf.a/2 #cocol of kids resulting from adult
  #}else{
  #  dur.cocol.a <- DurInf.k/2 
  #  dur.cocol.k <- DurInf.k #cocol of kids resulting from adult
  #}
  
  comm.infect.rate.k <- exp(a1)/(1+exp(a1)) #ilogit
 # comm.infect.rate.a <- exp(b1)/(1+exp(b1)) #ilogit
  comm.infect.rate.a <- 1/7000 #FIX this based on HH without kids

  hh.infect.rate.k <- exp(c1)/(1+exp(c1)) #ilogit
  hh.infect.rate.a <- exp(d1)/(1+exp(d1)) #ilogit
  
  #For kids and adult, constrain the HH txn rate to be > community rate
  a1 ~dnorm(-3, 1/4)
  a1.diff ~ dnorm(0, 1)
  c1 <- a1 + exp(a1.diff)

  b1 <- logit(comm.infect.rate.a)
 # b1 ~ dnorm(logit(1/7000),1/4)
  b1.diff ~ dnorm(0, 1)
  d1 <- b1 + exp(b1.diff) 

  #Observation model
  for (j in 1:n.pairs){
    for(i in 1:n.groups){
      logit(p[i,j]) = p0[i]    
      y[j,i] ~ dbern(Xcat[x[j],i]*p[i,j]) #Xcat is a mat with 2^i rows and i columns
    }
  }
  
  #Priors
  for(i in 1:n.groups){
    p0[i] ~ dnorm(0,1e-1)
  }

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
                       data=list('y' =t(Y.dat),
                                 'Xcat'=Xcat,
                                 'n.pairs' =ncol(Y.dat),
                                 'n.groups'=2),
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
