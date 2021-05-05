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

post_means<-apply(posterior_samples.all, 2, median)
sample.labs<-names(post_means)
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
ci<-matrix(sprintf("%.1f",round(ci,1)), ncol=2)
row.names(ci)<-sample.labs
#post_means<-sprintf("%.1f",round(post_means,1))
#names(post_means)<-sample.labs
