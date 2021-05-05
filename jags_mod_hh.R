biv.mod.func <- function(Y1, Y2){

  biv_mod.hh <- "model{
  for(i in 1:N){
    COL1[i] ~ dbern (P1[i]) 
    COL2[i] ~ dbern (P2[i])
    COL12[i] ~ dbern (P3[i])
    
    
    logit(P1[i]) <- alpha[1]  #+phi[i]    
    logit(P2[i]) <- alpha[2]  #+phi[i]    
    
    logit(P3[i]) <-  alpha[1]  + alpha[2] + alpha[3] #+phi[i]

     phi[i] ~dnorm(0, tau.disp)
   
  }
  
  
  

  # Prior distribution of model parameters
  for(k in 1:3){
   alpha[k] ~ dnorm(0,0.001)
  }
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
model_spec<-textConnection(biv_mod.hh)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('COL1' = Y1,
                                 'COL2' = Y2,
                                'COL12' = Y1*Y2,
                                 'N'=length(Y1)),
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

cocol.prob <- ilogit(post_means['alpha[1]']+post_means['alpha[2]']+post_means['alpha[3]'])

#How many co-colonization events result from 
prev.txn <- cocol.prob - ilogit(post_means['alpha[1]'])*ilogit(post_means['alpha[1]'] )


out.list <- list('post_means'=post_means,'ci'=ci,'cocol.prob'=cocol.prob)
return(out.list)
}
