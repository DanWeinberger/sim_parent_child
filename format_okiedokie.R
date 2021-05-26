library(readxl)
library(reshape2)
library(dplyr)
library(lme4)
library(rjags)
library(HDInterval)

a1a <- read_excel('./DONOTSYNC/20210311_overview.xlsx', sheet='OK-2')
a1a$study= 'OK2'

a1b <- read_excel('./DONOTSYNC/20210311_overview.xlsx', sheet='OK-3')
a1b$study= 'OK3'

a1c <- read_excel('./DONOTSYNC/20210311_overview.xlsx', sheet='OK-4')
a1c$study= 'OK4'

a1 <- bind_rows(a1a,a1b,a1c)

a1$shared_carriage <-  NULL
a1$shared_serotype <- NULL
a1$serotype <- NULL
a1$Household <- paste(1:nrow(a1), a1$study,sep='_')
a1$study <- NULL
a1.m <- melt(a1, id.vars = 'Household')
a1.m$parent <- 0
a1.m$parent[substr(a1.m$variable,1,1) == 'p'] <- 1
a1.m$variable <-  sub("_[^_]+$", "", a1.m$variable)

a1.m <- a1.m[!(a1.m$variable %in% c('kCE_CS','pCE_CS')),]

a1.m$st <- substring(a1.m$variable,4)
a1.m <- a1.m[,c('Household','parent','st','value')]

a1.m$st[substr(a1.m$st,1,1)=='6'] <- '6'

a1.c <- acast(a1.m, Household ~ parent ~ st, fun.aggregate = sum, na.rm=T)
a1.c[a1.c>1] <- 1 #if a parent is positive for both OP and saliva, just count 1


#Estimate prevalence by serotype and age
prev <- t(apply(a1.c,c(2,3) , mean, na.rm=T))

#Estimate N colonized by serotype and age
ColonizedN <- t(apply(a1.c,c(2,3) , sum, na.rm=T))

#RESHAPED DATA FOR ANALYSIS
b1.m <- melt(a1.c)
b1.c <- dcast(b1.m, Var3+Var1  ~ Var2)
names(b1.c) <- c('st','HH','child','parent')

#This tests effect of child colonization on parent colonization, by serotype
#There is probably a more direct approach to test this with a multivariate outcome
mod1 <- glmer(parent ~ (child|st), family='binomial', data=b1.c)
summary(mod1)
rand1 <- ranef(mod1)$st
rand1 <- rand1[order(rand1$child),]
View(rand1)

#As per Josh's suggestion, try a bivariate logistic regression with shared random intercept
mod2 <- "model{
  for (i in 1:n.hh){
    for (j in 1:n.sts){

  #Likelihood for observed data from kids and adults
  y[i,1,j] ~ dbern(mu[i,1,j])
  y[i,2,j] ~ dbern(mu[i,2,j])

    #Alpha0=global intercept
    #alpha1=global effect of being an adult
    #kappa is serotype-specific prevalence
    #delta is serotype-specific effect of being an adult
    #theta= random intercept to capture individual-level
 # logit(mu[i,1,j]) <- (alpha0 +          kappa[j] + theta[i,j] ) #Child
 # logit(mu[i,2,j]) <- (alpha0 + alpha1 + kappa[j] + delta[j] + theta[i,j] ) #adult
 # theta[i,j] ~ dnorm(0,prec.theta) #should this by centered on 0

  logit(mu[i,1,j]) <- (kappa[j] + theta[i,j]  ) #Child
  logit(mu[i,2,j]) <- (kappa[j] + delta[j] + theta[i,j]   ) #adult

  theta[i,j] ~ dnorm(0,prec.theta) #effect for seeing same serotype in same HH


    }
  }

  for (j in 1:n.sts){
    kappa[j] <- alpha0 + kappa.disp[j] #serotype-specific intercept
   kappa.disp[j] ~ dnorm(0, prec.kappa)
   delta.disp[j] ~ dnorm(0, prec.delta)
   delta[j] <- alpha1 + delta.disp[j]
  }
  
  alpha0 ~dnorm(0,1e-4)
  alpha1 ~dnorm(0,1e-4)

  prec.kappa ~ dgamma(0.001, 0.001)
  prec.theta ~ dgamma(0.001, 0.001)
  prec.delta ~ dgamma(0.001, 0.001)

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
model_spec<-textConnection(mod2)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1),
                       data=list('y' = a1.c,
                                 'n.hh' =dim(a1.c)[[1]],
                                 'n.sts' =dim(a1.c)[[3]]
                                 ),
                       n.adapt=1000, 
                       n.chains=1)

params<-c('kappa', 'delta','theta','prec.theta','prec.kappa','prec.delta')
##############################################
#Posterior Sampling
##############################################
posterior_samples<-coda.samples(model_jags, 
                                params, 
                                n.iter=10000)
posterior_samples.all<-do.call(rbind,posterior_samples)
saveRDS(posterior_samples.all,'posterior_samples_jags.rds')



plot(posterior_samples.all[,'kappa[1]'], posterior_samples.all[,'delta[1]'])
plot(posterior_samples.all[,'kappa[1]'], type='l')
plot(posterior_samples.all[,'delta[1]'], type='l')
plot(posterior_samples.all[,'alpha0'], type='l')
plot(posterior_samples.all[,'alpha1'], type='l')

#If we combine Alpha0 + kappa, get a nicely mixed intercept!
plot(posterior_samples.all[,'alpha0']+posterior_samples.all[,'kappa[1]'], type='l')

#If we combine alpha1 + delta, get a nicely mixed serotype-specific effect of adult prev!
plot(posterior_samples.all[,'alpha1']+posterior_samples.all[,'delta[1]'], type='l')
plot(posterior_samples.all[,'alpha1']+posterior_samples.all[,'delta[2]'], type='l')


plot(posterior_samples.all[,'prec.theta'], type='l')


post_means<-apply(posterior_samples.all, 2, median)
sample.labs<-names(post_means)
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
#ci<-matrix(sprintf("%.1f",round(ci,1)), ncol=2)
row.names(ci)<-sample.labs

#Serotype effect for adults--how much more common is the serotype in adults vs kids
delta.ci <- ci[grep('delta[', sample.labs, fixed=T),]
row.names(delta.ci) <- dimnames(a1.c)[[3]]

#Serotype-specific intercept--just gives estimate of relative prevalence, in kids
kappa.ci <- ci[grep('kappa[', sample.labs, fixed=T),]
row.names(kappa.ci) <- dimnames(a1.c)[[3]]

#this shows a bimodal distributio of theta--first is tightly clustered around
#0, presumably these are HH with no infections, second distribution
#is around ~1.5, presumably HH withat least one infection
#Could filter this out and estimate theta if kid colonized or adult colonized
post.means.theta <-post_means[grep('theta[', sample.labs, fixed=T)]
hist(post.means.theta)
hist(post.means.theta[post.means.theta>0.5])
post.means.theta.mat <- matrix(post.means.theta, nrow=dim(a1.c)[[1]])
pos.child <- a1.c[,1,]
pos.adult <- a1.c[,2,]
hist(post.means.theta.mat[pos.child==1])
hist(post.means.theta.mat[pos.adult==1])

#For each serotype calculate mean
##Overall
mean(post.means.theta.mat[pos.adult==1])
mean(post.means.theta.mat[pos.child==1])
pos.mean.st <- t(sapply(1:ncol(pos.child), function(x){
  pos.mean.adult.st <- mean(post.means.theta.mat[pos.adult[,x]==1,x])
  pos.kid.adult.st <- mean(post.means.theta.mat[pos.child[,x]==1,x])
  out.list=c('pos.mean.adult.st'=pos.mean.adult.st,'pos.kid.adult.st'=pos.kid.adult.st)
  out.list[is.nan(out.list)] <- 999
  return(out.list)
  }))
pos.mean.st <- cbind.data.frame(pos.mean.st,'st'=dimnames(a1.c)[[3]])

plot(pos.mean.st[,1], pos.mean.st[,2], xlim=c(0,5), ylim=c(0,5))
abline(a=0,b=1)