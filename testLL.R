#test.ll

LL.hh.test <- function(par,obs){
  betas <-rnorm(2,0,4)
  lambdas <-rnorm(2,0,4)
  simDat1 <-  t(replicate(100, sirHH(time=101, burn.days=100,
                                        logit.beta=betas,
                                        logit.lambda=lambdas  
  ), simplify='array'))
  Y.hat <- as.vector(table(factor(simDat1[,1],levels=c('0','1')), 
                           factor(simDat1[,2],levels=c('0','1'))))
  if(sum(Y.hat)>0 ){
    pi <- Y.hat/sum(Y.hat)  # + 1e-7 #add small constant to prevent infinite LL if 0 count is found
    LL <- dmultinom(x=Y,  prob=pi, log = TRUE) #fits obs Y|pi
  }else{
    LL=0
  }
  out.list=list('negLL'=-LL,'logit.betas'=betas,'logit.lambda'=lambdas)
  return(out.list)
}

test1 <- pbreplicate(100, LL.hh.test(), simplify=F)

ll <- sapply(test1,'[[', 'negLL')

betas <- as.data.frame(t(ilogit(sapply(test1,'[[', 'logit.betas'))))
names(betas) <- c('beta1','beta2')
  
lambdas <- as.data.frame(t(ilogit(sapply(test1,'[[', 'logit.lambda'))))
names(lambdas) <- c('lambda1','lambda2')

compare.parms <- cbind.data.frame(ll, betas, lambdas)

#This makes it looks like beta1 and beta2 are not identifiable and 
#la
compare.parms$comb1 <- compare.parms$beta1*compare.parms$lambda1
compare.parms$comb2 <- compare.parms$beta2*compare.parms$lambda2
compare.parms$comb3 <- compare.parms$lambda1*compare.parms$lambda2
compare.parms$comb4 <- compare.parms$beta1*compare.parms$beta2






