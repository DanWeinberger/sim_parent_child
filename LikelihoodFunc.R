LL.hh <- function(par,obs){
  simDat1 <-  replicate(1000, sirHH(time=100, burn.days=50,
                                    logit.beta=par[c(1,2)] ,

  ), simplify='array')
  
  
  pi1 <- apply(simDat1,2, Tab.fun) #prob at each time points
  pi2 <- apply(pi1,1,mean) #ave prob
  
  LL <- dmultinom(x=Y,  prob=pi2, log = TRUE) #fits obs Y|pi
  return(-LL)
}


Tab.fun <-  function(x){
  Y.hat <- as.vector(table(factor(x[1,],levels=c('0','1')), 
                           factor(x[2,],levels=c('0','1'))))
  pi <- Y.hat/sum(Y.hat)   + 1e-7 #add small constant to prevent infinite LL if 0 count is found
  return(pi)
}
