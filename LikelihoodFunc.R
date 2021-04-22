LL.hh <- function(par,obs){
  simDat1 <-  t(replicate(10000, sirHH(time=101, burn.days=100,
                                        logit.beta=par[c(1,2)] ,
                                        logit.lambda=par[c(3,4)]  
  ), simplify='array'))
  Y.hat <- as.vector(table(factor(simDat1[,1],levels=c('0','1')), 
                           factor(simDat1[,2],levels=c('0','1'))))

  pi <- Y.hat/sum(Y.hat)  # + 1e-7 #add small constant to prevent infinite LL if 0 count is found
  LL <- dmultinom(x=Y,  prob=pi, log = TRUE) #fits obs Y|pi
  return(-LL)
}


