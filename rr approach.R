source('SimFunc.R')
set.seed(123)
#Many of these parameters baed on Melegaro https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2870123/

#Generates simulated data from a Markov model for 365 time points
# At each time point, can be 1: neither colonized; 2: child colonized; 3: adult colonized; 4: both colonized
#Each column is a household
sim.data1 <- replicate(100,gen.pair.data(acq.rate.kid= 0.012,
                                        acq.rate.adult= 0.004,
                                        clear.rate.kid= 1/51,
                                        clear.rate.adult= 1/19,
                                        DurImmKid=90,
                                        DurKid=51,
                                        DurAdult=19,
                                        DurImmAdult=90), simplify=F)
#V2: shorter duration
sim.data2 <- replicate(100,gen.pair.data(acq.rate.kid= 0.012,
                                         acq.rate.adult= 0.004,
                                         clear.rate.kid= 1/19,
                                         clear.rate.adult= 1/19,
                                         DurKid=19,
                                         DurAdult=19,                                        
                                         DurImmKid=90,
                                         DurImmAdult=90), simplify=F)
sim.data3 <- replicate(100,gen.pair.data(acq.rate.kid= 0.012,
                                         acq.rate.adult= 0.004,
                                         clear.rate.kid= 1/100,
                                         clear.rate.adult= 1/19,
                                         DurImmKid=90,
                                         DurKid=100,
                                         DurAdult=19,
                                         DurImmAdult=90), simplify=F)
sim.data4 <- replicate(100,gen.pair.data(acq.rate.kid= 0.012,
                                         acq.rate.adult= 0.004,
                                         clear.rate.kid= 1/75,
                                         clear.rate.adult= 1/19,
                                         DurImmKid=90,
                                         DurKid=75,
                                         DurAdult=19,
                                         DurImmAdult=90), simplify=F)

hh.data1 <- sapply(sim.data1,'[[','HH.States')
hh.data2 <- sapply(sim.data2,'[[','HH.States')
hh.data3 <- sapply(sim.data3,'[[','HH.States')
hh.data4 <- sapply(sim.data4,'[[','HH.States')

rr.func <- function(x){
  rr<- ((sum(x==4)+0.5)/sum(x %in% c(2,4)))/((sum(x==3)+0.5)/sum(x %in% c(1,3)))
  return(rr)
}

rr.data1 <- apply(hh.data1,1,rr.func)
mean.rr.data1 <- exp(mean(log(rr.data1)))

rr.data2 <- apply(hh.data2,1,rr.func)
mean.rr.data2 <- exp(mean(log(rr.data2)))

rr.data3 <- apply(hh.data3,1,rr.func)
mean.rr.data3 <- exp(mean(log(rr.data3)))

rr.data4 <- apply(hh.data4,1,rr.func)
mean.rr.data4 <- exp(mean(log(rr.data4)))


mean.rr.data1
mean.rr.data2
mean.rr.data3
mean.rr.data4

dur <- c(51,19,100,75)
rrs <- c(mean.rr.data1,mean.rr.data2,mean.rr.data3, mean.rr.data4)
plot(dur, rrs)


