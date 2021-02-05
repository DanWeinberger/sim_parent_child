source('SimFunc.R')
set.seed(123)
#Many of these parameters baed on Melegaro https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2870123/

#these 2 values just initialize data but aren't critical to get right
init.prev.kid=0.5 #initial prevalence for the strain of interest
init.prev.adult= init.prev.kid/4

#Probability of transmitting to HH members if colonized
prob.transmit.kid <- 0.5 #prob kid transmits to HH member
prob.transmit.adult <- prob.transmit.kid/10 #prob adult transmits to HH member


#Generates simulated data from a Markov model for 365 time points
# At each time point, can be 1: neither colonized; 2: child colonized; 3: adult colonized; 4: both colonized
#Each column is a household
rr.func <- function(x){
  rr<- ((sum(x==4)+0.5)/sum(x %in% c(2,4)))/((sum(x==3)+0.5)/sum(x %in% c(1,3)))
  return(rr)
}

try.dur <- c(2,5,8,10,15,20,25,30,40,50,57,60,65,70,80,90,100,200,300)

rr.all <- sapply( try.dur, function(x) { 
   res1<- replicate(1000,gen.pair.data(acq.rate.kid= 0.012,
                                        acq.rate.adult= 0.004,
                                        clear.rate.kid= 1/x,
                                        clear.rate.adult= 1/19,
                                        DurImmKid=90,
                                        DurKid=x,
                                        DurAdult=19,
                                        DurImmAdult=90), simplify=F)
   hh.data <- sapply(res1,'[[','HH.States')
   
   rr.data1 <- apply(hh.data,1,rr.func)
   mean.rr.data1 <- exp(mean(log(rr.data1)))
   return(mean.rr.data1)
})

#calculate true value for how much more common acquisition is for adults when a kid is present
total.acq.adult <- (0.004 + 0.012*prob.transmit.kid)
true.acq.rr <- total.acq.adult/0.004

plot(try.dur, rr.all, xlab='Duration', ylab='Est RR', type='l')
abline(h=true.acq.rr)
abline(v=67, lty=2, col='gray')
abline(v=57, lty=2, col='gray')

#"true" RR in adults  is 2.5 (adults with kid in HH have 2.5x more acquisitions)
#Where the RR estimate matches the true, the duration in kids is ~67 (duration in adults set to 19)
0.004/19*10000 #adults duration vs acquisiton rate
0.012/65*10000
0.012/57*10000 #duration in kids 3x duration in adults
##this doesn't quite fit --why?