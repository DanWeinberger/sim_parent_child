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
                                        DurImmAdult=90), simplify=F)
sim.data2 <- replicate(100,gen.pair.data(acq.rate.kid= 0.012,
                                         acq.rate.adult= 0.004,
                                         clear.rate.kid= 1/51,
                                         clear.rate.adult= 1/19,
                                         DurImmKid=90,
                                         DurImmAdult=90), simplify=F)
