#Generate data from a simple transition model
gen.pair.data<-function(acq.rate.kid= set.acq.rate.kid,
                        acq.rate.adult= set.acq.rate.adult,
                        clear.rate.kid= 1/DurKid,
                        clear.rate.adult= 1/DurAdult
                        ){
  prob.transmit.kid.day<-prob.transmit.kid/DurKid
  prob.transmit.adult.day<-prob.transmit.kid/DurAdult
    
  ntimes=365
    state<-matrix(NA, ncol=ntimes, nrow=2)
    
    #Row1=kid
    #Row2=parent
    
    state[1,1] <- rbinom(n=1, size=1, prob=init.prev.kid)
    state[2,1] <- rbinom(n=1, size=1, prob=init.prev.adult)
    
    for(t in 2: ncol(state)){
      
      if(state[1,t-1]==0 & state[2,t-1]==0 ){ #both uncolonized
        state[1,t] <- rbinom(n=1,size=1, prob=(acq.rate.kid ) )
        state[2,t] <- rbinom(n=1,size=1, prob=(acq.rate.adult ) )
        
      }else if (state[1,t-1]==1 & state[2,t-1]==0 ){ #kid colonized, adult uncolonized
        state[1,t] <- rbinom(n=1,size=1, prob=(1-clear.rate.kid ) )
        state[2,t] <- rbinom(n=1,size=1, prob=(prob.transmit.kid.day + acq.rate.adult ) )
       
      }else if (state[1,t-1]==0 & state[2,t-1]==1 ){ #adult colonized, kid uncolonized
        state[1,t] <- rbinom(n=1,size=1, prob=( prob.transmit.adult.day + acq.rate.kid ) )
        state[2,t] <- rbinom(n=1,size=1, prob=(1-clear.rate.adult ) )
          
      }else if (state[1,t-1]==1 & state[2,t-1]==1 ){ #kid colonized, adult colonized
        state[1,t] <- rbinom(n=1,size=1, prob=(1-clear.rate.kid ) )
        state[2,t] <- rbinom(n=1,size=1, prob=(1 - clear.rate.adult ) )
      }
    }
    
    class <- apply(state,2, function(x){
      if(x[1]==0 & x[2]==0){ z=1}
      if(x[1]==1 & x[2]==0){ z=2}
      if(x[1]==0 & x[2]==1){ z=3}
      if(x[1]==1 & x[2]==1){ z=4}
      return(z)
      })
  return(class)
}
