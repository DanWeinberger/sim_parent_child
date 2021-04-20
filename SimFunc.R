#Generate data from a simple transition model
##States: 0 uncolonized; 1=colonized (non-HH); 2= uncolonized,immune 3=colonized (HH acquistions)
gen.pair.data<-function(logit.acq.rate.kid= set.logit.acq.rate.kid,
                        logit.acq.rate.adult= set.logit.acq.rate.adult,
                        clear.rate.kid= 1/DurKid,
                        clear.rate.adult= 1/DurAdult,
                        DurKid=51,
                        DurAdult=19,
                        DurImmKid=90,
                        DurImmAdult=90, 
                        ntimes=365,
                        burn.days =100
                        ){
  
  acq.rate.kid <- exp(logit.acq.rate.kid )/ (1+exp(logit.acq.rate.kid))
  acq.rate.adult <- exp(logit.acq.rate.adult )/ (1+exp(logit.acq.rate.adult))
  
  prob.transmit.kid.day<- prob.transmit.kid/DurKid
  prob.transmit.adult.day<- prob.transmit.adult/DurAdult
    
  WaneRateKid  <-  1/DurImmKid
  WaneRateAdult <- 1/DurImmAdult
    state<-matrix(NA, ncol=ntimes, nrow=2)
    
    transmit.kid.parent <- rep(0, ntimes)
    transmit.parent.kid <- rep(0, ntimes)
    
    state[1,1] <- rbinom(n=1, size=1, prob=init.prev.kid)
    state[2,1] <- rbinom(n=1, size=1, prob=init.prev.adult)
    
    for(t in 2: ncol(state)){
      
      if(state[1,t-1]==0 & state[2,t-1]==0 ){ #both uncolonized (box 1)
        state[1,t] <- rbinom(n=1,size=1, prob=(acq.rate.kid ) )
        state[2,t] <- rbinom(n=1,size=1, prob=(acq.rate.adult ) )
        
      }else if (state[1,t-1]%in% c(1,3) & state[2,t-1]==0 ){ #kid colonized, adult uncolonized
        state[1,t] <- state[1,t-1]*rbinom(n=1,size=1, prob=(1-clear.rate.kid ) )
        if(state[1,t]==0){state[1,t]==2} #if transition from col to uncol, go to immune class
        state[2,t] <- rbinom(n=1,size=1, prob=(prob.transmit.kid.day + acq.rate.adult ) )
        
        if(state[2,t]==1){state[2,t] <-3 } #HH acquistion
        
      }else if (state[1,t-1]==0 & state[2,t-1]%in% c(1,3) ){ #adult colonized, kid uncolonized
        state[1,t] <- rbinom(n=1,size=1, prob=( prob.transmit.adult.day + acq.rate.kid ) )
        state[2,t] <- state[2,t-1]*rbinom(n=1,size=1, prob=(1-clear.rate.adult ) )
        if(state[2,t]==0){state[2,t]==2} #if transition from col to uncol, go to immune class
        
        if(state[1,t]==1){state[1,t] <- 3 } #HH acquistion
        
        
      }else if (state[1,t-1]%in% c(1,3) & state[2,t-1]%in% c(1,3) ){ #kid colonized, adult colonized
        state[1,t] <- state[1,t-1]*rbinom(n=1,size=1, prob=(1-clear.rate.kid ) )
        if(state[1,t]==0){state[1,t]==2} #if transition from col to uncol, go to immune class
        
        state[2,t] <- state[2,t-1]*rbinom(n=1,size=1, prob=(1 - clear.rate.adult ) )
        if(state[2,t]==0){state[2,t]==2} #if transition from col to uncol, go to immune class
      
      ##IMMUNE CATEGORIES  
      }else if (state[1,t-1]==2 & state[2,t-1] %in% c(1,3) ){ #kid immune, adult colonized
        recover.kid <- rbinom(n=1,size=1, prob=(1 - WaneRateKid ) )
        if(recover.kid){
          state[1,t]=0 #go back to uncolonized
        }else{
          state[1,t]=2
        }
        state[2,t] <- rbinom(n=1,size=1, prob=(1 - clear.rate.adult ) )
        if(state[2,t]==0){state[2,t]==2} #if transition from col to uncol, go to immune class
  
      }else if (state[1,t-1]%in% c(1,3) & state[2,t-1]==2 ){ #adult immune, kid colonized
      recover.adult <- rbinom(n=1,size=1, prob=(1 - WaneRateAdult ) )
      if(recover.adult){
        state[2,t]=0 #go back to uncolonized
      }else{
        state[2,t]=2
      }
      state[1,t] <- rbinom(n=1,size=1, prob=(1 - clear.rate.kid ) )
      if(state[1,t]==0){state[1,t]==2} #if transition from col to uncol, go to immune class
    
    }else if (state[1,t-1]==2 & state[2,t-1]==2 ){ #both immune
      recover.kid <- rbinom(n=1,size=1, prob=(1 - WaneRateKid ) )
      if(recover.kid){
        state[1,t]=0 #go back to uncolonized
      }else{
        state[1,t]=2
      }
      recover.adult <- rbinom(n=1,size=1, prob=(1 - WaneRateAdult ) )
      if(recover.adult){
        state[2,t]=0 #go back to uncolonized
      }else{
        state[2,t]=2
      }
      
    }else if (state[1,t-1]==2 & state[2,t-1]==0 ){ #kid immune, adult uncolonized
      recover.kid <- rbinom(n=1,size=1, prob=(1 - WaneRateKid ) )
      if(recover.kid){
        state[1,t]=0 #go back to uncolonized
      }else{
        state[1,t]=2
      }
      state[2,t] <- rbinom(n=1,size=1, prob=( acq.rate.adult ) )
    }else if (state[1,t-1]==2 & state[2,t-1]==0 ){ #adult immune, kid uncolonized
      state[1,t] <- rbinom(n=1,size=1, prob=( acq.rate.kid ) )
      
      recover.adult <- rbinom(n=1,size=1, prob=(1 - WaneRateAdult ) )
      if(recover.adult){
        state[2,t]=0 #go back to uncolonized
      }else{
        state[2,t]=2
      }
    }
    }

    Individual_States <- state[,-(1:burn.days)]
    Individual_States[Individual_States==3] <- 1
  return(Individual_States)
}
