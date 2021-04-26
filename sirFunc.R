#HH transitions
sirHH <- function( time=180, 
                   logit.beta, #=logit(c(1/100, 1/300)) , #comunity infection rate for kids and adults
                      durR=c(60,60), #duration of protection from subsequent infection
                      durInf=c(60,21), #duration for ids and adults
                      burn.days=100
)  {
  
  beta <- ilogit(logit.beta)
  
  delta <- 1/durInf
  mu <- 1/durR
  
  X <- array(NA, dim=c(2,3,time))
  dimnames(X)[[1]] <- c('kid','adult')
  dimnames(X)[[2]] <- c('S','I','R')
  
  prev.init <- beta*durInf
  prev.init[prev.init>1] <- 1

  X[,c('I'),1] <- rbinom(2,1,prev.init) #initialize at inidence*duration (prev)
  X[,c('S'),1] <- 1 - X[,c('I'),1]
  X[,c('R'),1] <- 0
  #X[,,1]
  
  for(i in 2:time){
    for(j in 1:2){
      
      other.person <- ifelse(j==1,2,1)
      
      #probability of escaping infection for a susceptible person depends on community transmision rate and HH rate
      prob_escape_infect <- ( X[j,c('S'),i-1]*(1-beta[j]) )
      
      # Prob a susceptible person stays susceptible: S[i-1]*(1-beta[j])*(1-lambda[j])^other.person
      prob_S <- prob_escape_infect + #Susceptible person Escape infection from HH if HH member is infected
        X[j,c('R'),i-1]*mu[j] #If perosn was R on previous, what is probability of being S
      
      # Prob a susceptible person is infected: S[i-1]*(1- (1-beta[j])*(1-lambda[j])^other.person)
      prob_I <- X[j,c('S'),i-1]*(1-prob_escape_infect) +  #prob a susceptible person becomes infected 
        X[j,c('I'),i-1]*(1-delta[j]) #Prob infected person Stay in I
      
      prob_R <- X[j,c('I'),i-1]*delta[j] + #uninfected but immune
        X[j,c('R'),i-1]*(1-mu[j]) #stay recovered
      
      probs <- c(prob_S,prob_I,prob_R)
      
      if(is.na(sum(probs))){
        X[j,,i] <- X[j,,i-1]
      }else{
        X[j,,i] <- rmultinom(1,1,probs) #which state is person in?
      }
    }
    
    
  }
  return(X[,'I',-c(1:burn.days)])
}

