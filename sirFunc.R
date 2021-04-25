#HH transitions
sirHH <- function( times=60, 
                   logit.beta, #=logit(c(1/100, 1/300)) , #comunity infection rate for kids and adults
                   logit.lambda, #= logit(c(1/600,1/60)), ##H infection rate for adult-kid and kid-adults
                   logit.mu=logit(c(1/10000,1/10000)), #waning of protection from subsequent infection
                   logit.delta=logit(c(1/60,1/21)), #1/duration for ids and adults
                   burn.days=0
)  {

  beta <- ilogit(logit.beta)
  lambda <- ilogit(logit.lambda)
  delta <- ilogit(logit.delta)
  mu <- ilogit(logit.mu)
  
  X <- array(NA, dim=c(2,3,times))
  probs <- array(NA, dim=c(2,3,times))
  
    dimnames(X)[[1]] <- c('kid','adult')
  dimnames(X)[[2]] <- c('S','I','R')
  
  dimnames(probs)[[1]] <- c('kid','adult')
  dimnames(probs)[[2]] <- c('S','I','R')

   X[,c('I'),1] <- rbinom(2,1,beta)
   X[,c('S'),1] <- 1 - X[,c('I'),1]
   X[,c('R'),1] <- 0

  for(i in 2:times){
   for(j in 1:2){
    
    other.person <- ifelse(j==1,2,1)
    
    #probability of escaping infection for a susceptible person depends on community transmision rate and HH rate
    prob_escape_infect <- ( X[j,c('S'),i-1]*(1-beta[j]) )  #Susceptible person escape infection from community

     # Prob a susceptible person stays susceptible: S[i-1]*(1-beta[j])*(1-lambda[j])^other.person
    prob_S <- prob_escape_infect  #Susceptible person Escape infection from HH if HH member is infected

    # Prob a susceptible person is infected: S[i-1]*(1- (1-beta[j])*(1-lambda[j])^other.person)
    prob_I <- X[j,c('S'),i-1]*(1-prob_escape_infect) +  #prob a susceptible person becomes infected 
                        X[j,c('I'),i-1]*(1-delta[j]) #Prob infected person Stay in I
 
    prob_R <- X[j,c('I'),i-1]*(delta[j]) +
                 X[j,c('R'),i-1] #if R=1 at previous period, stay in 
    
    probs[j,,i] <- c(prob_S,prob_I, prob_R)
    
    if(is.na(sum(probs[j,,i])) | sum(probs[j,,i]) ==0   ){
      X[j,,i] <-  0 #X[j,,i-1]
    }else{
    X[j,,i] <- rbinom(1,1,probs[j,,i]) #which state is person in?
   }
 
   }
    
    #Once a person is infected or immune, they are censored--shouldn't contribute to S
    #censorR <- which(X[,'R',]==1) 
    #censorI <- which(X[,'I',]==1) 
    
    # X[,'S',][censorR] <- NA
    # X[,'S',][censorI] <- NA
    # X[,'I',][censorR] <- NA
    
  } 
  return(X[,'I',-c(1:burn.days)])
}


