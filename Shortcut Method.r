

# This file contains functions to estimate the infection rate and
# the recovery rate in a partially known network using the Shortcut method.

# Last modified date: 25 September 2023

# Authors: Omar De La Cruz Cabrera and Razan Alsehibani.



## Simulation parameters
n = 200  # Number of individuals
m = 200  # Number of time steps
u = 100    # Number of repetitions 


## SIS parameters
beta  =  0.04 # Rate of infection
gamma =  0.02 # Rate of recovery


## UK parameters
alpha = 0.007 # Rate of contact at random
delta = 0.006 # Rate of contact through network
rao = alpha + delta # Rate of contact at random or through network 

start_time <- Sys.time()

library(igraph)

AB2 <- readMM(file="./matrixData/AdjERSimple_Edges_0.1.mtx")   

AE <- AB2*1 
AAB <- graph.adjacency(AE, mode="undirected", weighted=NULL) 

listGammaHatE <- vector("integer", u)
candBE = numeric(u)
valueBE = numeric(u)
for (h in 1:u){
  
  
  sE = c(1,rep(0,n-1) )  # initial status (1:Infected, 0:Susceptible)
  kE = c(1, rep(0,n-1) ) # initial status (1:Known, 0:Unknown)
  
  
  KE = matrix(0,nrow = m,ncol = n)
  KE[1,] = kE
  
  
  
  
  SE = matrix(0,nrow = m,ncol = n)
  SE[1,] = sE
  
  
  
  for (t in 2:m){   ## loop over time
    
    
    for (i in 1:n){  ## loop over nodes
      
      if (sE[i] == 0 && kE[i] == 0){
        probOfNoInfection = exp(-beta*sum(AE[i,]*ifelse(sE==1,1,0)))  # "ifelse" is not needed for SIS; but it would be needed for SIR
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size=1,
                       prob = c(probOfNoInfection,1-probOfNoInfection))
        
        probOfNotContacted = exp(-rao*sum(AE[i,]))
        kE[i] = sample(x = 0:1,
                       size=1,replace = FALSE, 
                       prob = c(probOfNotContacted,1-probOfNotContacted))
        
      }
      
      else    if (sE[i] == 1 && kE[i] == 0 ){
        probOfNoRecovery = exp(-gamma)
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size = 1,
                       prob = c(1-probOfNoRecovery,probOfNoRecovery))
        
        probOfNotContacted = exp(-rao*sum(AE[i,]))
        kE[i] = sample(x = 0:1,
                       size=1,
                       prob = c(probOfNotContacted,1-probOfNotContacted))
        
      }
      
      else if (sE[i] == 0 && kE[i] == 1){
        probOfNoInfection = exp(-beta*sum(AE[i,]*ifelse(sE==1,1,0)))  # "ifelse" is not needed for SIS; but it would be needed for SIR
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size=1,
                       prob = c(probOfNoInfection,1-probOfNoInfection))
        
      }
      
      else    if (sE[i] == 1 && kE[i] == 1 ){
        probOfNoRecovery = exp(-gamma)
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size = 1,
                       prob = c(1-probOfNoRecovery,probOfNoRecovery))
        
      }
      
      
    }  ## i loop
    
    
    SE[t,] = sE
    KE[t,] = kE
    print(t)
    
    
  }  ## t loop 
  
  
  
  
  
  order(colSums(KE)) ## The order of individuals based on the discovery date
  KsE=subset(KE,select = order(colSums(KE))) ## Sort K based on the order above
  SsE=subset(SE,select = order(colSums(KE))) ## Sort S based on the order above
  Ks2E <- replace(KsE, KsE == 0, -1) ## replace any unknown(0) with (-1)
  final2E <- replace(SsE, Ks2E == -1, -1) ## replace the infectious status for any unknown individual in Ss with (-1)
  
  niE = length(which(final2E==1)) ## number of infected individuals in final2E 
  nsE = length(which(final2E==0)) ## number of susceptible individuals in final2E
  totE = niE+nsE ## total number of individuals with known status
  prevE = niE / totE ## the percentage of infected individuals 
  
  finalpE<-replace(final2E, final2E==-1, prevE) ## replace the unknown values with prevelance
  
  ## claculate gamma
  counterE = 0 ## counter of how many times an individual remains I given it was I.
  
  for(t in 3:m){
    
    for(i in 1:n){
      if (final2E[t-1,i]== 1 && final2E[t,i]== 1){
        counterE = counterE + 1
      }
      
      
    } # i loop
  } #t loop  
  
  
  gammaHatE = log(niE/counterE)
  listGammaHatE[[h]] <- gammaHatE
  
  RRE = matrix(0,nrow = m,ncol = n)  ## RRE is a matrix that contains the total infection risk from neighbours for every individual at each time step
  rrE = c(1,rep(0,n-1))
  RRE[1,] = rrE
  
  for (t in 2:m){   ## loop over time
    
    for (i in 1:n){  ## loop over nodes
      
      rrE =  finalpE[t,]%*%AE
      RRE[t,] = rrE[1,]
      
      
    }
  }
  
  
  
  
  EstiBetaE = function(beta){
    prob = 0
    sumlog = 0
    
    for (t in 2:m){   ## loop over time
      
      for (i in 1:n){  ## loop over nodes
        
        
        if (final2E[t-1,i]== 0 && final2E[t,i]== 0){
          prob = exp(-beta*rrE[,i])
          sumlog = sumlog + log(prob)
          
        }
        
        else if (final2E[t-1,i]== 0 && final2E[t,i]== 1){
          prob = (1- exp(-beta*rrE[,i]))
          sumlog = sumlog + log(prob)
          
        }
      }# i loop 
    } #t loop 
    
    return(sumlog)
  }
 
  
  opBE = optimize(EstiBetaE,c(0,1), maximum = TRUE)
  candBE[h]= opBE$maximum
  valueBE[h]= opBE$objective
  
} ## h loop


MLE = candBE[which.max(valueBE)]
plot(candBE)

end_time <- Sys.time()
Total_Time = end_time - start_time
