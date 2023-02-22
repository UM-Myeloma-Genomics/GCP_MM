
#
# Probabilities of states in Phase 1
#

calculatePODInduction = function(TimeP1,SDeath,SP2){
  L = length(TimeP1)
  fd = rep(0,L) 
  for (j in 1:(L-1)){
    fd[j] = SDeath[j] - SDeath[j+1]
  }
  Pd = cumsum(fd*SP2)
  return (Pd)
} 


calculateMToP2 = function(TimeP1,SP2,SDeath){
  L = length(TimeP1)
  fp2 = rep(0,L) 
  for (j in 1:(L-1)){
    fp2[j] = SP2[j] - SP2[j+1]
  }
  Pp2 = cumsum(fp2*SDeath)
  return (Pp2)
} 


claculateAliveInInduction = function(PDeath,PendInduction){
  return (1-PDeath-PendInduction)
}



calculateDeathInPODInduction =function(TimeP4,SPODinIND,TimeP1,PrmoveP2){
  
  
  TimeY = seq(1,max(TimeP1)+max(TimeP4))
  L = length(TimeY)
  fp2 = rep(0,L) 
  for (j in 1:(length(TimeP1)-1)){
    fp2[j] = PrmoveP2[j+1] - PrmoveP2[j]
  }
  FProgrP2X = 1-SPODinIND
  FProgrP2 = rep(FProgrP2X[length(TimeP4)],L) 
  FProgrP2[1:length(TimeP4)] = FProgrP2X 
  
  PDeath = rep(0,L)
  
  for (t in 1:L){
    PDeath[t] = sum(fp2[1:t]*FProgrP2[t:1])
  }
  
  return (PDeath)
  
}


#
# Helper functions for competing evenets at P2
#


calculateDeathInP2Comp = function(TimeP1,SDeathP2,SProgr){
  L = length(TimeP1)
  fd = rep(0,L) 
  for (j in 1:(L-1)){
    fd[j] = SDeathP2[j] - SDeathP2[j+1]
  }
  Pd2 = cumsum(fd*SProgr)
  return (Pd2)
} 


calculateToProgressionComp = function(TimeP1,SProgr,SDeathP2){
  L = length(TimeP1)
  fp2 = rep(0,L) 
  for (j in 1:(L-1)){
    fp2[j] = SProgr[j] - SProgr[j+1]
  }
  Ppr = cumsum(fp2*SDeathP2)
  return (Ppr)
} 





#
# Probabilities of states in Phase 2
#


calculateDeathInP2 = function(TimeP2,SDeathP2,SProgr,TimeP1,SP2){
  TimeP3 = seq(1,max(TimeP1)+max(TimeP2))
  L = length(TimeP3)
  fp2 = rep(0,L) 
  for (j in 1:(length(TimeP1)-1)){
    fp2[j] = SP2[j+1] - SP2[j]
  }
  FDeathP2 = calculateDeathInP2Comp(TimeP2,SDeathP2,SProgr)
  FDeathP2st = rep(FDeathP2[length(TimeP2)],L) 
  FDeathP2st[1:length(TimeP2)] = FDeathP2 

  PDeathinP2 = rep(0,L)
  
  for (t in 1:L){
    PDeathinP2[t] = sum(fp2[1:t]*FDeathP2st[(t:1)])
  }
  
  return (PDeathinP2)
} 


calculateToProgression =function(TimeP2,SDeathP2,SProgr,TimeP1,SP2){

  TimeP3 = seq(1,max(TimeP1)+max(TimeP2))
  L = length(TimeP3)
  fp2 = rep(0,L) 
  for (j in 1:(length(TimeP1)-1)){
    fp2[j] = SP2[j+1] - SP2[j]
  }
  FDeathP2 = calculateToProgressionComp(TimeP2,SProgr,SDeathP2)
  FDeathP2st = rep(FDeathP2[length(TimeP2)],L) 
  FDeathP2st[1:length(TimeP2)] = FDeathP2 
  
  PProgression = rep(0,L)
  
  for (t in 1:L){
    PProgression[t] = sum(fp2[1:t]*FDeathP2st[(t:1)])
  }
  
  return (PProgression)
  
} 


claculateAliveInP2 = function(PrmoveP2,PDeathP2,PProgression){
  L = length(PrmoveP2)
  PrmoveP2X = rep(PrmoveP2[L],length(PDeathP2))
  PrmoveP2X[1:L] = PrmoveP2
  return (PrmoveP2X-PDeathP2-PProgression)
}

calculateDeathfromProgression = function(TimeP3,SDeathPr,TimeP2,SDeathP2,SProgr,TimeP1,SP2){
  TimeX = seq(1,max(TimeP1)+max(TimeP2))
  L = length(TimeX)
  fp2 = rep(0,L) 
  for (j in 1:(length(TimeP1)-1)){
    fp2[j] = SP2[j+1] - SP2[j]
  }
  FDeathP2 = calculateToProgressionComp(TimeP2,SProgr,SDeathP2)
  FDeathP2st = rep(FDeathP2[length(TimeP2)],L) 
  FDeathP2st[1:length(TimeP2)] = FDeathP2 
  
  PProgression = rep(0,L)
  
  for (t in 1:L){
    PProgression[t] = sum(fp2[1:t]*FDeathP2st[(t:1)])
  }
  
  
  TimeY = seq(1,max(TimeX)+max(TimeP3))
  L = length(TimeY)
  fp2 = rep(0,L) 
  for (j in 1:(length(TimeX)-1)){
    fp2[j] = PProgression[j+1] - PProgression[j]
  }
  FProgrP2X = 1-SDeathPr
  FProgrP2 = rep(FProgrP2X[length(TimeP3)],L) 
  FProgrP2[1:length(TimeP3)] = FProgrP2X 
  
  PDeath = rep(0,L)
  
  for (t in 1:L){
    PDeath[t] = sum(fp2[1:t]*FProgrP2[t:1])
  }
  
  return (PDeath)
}

