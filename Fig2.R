#setwd("set to current folder")


rm(list = ls())

#### Viral dynamics model


SARS_COV_model <- function(t,x,params){
  with(as.list(x),{   
    
    ddt_S = -beta*V*S   # susceptible target cells
    ddt_I = beta*V*S - delta*I^k*I - m*E^r*I/(E^r+E50^r)   # Infected cells
    ddt_V = p*I-c*V  # Viral load
    
    ddt_M1 = w*I*M1-q*M1
    ddt_M2 = q*(M1-M2)# M's are the intermediate compartments of the immune response that gives rise to effector cells (E) in response to viral infection 
    ddt_E = q*M2 - de*E # Effector cells that kill the infected cells
    
    der <- c(ddt_S,ddt_I,ddt_V,ddt_M1,ddt_M2,ddt_E)
    
    list(der)
  })       
}




library(deSolve)
library(DEoptim)


## Predifining colors ID for simulations
coloresID = c("steelblue4","slateblue3","navy","lightsteelblue4","dodgerblue2", 
              "cadetblue2","cornflowerblue","deepskyblue3","royalblue1","skyblue3","midnightblue",
              
              "salmon4","indianred3","brown2","coral1","chocolate2",
              "firebrick1","lightcoral","tomato2","darkred",
              
              "magenta4",
              
              "aquamarine3","darkolivegreen3","cadetblue4","chartreuse3","forestgreen","seagreen3",
              "palegreen3","springgreen4","darkseagreen4","green4","greenyellow","lightgreen")

### Specifying the marker type for each patient to be used in the simulation 
pchID = c(21,22,23,24,25,
          21,22,23,24,25,
          21,22,23,24,25,
          21,22,23,24,25,
          21,22,23,24,25)

#### Load Viral load data file
data <- read.csv("Viral_Loads.csv", header=TRUE,stringsAsFactors=FALSE) # Viral_Loads.csv is the file that records viral loads observed for each patient over the course of their infection
datanew = data

## Determine all unique Patient IDs in the dataset
IDs = unique(data[,"ID"])



### storing data in a temp file datanew
datanew = data


#### Defining final fig structure
nrows = 1
ncolumns =2

par(mfrow = c(nrows,ncolumns),     
    oma = c(0, 0, 0, 0), # 
    mar = c(3.5, 3.5, 0.0, 0.0), 
    mgp = c(2.5, 1, 0),    
    xpd = FALSE)
ld= round(10^min(data$VL))


#### Loading Parameter set that fits each patient
PAR_early_all=read.csv("Estimated_params.csv", header=TRUE,stringsAsFactors=FALSE)
# Estimated_params.csv is the file that has columns -  tzero,	log10beta,	delta,	k,	log10p,	m	,log10w,	E50,	r,	q,	de,	vx,	c,	cov_CvsNC,	cov_study





############################################################################################################
########################################################################   FIG A
############################################################################################################


tend = max(data$dao)
ijk=1
for(ID in IDs){   # [c(2,3)]){
  i = which(IDs %in% ID)
  print(i)
  
  PAR_early = PAR_early_all[i,] # Finding parameters estimated for a patient with identifier ID
  
  #Loading all parameters for plotting simulated data
  tzero=PAR_early$tzero                 # time of infection (before first positive when t=0)	
  beta=10^PAR_early$log10beta           # Virus Infectivity 
  delta=PAR_early$delta                 # density dependent death coefficient   
  k=PAR_early$k                         # Power density dependent rate   
  p=10^PAR_early$log10p                 # Virus production rate 
  w=10^PAR_early$log10w                 # Effector Precursor proliferation rate
  E50=PAR_early$E50                     # maximum killing rate by effector cells 
  r=PAR_early$r                         # Hill coefficient - killing rate
  q=PAR_early$q                         # Precursors to effector "differentiation" rate 
  de=PAR_early$de                       # Death rate of effector cells
  c=PAR_early$c                         # Virus Clearance rate
  
  
  # We also predefined the observed patients in two categories based upon whether the viral suppression is observed in the observed duration or not
  # If teh clearance is observed then Clearance_or_not=1 (which allows Effector cells to act upon infection and clear it) else Clearance_or_no=0 (means effector cell population was either not stimuated or they lacked the power to kill infected cells resolving infection)
  Clearance_or_not=PAR_early$cov_CvsNC
  
  if(Clearance_or_not==1){
    m=PAR_early$m
  }else{
    m=0
  }
  
  # As different studies analyzed different samples, we accounted the difference in the volume of sample using parameter "vx" d
  Study_no=PAR_early$cov_study  # Study ID = 1 for singapore (NASAL), 2 for German (SPETUM), 3 for Korea (NASAL) , 4 for France(Europe)  (NASAL)
  
  if (Study_no==2){
    vx=10^PAR_early$vx  # German patients (Study_no=2) had their spetum analyzed and we assume that the vol of spetum is vx times the volume of nasal swabs
  }
  else{
    vx=1  # nasal swabs were employed in studies 1,3 and 4, and we use the volume of their sample to be 1
  }
  
  
  #####  Initial conditions for Viral dynamics model
  S_0 = 1e7
  I_0 = 1
  V_0 = p*I_0/c
  M1_0 = 1
  M2_0 = 0
  M3_0 = 0
  E_0 = 0
  
  # Solving ODE
  names(S_0)=""
  names(I_0)=""
  init.x <- c(S=S_0,I=I_0,V=V_0,M1=M1_0, M2=M2_0, E=E_0)
  t.out <- seq(0,30,by=0.01)  
  
  params=c()
  # Solving Viral Dynamics ODE MODEL  
  
  out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
  
  
  # Plotting
  I = out$I
  E = out$E
  
  pkI = delta*I^k ## Innate immune response
  pkE = m*E^r/(E^r+E50^r)  ## Late (acquired) immune response
  
  if(ijk==1 ){ ## Innate immune response for the first patient
    
    plot(out$time,pkI, bty="n", type="l",xaxt="n",yaxt="n",
         lwd=4,col="steelblue3",
         xlim=c(0,25),
         ylim=c(0,12),
         cex.main=1.3,cex.axis=1.3,cex.lab=1.3,
         xlab="Days post infection",
         ylab="Killing rates (1/day)",
    )
    axis(1,at=axTicks(1),labels=as.character(axTicks(1)),cex.axis=1.3,font=1,lwd=1,las=1)
    axis(2,at=seq(0,12,length.out = 4),labels=as.character(seq(0,12,length.out = 4)),
         cex.axis=1.3,font=1,lwd=1,las=1)

    
    text(5.4,10,"Mediated by:",adj=0,cex=1.0)
    legend(5,9.5,c("Early (innate) response","Late (acquired) response"),bty="n",
           col=c("steelblue3","aquamarine3"),lwd=4,cex=1.0)
    mtext("a)",adj=-0.2,line=0.2,font=1,outer=FALSE,cex=1.1)
   
  }
  else{  ## Innate immune response for remaining patients
    lines(out$time,pkI,col="steelblue3",
          lwd=4)
  }
  
  
  ijk=ijk+1
}
for(ID in IDs){   ### Late (acquired) immune response for all patients
  i = which(IDs %in% ID)
  print(i)
  
  PAR_early = PAR_early_all[i,] # Finding parameters estimated for a patient with identifier ID
  
  #Loading all parameters for plotting simulated data
  tzero=PAR_early$tzero                 # time of infection (before first positive when t=0)	
  beta=10^PAR_early$log10beta           # Virus Infectivity 
  delta=PAR_early$delta                 # density dependent death coefficient   
  k=PAR_early$k                         # Power density dependent rate   
  p=10^PAR_early$log10p                 # Virus production rate 
  w=10^PAR_early$log10w                 # Effector Precursor proliferation rate
  E50=PAR_early$E50                     # maximum killing rate by effector cells 
  r=PAR_early$r                         # Hill coefficient - killing rate
  q=PAR_early$q                         # Precursors to effector "differentiation" rate 
  de=PAR_early$de                       # Death rate of effector cells
  c=PAR_early$c                         # Virus Clearance rate
  
  
  # We also predefined the observed patients in two categories based upon whether the viral suppression is observed in the observed duration or not
  # If teh clearance is observed then Clearance_or_not=1 (which allows Effector cells to act upon infection and clear it) else Clearance_or_no=0 (means effector cell population was either not stimuated or they lacked the power to kill infected cells resolving infection)
  Clearance_or_not=PAR_early$cov_CvsNC
  
  if(Clearance_or_not==1){
    m=PAR_early$m
  }else{
    m=0
  }
  
  # As different studies analyzed different samples, we accounted the difference in the volume of sample using parameter "vx" d
  Study_no=PAR_early$cov_study  # Study ID = 1 for singapore (NASAL), 2 for German (SPETUM), 3 for Korea (NASAL) , 4 for France(Europe)  (NASAL)
  
  if (Study_no==2){
    vx=10^PAR_early$vx  # German patients (Study_no=2) had their spetum analyzed and we assume that the vol of spetum is vx times the volume of nasal swabs
  }
  else{
    vx=1  # nasal swabs were employed in studies 1,3 and 4, and we use the volume of their sample to be 1
  }
  
  
  #####  Initial conditions for Viral dynamics model
  S_0 = 1e7
  I_0 = 1
  V_0 = p*I_0/c
  M1_0 = 1
  M2_0 = 0
  M3_0 = 0
  E_0 = 0
  
  # Solving ODE
  names(S_0)=""
  names(I_0)=""
  init.x <- c(S=S_0,I=I_0,V=V_0,M1=M1_0, M2=M2_0, E=E_0)
  t.out <- seq(0,30,by=0.01)  
  
  params=c()
  # Solving Viral Dynamics ODE MODEL  
  
  out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
  
  
  # Plotting
  I = out$I
  E = out$E
  
  pkI = delta*I^k ## Innate immune response
  pkE = m*E^r/(E^r+E50^r)  ## Late (acquired) immune response

  lines(out$time,pkE,col="aquamarine3",
          lwd=4)

  ijk=ijk+1
}








############################################################################################################
########################################################################   FIG B
############################################################################################################



ijk=1
for(ID in IDs){   
  i = which(IDs %in% ID)
  print(i)
  
  PAR_early = PAR_early_all[i,] # Finding parameters estimated for a patient with identifier ID
  
  #Loading all parameters for plotting simulated data
  tzero=PAR_early$tzero                 # time of infection (before first positive when t=0)	
  beta=10^PAR_early$log10beta           # Virus Infectivity 
  delta=PAR_early$delta                 # density dependent death coefficient   
  k=PAR_early$k                         # Power density dependent rate   
  p=10^PAR_early$log10p                 # Virus production rate 
  w=10^PAR_early$log10w                 # Effector Precursor proliferation rate
  E50=PAR_early$E50                     # maximum killing rate by effector cells 
  r=PAR_early$r                         # Hill coefficient - killing rate
  q=PAR_early$q                         # Precursors to effector "differentiation" rate 
  de=PAR_early$de                       # Death rate of effector cells
  c=PAR_early$c                         # Virus Clearance rate
  
  
  # We also predefined the observed patients in two categories based upon whether the viral suppression is observed in the observed duration or not
  # If teh clearance is observed then Clearance_or_not=1 (which allows Effector cells to act upon infection and clear it) else Clearance_or_no=0 (means effector cell population was either not stimuated or they lacked the power to kill infected cells resolving infection)
  Clearance_or_not=PAR_early$cov_CvsNC
  
  if(Clearance_or_not==1){
    m=PAR_early$m
  }else{
    m=0
  }
  
  # As different studies analyzed different samples, we accounted the difference in the volume of sample using parameter "vx" d
  Study_no=PAR_early$cov_study  # Study ID = 1 for singapore (NASAL), 2 for German (SPETUM), 3 for Korea (NASAL) , 4 for France(Europe)  (NASAL)
  
  if (Study_no==2){
    vx=10^PAR_early$vx  # German patients (Study_no=2) had their spetum analyzed and we assume that the vol of spetum is vx times the volume of nasal swabs
  }
  else{
    vx=1  # nasal swabs were employed in studies 1,3 and 4, and we use the volume of their sample to be 1
  }
  
  
  #####  Initial conditions for Viral dynamics model
  S_0 = 1e7
  I_0 = 1
  V_0 = p*I_0/c
  M1_0 = 1
  M2_0 = 0
  M3_0 = 0
  E_0 = 0
  
  # Solving ODE
  names(S_0)=""
  names(I_0)=""
  init.x <- c(S=S_0,I=I_0,V=V_0,M1=M1_0, M2=M2_0, E=E_0)
  t.out <- seq(0,30,by=0.01)  
  
  params=c()
  # Solving Viral Dynamics ODE MODEL  
  
  out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
  
  
  # Plotting
  I = out$I
  E = out$E
  
  pkI = delta*I*I^k ## Innate immune response
  pkE = m*I*E^r/(E^r+E50^r)  ## Late (acquired) immune response
  
  if(ijk==1 ){
    


    plot(out$time,pkI, bty="n", type="l",xaxt="n",yaxt="n",log="y",
         lwd=4,col="steelblue3",#"steelblue3",
         xlim=c(0,25),
         #ylim=c(1,25),
         ylim=c(1,1e8),
         cex.main=1.3,cex.axis=1.3,cex.lab=1.3,
         xlab="Days post infection",
         ylab="Infected cells cleared per day"
    )
    axis(1,at=axTicks(1),labels=as.character(axTicks(1)),cex.axis=1.3,font=1,lwd=1,las=1)

    axis(2,at=c(1,1e2,1e4,1e6,1e8),labels=expression("1","10"^"2","10"^"4","10"^"6","10"^"8"),
         cex.axis=1.3,font=1,lwd=1,las=1)

    text(6.4,1.5e7,"Mediated by:",adj=0,cex=1.0)
    legend(6.2,1e7,c("Early (innate) response","Late (acquired) response"),bty="n",
           col=c("steelblue3","aquamarine3"),lwd=4,cex=1.0)
    # 
    mtext("b)",side=3,adj=-0.3,line=0.01,font=1,outer=FALSE,cex=1.1)
    # 
    
  }
  else{
    lines(out$time,pkI,col="steelblue3",
          lwd=4)
  }
  
  
  ijk=ijk+1
}

ijk=1
for(ID in IDs){   # [c(2,3)]){
  i = which(IDs %in% ID)
  print(i)
  
  PAR_early = PAR_early_all[i,] # Finding parameters estimated for a patient with identifier ID
  
  #Loading all parameters for plotting simulated data
  tzero=PAR_early$tzero                 # time of infection (before first positive when t=0)	
  beta=10^PAR_early$log10beta           # Virus Infectivity 
  delta=PAR_early$delta                 # density dependent death coefficient   
  k=PAR_early$k                         # Power density dependent rate   
  p=10^PAR_early$log10p                 # Virus production rate 
  w=10^PAR_early$log10w                 # Effector Precursor proliferation rate
  E50=PAR_early$E50                     # maximum killing rate by effector cells 
  r=PAR_early$r                         # Hill coefficient - killing rate
  q=PAR_early$q                         # Precursors to effector "differentiation" rate 
  de=PAR_early$de                       # Death rate of effector cells
  c=PAR_early$c                         # Virus Clearance rate
  
  
  # We also predefined the observed patients in two categories based upon whether the viral suppression is observed in the observed duration or not
  # If teh clearance is observed then Clearance_or_not=1 (which allows Effector cells to act upon infection and clear it) else Clearance_or_no=0 (means effector cell population was either not stimuated or they lacked the power to kill infected cells resolving infection)
  Clearance_or_not=PAR_early$cov_CvsNC
  
  if(Clearance_or_not==1){
    m=PAR_early$m
  }else{
    m=0
  }
  
  # As different studies analyzed different samples, we accounted the difference in the volume of sample using parameter "vx" d
  Study_no=PAR_early$cov_study  # Study ID = 1 for singapore (NASAL), 2 for German (SPETUM), 3 for Korea (NASAL) , 4 for France(Europe)  (NASAL)
  
  if (Study_no==2){
    vx=10^PAR_early$vx  # German patients (Study_no=2) had their spetum analyzed and we assume that the vol of spetum is vx times the volume of nasal swabs
  }
  else{
    vx=1  # nasal swabs were employed in studies 1,3 and 4, and we use the volume of their sample to be 1
  }
  
  
  #####  Initial conditions for Viral dynamics model
  S_0 = 1e7
  I_0 = 1
  V_0 = p*I_0/c
  M1_0 = 1
  M2_0 = 0
  M3_0 = 0
  E_0 = 0
  
  # Solving ODE
  names(S_0)=""
  names(I_0)=""
  init.x <- c(S=S_0,I=I_0,V=V_0,M1=M1_0, M2=M2_0, E=E_0)
  t.out <- seq(0,30,by=0.01)  
  
  params=c()
  # Solving Viral Dynamics ODE MODEL  
  
  out <- as.data.frame(lsodar(init.x,t.out,SARS_COV_model,params))
  
  
  # Plotting
  I = out$I
  E = out$E
  
  pkI = delta*I*I^k ## Innate immune response
  pkE = m*I*E^r/(E^r+E50^r)  ## Late (acquired) immune response
  
  lines(out$time,pkE,col="aquamarine3",
          lwd=4)

  
  
  ijk=ijk+1
}


