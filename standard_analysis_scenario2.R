
#---------------------
#libraries
#---------------------

library(flexsurv)

#---------------------
#set working directory
#---------------------

setwd("H:/CF_UK_lifeexp/BARS_sim")

#---------------------
#set number of simulations
#---------------------

nsim=1000

#---------------------
#set simulation scenario (1,2,3,4)
#---------------------

scenario=2

#---------------------
#storage of simulation results
#---------------------

#save aic from 11 models
flexmod.aic=matrix(nrow=nsim,ncol=11)

#save knot locations from 10 models
flexmod.knots=array(dim=c(nsim,10,12))

#save parameters estimates from 11 models
#NOTE: the maximum number of parameters is 12 (from flexmod.10knot)
flexmod.params=array(dim=c(nsim,11,12))

#save covariance matrix from 11 models
flexmod.varcov=array(dim=c(nsim,11,12,12))

#---------------------------------------
#generate random seeds
#---------------------------------------
set.seed=3010
seeds=sample(1:1e+7,size=nsim,replace=F)

#---------------------------------------
#---------------------------------------
#---------------------------------------
#start simulation loop
#---------------------------------------
#---------------------------------------
#---------------------------------------

for(sim in 1:nsim){
  
  print(sim)
  
  set.seed=seeds[sim]

  #---------------------
  #simulate data
  #---------------------
  source("./sim_data.R")
  
  #---------------------
  #fit flexible parametric survival models 
  #with pre-defined number of internal knots (nknot) and knot locations
  #---------------------
  
  for(nknot in 0:10){
    eval(parse(text=paste0("flexmod.",nknot,"knot=flexsurvspline(Surv(t,d)~1,k=nknot,scale=\"hazard\", timescale=\"log\",method=\"BFGS\")")))
  }

  #---------------------
  #save AICs from each model
  #---------------------
  
  flexmod.aic[sim,11]=flexmod.0knot$AIC
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.aic[sim,nknot]=flexmod.",nknot,"knot$AIC")))
  }
  
  #---------------------
  #save knots from each model
  #---------------------
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.knots[sim,nknot,1:(nknot+2)]=flexmod.",nknot,"knot$knots")))
  }
  
  #---------------------
  #save parameters from each model
  #---------------------
  
  flexmod.params[sim,11,1:2]=flexmod.0knot$coefficients
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.params[sim,nknot,1:(nknot+2)]=flexmod.",nknot,"knot$coefficients")))
  }
  
  #---------------------
  #save covariance matrix from each model
  #---------------------
  
  flexmod.varcov[sim,11,1:2,1:2]=solve(flexmod.0knot$opt$hessian)
  #flexmod.varcov[sim,11,1:2,1:2]=flexmod.0knot$cov
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.varcov[sim,nknot,1:(nknot+2),1:(nknot+2)]=solve(flexmod.",nknot,"knot$opt$hessian)")))
    #eval(parse(text=paste0("flexmod.varcov[sim,nknot,1:(nknot+2),1:(nknot+2)]=flexmod.",nknot,"knot$cov")))
  }

}

save(flexmod.aic,file="./sim_results/flexmod_aic_scenario2.RData")

save(flexmod.knots,file="./sim_results/flexmod_knots_scenario2.RData")

save(flexmod.params,file="./sim_results/flexmod_params_scenario2.RData")

save(flexmod.varcov,file="./sim_results/flexmod_varcov_scenario2.RData")


