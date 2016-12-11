
#---------------------
#libraries
#---------------------

library(flexsurv)
library(nlme)

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

scenario=1

#---------------------
#storage of simulation results
#---------------------

#save aic from 11 models
flexmod.aic=matrix(nrow=nsim,ncol=11)
optim.aic=matrix(nrow=nsim,ncol=11)

#save knot locations from models with 1-10 internal knots
#the knots from flexsurv and optim should be the same - this is a check
flexmod.knots=array(dim=c(nsim,10,12))
optim.knots=array(dim=c(nsim,10,12))

#save parameters estimates from 11 models
#NOTE: the maximum number of parameters is 12 (from flexmod.10knot)
flexmod.params=array(dim=c(nsim,11,12))
optim.params=array(dim=c(nsim,11,12))

#save covariance matrix from 11 models
flexmod.varcov=array(dim=c(nsim,11,12,12))
optim.varcov=array(dim=c(nsim,11,12,12))

#NOTE: results from the model with 0 internal knots (i.e. the weibull model) are stored in position 11. 
#Results from models with 1 to 10 internal knots are stored in positions 1-10

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
    eval(parse(text=paste0("source(\"optim_",nknot,"knots.R\")")))
    }
  
  #---------------------
  #save AICs from each model
  #---------------------
  
  flexmod.aic[sim,11]=flexmod.0knot$AIC
  optim.aic[sim,11]=flexmod.0knot$AIC
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.aic[sim,nknot]=flexmod.",nknot,"knot$AIC")))
    eval(parse(text=paste0("optim.aic[sim,nknot]=2*optim.fit.",nknot,"knot$value+2*(nknot+2)")))
  }
  
  #---------------------
  #save knots from each model
  #---------------------
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.knots[sim,nknot,1:(nknot+2)]=flexmod.",nknot,"knot$knots")))
    eval(parse(text=paste0("optim.knots[sim,nknot,1:(nknot+2)]=optim.knots.",nknot,"knot")))
  }
  
  #---------------------
  #save parameters from each model
  #---------------------
  
  flexmod.params[sim,11,1:2]=flexmod.0knot$coefficients
  optim.params[sim,11,1:2]=optim.0knot$coefficients
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.params[sim,nknot,1:(nknot+2)]=flexmod.",nknot,"knot$coefficients")))
    eval(parse(text=paste0("optim.params[sim,nknot,1:(nknot+2)]=optim.fit.",nknot,"knot$par")))
  }
  
  #---------------------
  #save covariance matrix from each model
  #---------------------
  
  flexmod.varcov[sim,11,1:2,1:2]=solve(flexmod.0knot$opt$hessian)
  #flexmod.varcov[sim,11,1:2,1:2]=flexmod.0knot$cov
  optim.varcov[sim,11,1:2,1:2]=varcov.myhess.0knot
  
  for(nknot in 1:10){
    eval(parse(text=paste0("flexmod.varcov[sim,nknot,1:(nknot+2),1:(nknot+2)]=solve(flexmod.",nknot,"knot$opt$hessian)")))
    #eval(parse(text=paste0("flexmod.varcov[sim,nknot,1:(nknot+2),1:(nknot+2)]=flexmod.",nknot,"knot$cov")))
    eval(parse(text=paste0("optim.varcov[sim,nknot,1:(nknot+2),1:(nknot+2)]=varcov.myhess.",nknot,"knot")))
  }

}

save(flexmod.aic,file="./sim_results/flexmod_aic_scenario1.RData")
save(optim.aic,file="./sim_results/optim_aic_scenario1.RData")

save(flexmod.knots,file="./sim_results/flexmod_knots_scenario1.RData")
save(optim.knots,file="./sim_results/optim_knots_scenario1.RData")

save(flexmod.params,file="./sim_results/flexmod_params_scenario1.RData")
save(optim.params,file="./sim_results/optim_params_scenario1.RData")

save(flexmod.varcov,file="./sim_results/flexmod_varcov_scenario1.RData")
save(optim.varcov,file="./sim_results/optim_varcov_scenario1.RData")


