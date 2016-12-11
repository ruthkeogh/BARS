
#log likelihood function

loglik=function(params){
  gamma0=params[1]
  gamma1=params[2]
  
  loglog.survfunc=gamma0+gamma1*log(t)
  
  
  survfunc=exp(-exp(loglog.survfunc))
  
  hazfunc=gamma1*(1/t)*exp(gamma0+gamma1*log(t)) 
  
  #negative loglik
  -sum((1-d)*log(survfunc)+d*log(hazfunc*survfunc))
}

optim.fit.0knot=optim(flexmod.0knot$coefficients,loglik,method="BFGS",hessian=F)

myhess.0knot=fdHess(flexmod.0knot$coefficients,loglik)

varcov.myhess.0knot=solve(myhess.0knot$Hessian)
