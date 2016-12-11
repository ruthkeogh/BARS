

#knot positions
optim.knots.1knot=quantile(log(t[d==1]),c(seq(0,1-1/2,1/2),1))
optim.myknots=optim.knots.1knot #I rename here just to save myself some copy and pasting in the loglik function below. 

maxknot=3 #total number of knots (including boundary knots)

#spline basis function

v=function(t,kmin,kmax,k){
  (pmax(log(t)-k,0))^3-((kmax-k)/(kmax-kmin))*(pmax(log(t)-kmin,0))^3-
    ((k-kmin)/(kmax-kmin))*(pmax(log(t)-kmax,0))^3
}

#derivative of spline basis function (used in the hazard function)

v.deriv=function(t,kmin,kmax,k){
  (3/t)*(pmax((log(t)-k),0))^2-((kmax-k)/(kmax-kmin))*(3/t)*(pmax((log(t)-kmin),0))^2-
    ((k-kmin)/(kmax-kmin))*(3/t)*(pmax((log(t)-kmax),0))^2
}

#log likelihood function

loglik=function(params){
  gamma0=params[1]
  gamma1=params[2]
  gamma2=params[3]
  
  loglog.survfunc=gamma0+gamma1*log(t)+
    gamma2*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[2]))
  
  
  survfunc=exp(-exp(loglog.survfunc))
  
  hazfunc=(gamma1*(1/t)+
             gamma2*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[2])))*
    exp(gamma0+gamma1*log(t)+
          gamma2*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[2]))) 
  
  #negative loglik
  -sum((1-d)*log(survfunc)+d*log(hazfunc*survfunc))
}

optim.fit.1knot=optim(flexmod.1knot$coefficients,loglik,method="BFGS",hessian=F)

myhess.1knot=fdHess(flexmod.1knot$coefficients,loglik)

varcov.myhess.1knot=solve(myhess.1knot$Hessian)
