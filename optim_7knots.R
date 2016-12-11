

#knot positions
optim.knots.7knot=quantile(log(t[d==1]),c(seq(0,1-1/8,1/8),1))
optim.myknots=optim.knots.7knot #I rename here just to save myself some copy and pasting in the loglik function below. 

maxknot=9 #total number of knots (including boundary knots)

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
  gamma3=params[4]
  gamma4=params[5]
  gamma5=params[6]
  gamma6=params[7]
  gamma7=params[8]
  gamma8=params[9]
  
  loglog.survfunc=gamma0+gamma1*log(t)+
    gamma2*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[2]))+
    gamma3*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[3]))+
    gamma4*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[4]))+
    gamma5*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[5]))+
    gamma6*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[6]))+
    gamma7*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[7]))+
    gamma8*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[8]))
  
  
  survfunc=exp(-exp(loglog.survfunc))
  
  hazfunc=(gamma1*(1/t)+
             gamma2*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[2]))+
             gamma3*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[3]))+
             gamma4*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[4]))+
             gamma5*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[5]))+
             gamma6*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[6]))+
             gamma7*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[7]))+
             gamma8*sapply(t,FUN=function(t)v.deriv(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[8])))*
    exp(gamma0+gamma1*log(t)+
          gamma2*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[2]))+
          gamma3*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[3]))+
          gamma4*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[4]))+
          gamma5*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[5]))+
          gamma6*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[6]))+
          gamma7*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[7]))+
          gamma8*sapply(t,FUN=function(t)v(t,optim.myknots[1],optim.myknots[maxknot],optim.myknots[8]))) 
  
  #negative loglik
  -sum((1-d)*log(survfunc)+d*log(hazfunc*survfunc))
}

optim.fit.7knot=optim(flexmod.7knot$coefficients,loglik,method="BFGS",hessian=F)

myhess.7knot=fdHess(flexmod.7knot$coefficients,loglik)

varcov.myhess.7knot=solve(myhess.7knot$Hessian)
