#-------------------------
#working directory
#-------------------------

setwd("H:/CF_UK_lifeexp/BARS_sim")

#-------------------------
#libraries
#-------------------------

library(xtable)
library(ggplot2)
library(gridExtra)
library(plyr)

#-------------------------
#load results
#-------------------------

load(file="./sim_results/flexmod_aic_scenario1.RData") #flexmod.aic
load(file="./sim_results/flexmod_knots_scenario1.RData") #flexmod.knots
load(file="./sim_results/flexmod_params_scenario1.RData") #flexmod.params
load(file="./sim_results/flexmod_varcov_scenario1.RData") #flexmod.varcov
flexmod.aic1=flexmod.aic
flexmod.knots1=flexmod.knots
flexmod.params1=flexmod.params
flexmod.varcov1=flexmod.varcov

load(file="./sim_results/flexmod_aic_scenario2.RData") #flexmod.aic
load(file="./sim_results/flexmod_knots_scenario2.RData") #flexmod.knots
load(file="./sim_results/flexmod_params_scenario2.RData") #flexmod.params
load(file="./sim_results/flexmod_varcov_scenario2.RData") #flexmod.varcov
flexmod.aic2=flexmod.aic
flexmod.knots2=flexmod.knots
flexmod.params2=flexmod.params
flexmod.varcov2=flexmod.varcov

load(file="./sim_results/flexmod_aic_scenario3.RData") #flexmod.aic
load(file="./sim_results/flexmod_knots_scenario3.RData") #flexmod.knots
load(file="./sim_results/flexmod_params_scenario3.RData") #flexmod.params
load(file="./sim_results/flexmod_varcov_scenario3.RData") #flexmod.varcov
flexmod.aic3=flexmod.aic
flexmod.knots3=flexmod.knots
flexmod.params3=flexmod.params
flexmod.varcov3=flexmod.varcov

load(file="./sim_results/flexmod_aic_scenario4.RData") #flexmod.aic
load(file="./sim_results/flexmod_knots_scenario4.RData") #flexmod.knots
load(file="./sim_results/flexmod_params_scenario4.RData") #flexmod.params
load(file="./sim_results/flexmod_varcov_scenario4.RData") #flexmod.varcov
flexmod.aic4=flexmod.aic
flexmod.knots4=flexmod.knots
flexmod.params4=flexmod.params
flexmod.varcov4=flexmod.varcov

#------------------------------------------
#true survivor and hazard functions
#this corresponds to how the data were generated under the 4 simulation scenarios
#------------------------------------------

#true survivor function: scenario 1

true.surv.func1=function(t){
  exp(-0.6*t^(0.8))
}

#true hazard function: scenario 1

true.haz.func1=function(t){
  (0.6*0.8*(t^(0.8-1))*exp(-0.6*t^(0.8)))/
    (exp(-0.6*t^(0.8)))
}

#true survivor function: scenarios 2-4

surv.func=function(p,lambda1,lambda2,gamma1,gamma2,t){
  p*exp(-lambda1*t^(gamma1))+(1-p)*exp(-lambda2*t^(gamma2))
}

true.surv.func2=function(t){surv.func(0.2,0.2,1.6,0.8,1,t)}
true.surv.func3=function(t){surv.func(0.5,1,1,1.5,0.5,t)}
true.surv.func4=function(t){surv.func(0.7,0.03,0.3,1.9,2.5,t)}

#true hazard function: scenarios 2-4

haz.func=function(p,lambda1,lambda2,gamma1,gamma2,t){
  (p*lambda1*gamma1*(t^(gamma1-1))*exp(-lambda1*t^(gamma1))+(1-p)*lambda2*gamma2*(t^(gamma2-1))*exp(-lambda2*t^(gamma2)))/
    (p*exp(-lambda1*t^(gamma1))+(1-p)*exp(-lambda2*t^(gamma2)))
}

true.haz.func2=function(t){haz.func(0.2,0.2,1.6,0.8,1,t)}
true.haz.func3=function(t){haz.func(0.5,1,1,1.5,0.5,t)}
true.haz.func4=function(t){haz.func(0.7,0.03,0.3,1.9,2.5,t)}

#------------------------------------------
#spline basis functions
#------------------------------------------

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

#------------------------------------------
#survivor functions: for flex models with 1 to 10 internal knots
#------------------------------------------

#0 knots (weibull)

surv.func.0knot=function(t,flexmod.params.1,flexmod.params.2){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)))
}

#1 knot

surv.func.1knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)))
}

#2 knots

surv.func.2knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)))
}

#3 knots

surv.func.3knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)))
}

#4 knots

surv.func.4knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)))
}

#5 knots

surv.func.5knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                         flexmod.params.7,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
             flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)))
}

#6 knots

surv.func.6knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                          flexmod.params.7,flexmod.params.8,
                          flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                          flexmod.knots.7){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
             flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
             flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)))
}

#7 knots

surv.func.7knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                          flexmod.params.7,flexmod.params.8,flexmod.params.9,
                          flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                          flexmod.knots.7,flexmod.knots.8){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
             flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
             flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
             flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)))
}

#8 knots

surv.func.8knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                          flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,
                          flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                          flexmod.knots.7,flexmod.knots.8,flexmod.knots.9){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
             flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
             flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
             flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
             flexmod.params.10*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)))
}

#9 knots

surv.func.9knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                          flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,flexmod.params.11,
                          flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                          flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,flexmod.knots.10){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
             flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
             flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
             flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
             flexmod.params.10*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)+
             flexmod.params.11*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10)))
}


#10 knots

surv.func.10knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                         flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,flexmod.params.11,flexmod.params.12,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                         flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,flexmod.knots.10,flexmod.knots.11){
  exp(-exp(flexmod.params.1+flexmod.params.2*log(t)+
             flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
             flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
             flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
             flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
             flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
             flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
             flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
             flexmod.params.10*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)+
             flexmod.params.11*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10)+
             flexmod.params.12*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.11)))
}

#------------------------------------------
#hazard functions: for flex models with 1 to 10 internal knots
#------------------------------------------

#0 knot (weibull)

haz.func.0knot=function(t,flexmod.params.1,flexmod.params.2){
  (flexmod.params.2*(1/t))*
    exp(flexmod.params.1+flexmod.params.2*log(t))  
}

#1 knot

haz.func.1knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2))  
}

#2 knot

haz.func.2knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3))  
}

#3 knot

haz.func.3knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4))  
}

#4 knot

haz.func.4knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5))  
}

#5 knot

haz.func.5knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                        flexmod.params.7,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
     flexmod.params.7*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
          flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6))  
}

#6 knot

haz.func.6knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                        flexmod.params.7,flexmod.params.8,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                        flexmod.knots.7){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
     flexmod.params.7*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
     flexmod.params.8*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
          flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
          flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7))  
}

#7 knot

haz.func.7knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                        flexmod.params.7,flexmod.params.8,flexmod.params.9,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                        flexmod.knots.7,flexmod.knots.8){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
     flexmod.params.7*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
     flexmod.params.8*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
     flexmod.params.9*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
          flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
          flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
          flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8))  
}

#8 knot

haz.func.8knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                        flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,
                        flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                        flexmod.knots.7,flexmod.knots.8,flexmod.knots.9){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
     flexmod.params.7*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
     flexmod.params.8*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
     flexmod.params.9*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
     flexmod.params.10*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
          flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
          flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
          flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
          flexmod.params.10*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9))  
}

#9 knot

haz.func.9knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                         flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,flexmod.params.11,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                         flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,flexmod.knots.10){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
     flexmod.params.7*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
     flexmod.params.8*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
     flexmod.params.9*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
     flexmod.params.10*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)+
     flexmod.params.11*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
          flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
          flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
          flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
          flexmod.params.10*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)+
          flexmod.params.11*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10))  
}

#10 knot

haz.func.10knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                         flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,flexmod.params.11,flexmod.params.12,
                         flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                         flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,flexmod.knots.10,flexmod.knots.11){
  (flexmod.params.2*(1/t)+
     flexmod.params.3*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
     flexmod.params.4*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
     flexmod.params.5*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
     flexmod.params.6*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
     flexmod.params.7*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
     flexmod.params.8*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
     flexmod.params.9*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
     flexmod.params.10*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)+
     flexmod.params.11*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10)+
     flexmod.params.12*v.deriv(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.11))*
    exp(flexmod.params.1+flexmod.params.2*log(t)+
          flexmod.params.3*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)+
          flexmod.params.4*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)+
          flexmod.params.5*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)+
          flexmod.params.6*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)+
          flexmod.params.7*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)+
          flexmod.params.8*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)+
          flexmod.params.9*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)+
          flexmod.params.10*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)+
          flexmod.params.11*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10)+
          flexmod.params.12*v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.11))  
}

#------------------------------------------
#survivor functions as function of 't' and 'sim' for each simulation scenario
#------------------------------------------

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.0knot.scenario",scen,"=function(t,sim){
                         surv.func.0knot(t,flexmod.params",scen,"[sim,11,1],flexmod.params",scen,"[sim,11,2])
}")))
}

#surv.func.1knot.scenario1=function(t,sim){
#  surv.func.1knot(t,flexmod.params1[sim,1,1],flexmod.params1[sim,1,2],flexmod.params1[sim,1,3],
#                  flexmod.knots1[sim,1,1],flexmod.knots1[sim,1,max.knot],flexmod.knots1[sim,1,2])

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.1knot.scenario",scen,"=function(t,sim){
                         surv.func.1knot(t,flexmod.params",scen,"[sim,1,1],flexmod.params",scen,"[sim,1,2],flexmod.params",scen,"[sim,1,3],
                         flexmod.knots",scen,"[sim,1,1],flexmod.knots",scen,"[sim,1,3],flexmod.knots",scen,"[sim,1,2])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.2knot.scenario",scen,"=function(t,sim){
                         surv.func.2knot(t,flexmod.params",scen,"[sim,2,1],flexmod.params",scen,"[sim,2,2],flexmod.params",scen,"[sim,2,3],flexmod.params",scen,"[sim,2,4],
                         flexmod.knots",scen,"[sim,2,1],flexmod.knots",scen,"[sim,2,4],flexmod.knots",scen,"[sim,2,2],flexmod.knots",scen,"[sim,2,3])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.3knot.scenario",scen,"=function(t,sim){
                         surv.func.3knot(t,flexmod.params",scen,"[sim,3,1],flexmod.params",scen,"[sim,3,2],flexmod.params",scen,"[sim,3,3],flexmod.params",scen,"[sim,3,4],
                         flexmod.params",scen,"[sim,3,5],
                         flexmod.knots",scen,"[sim,3,1],flexmod.knots",scen,"[sim,3,5],flexmod.knots",scen,"[sim,3,2],flexmod.knots",scen,"[sim,3,3],
                         flexmod.knots",scen,"[sim,3,4])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.4knot.scenario",scen,"=function(t,sim){
                         surv.func.4knot(t,flexmod.params",scen,"[sim,4,1],flexmod.params",scen,"[sim,4,2],flexmod.params",scen,"[sim,4,3],flexmod.params",scen,"[sim,4,4],
                         flexmod.params",scen,"[sim,4,5],flexmod.params",scen,"[sim,4,6],
                         flexmod.knots",scen,"[sim,4,1],flexmod.knots",scen,"[sim,4,6],flexmod.knots",scen,"[sim,4,2],flexmod.knots",scen,"[sim,4,3],
                         flexmod.knots",scen,"[sim,4,4],flexmod.knots",scen,"[sim,4,5])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.5knot.scenario",scen,"=function(t,sim){
                         surv.func.5knot(t,flexmod.params",scen,"[sim,5,1],flexmod.params",scen,"[sim,5,2],flexmod.params",scen,"[sim,5,3],flexmod.params",scen,"[sim,5,4],
                         flexmod.params",scen,"[sim,5,5],flexmod.params",scen,"[sim,5,6],flexmod.params",scen,"[sim,5,7],
                         flexmod.knots",scen,"[sim,5,1],flexmod.knots",scen,"[sim,5,7],flexmod.knots",scen,"[sim,5,2],flexmod.knots",scen,"[sim,5,3],
                         flexmod.knots",scen,"[sim,5,4],flexmod.knots",scen,"[sim,5,5],flexmod.knots",scen,"[sim,5,6])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.6knot.scenario",scen,"=function(t,sim){
                         surv.func.6knot(t,flexmod.params",scen,"[sim,6,1],flexmod.params",scen,"[sim,6,2],flexmod.params",scen,"[sim,6,3],flexmod.params",scen,"[sim,6,4],
                         flexmod.params",scen,"[sim,6,5],flexmod.params",scen,"[sim,6,6],flexmod.params",scen,"[sim,6,7],flexmod.params",scen,"[sim,6,8],
                         flexmod.knots",scen,"[sim,6,1],flexmod.knots",scen,"[sim,6,8],flexmod.knots",scen,"[sim,6,2],flexmod.knots",scen,"[sim,6,3],
                         flexmod.knots",scen,"[sim,6,4],flexmod.knots",scen,"[sim,6,5],flexmod.knots",scen,"[sim,6,6],flexmod.knots",scen,"[sim,6,7])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.7knot.scenario",scen,"=function(t,sim){
                         surv.func.7knot(t,flexmod.params",scen,"[sim,7,1],flexmod.params",scen,"[sim,7,2],flexmod.params",scen,"[sim,7,3],flexmod.params",scen,"[sim,7,4],
                         flexmod.params",scen,"[sim,7,5],flexmod.params",scen,"[sim,7,6],flexmod.params",scen,"[sim,7,7],flexmod.params",scen,"[sim,7,8],
                         flexmod.params",scen,"[sim,7,9],
                         flexmod.knots",scen,"[sim,7,1],flexmod.knots",scen,"[sim,7,9],flexmod.knots",scen,"[sim,7,2],flexmod.knots",scen,"[sim,7,3],
                         flexmod.knots",scen,"[sim,7,4],flexmod.knots",scen,"[sim,7,5],flexmod.knots",scen,"[sim,7,6],flexmod.knots",scen,"[sim,7,7],
                         flexmod.knots",scen,"[sim,7,8])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.8knot.scenario",scen,"=function(t,sim){
                         surv.func.8knot(t,flexmod.params",scen,"[sim,8,1],flexmod.params",scen,"[sim,8,2],flexmod.params",scen,"[sim,8,3],flexmod.params",scen,"[sim,8,4],
                         flexmod.params",scen,"[sim,8,5],flexmod.params",scen,"[sim,8,6],flexmod.params",scen,"[sim,8,7],flexmod.params",scen,"[sim,8,8],
                         flexmod.params",scen,"[sim,8,9],flexmod.params",scen,"[sim,8,10],
                         flexmod.knots",scen,"[sim,8,1],flexmod.knots",scen,"[sim,8,10],flexmod.knots",scen,"[sim,8,2],flexmod.knots",scen,"[sim,8,3],
                         flexmod.knots",scen,"[sim,8,4],flexmod.knots",scen,"[sim,8,5],flexmod.knots",scen,"[sim,8,6],flexmod.knots",scen,"[sim,8,7],
                         flexmod.knots",scen,"[sim,8,8],flexmod.knots",scen,"[sim,8,9])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.9knot.scenario",scen,"=function(t,sim){
                         surv.func.9knot(t,flexmod.params",scen,"[sim,9,1],flexmod.params",scen,"[sim,9,2],flexmod.params",scen,"[sim,9,3],flexmod.params",scen,"[sim,9,4],
                         flexmod.params",scen,"[sim,9,5],flexmod.params",scen,"[sim,9,6],flexmod.params",scen,"[sim,9,7],flexmod.params",scen,"[sim,9,8],
                         flexmod.params",scen,"[sim,9,9],flexmod.params",scen,"[sim,9,10],flexmod.params",scen,"[sim,9,11],
                         flexmod.knots",scen,"[sim,9,1],flexmod.knots",scen,"[sim,9,11],flexmod.knots",scen,"[sim,9,2],flexmod.knots",scen,"[sim,9,3],
                         flexmod.knots",scen,"[sim,9,4],flexmod.knots",scen,"[sim,9,5],flexmod.knots",scen,"[sim,9,6],flexmod.knots",scen,"[sim,9,7],
                         flexmod.knots",scen,"[sim,9,8],flexmod.knots",scen,"[sim,9,9],flexmod.knots",scen,"[sim,9,10])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.func.10knot.scenario",scen,"=function(t,sim){
                         surv.func.10knot(t,flexmod.params",scen,"[sim,10,1],flexmod.params",scen,"[sim,10,2],flexmod.params",scen,"[sim,10,3],flexmod.params",scen,"[sim,10,4],
                         flexmod.params",scen,"[sim,10,5],flexmod.params",scen,"[sim,10,6],flexmod.params",scen,"[sim,10,7],flexmod.params",scen,"[sim,10,8],
                         flexmod.params",scen,"[sim,10,9],flexmod.params",scen,"[sim,10,10],flexmod.params",scen,"[sim,10,11],flexmod.params",scen,"[sim,10,12],
                         flexmod.knots",scen,"[sim,10,1],flexmod.knots",scen,"[sim,10,12],flexmod.knots",scen,"[sim,10,2],flexmod.knots",scen,"[sim,10,3],
                         flexmod.knots",scen,"[sim,10,4],flexmod.knots",scen,"[sim,10,5],flexmod.knots",scen,"[sim,10,6],flexmod.knots",scen,"[sim,10,7],
                         flexmod.knots",scen,"[sim,10,8],flexmod.knots",scen,"[sim,10,9],flexmod.knots",scen,"[sim,10,10],flexmod.knots",scen,"[sim,10,11])
}")))
}

#------------------------------------------
#hazard functions as function of 't' and 'sim' for each simulation scenario
#------------------------------------------

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.0knot.scenario",scen,"=function(t,sim){
                         haz.func.0knot(t,flexmod.params",scen,"[sim,11,1],flexmod.params",scen,"[sim,11,2])
}")))
}


for(scen in 1:4){
  eval(parse(text=paste0("haz.func.1knot.scenario",scen,"=function(t,sim){
                         haz.func.1knot(t,flexmod.params",scen,"[sim,1,1],flexmod.params",scen,"[sim,1,2],flexmod.params",scen,"[sim,1,3],
                         flexmod.knots",scen,"[sim,1,1],flexmod.knots",scen,"[sim,1,3],flexmod.knots",scen,"[sim,1,2])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.2knot.scenario",scen,"=function(t,sim){
                         haz.func.2knot(t,flexmod.params",scen,"[sim,2,1],flexmod.params",scen,"[sim,2,2],flexmod.params",scen,"[sim,2,3],flexmod.params",scen,"[sim,2,4],
                         flexmod.knots",scen,"[sim,2,1],flexmod.knots",scen,"[sim,2,4],flexmod.knots",scen,"[sim,2,2],flexmod.knots",scen,"[sim,2,3])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.3knot.scenario",scen,"=function(t,sim){
                         haz.func.3knot(t,flexmod.params",scen,"[sim,3,1],flexmod.params",scen,"[sim,3,2],flexmod.params",scen,"[sim,3,3],flexmod.params",scen,"[sim,3,4],
                         flexmod.params",scen,"[sim,3,5],
                         flexmod.knots",scen,"[sim,3,1],flexmod.knots",scen,"[sim,3,5],flexmod.knots",scen,"[sim,3,2],flexmod.knots",scen,"[sim,3,3],
                         flexmod.knots",scen,"[sim,3,4])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.4knot.scenario",scen,"=function(t,sim){
                         haz.func.4knot(t,flexmod.params",scen,"[sim,4,1],flexmod.params",scen,"[sim,4,2],flexmod.params",scen,"[sim,4,3],flexmod.params",scen,"[sim,4,4],
                         flexmod.params",scen,"[sim,4,5],flexmod.params",scen,"[sim,4,6],
                         flexmod.knots",scen,"[sim,4,1],flexmod.knots",scen,"[sim,4,6],flexmod.knots",scen,"[sim,4,2],flexmod.knots",scen,"[sim,4,3],
                         flexmod.knots",scen,"[sim,4,4],flexmod.knots",scen,"[sim,4,5])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.5knot.scenario",scen,"=function(t,sim){
                         haz.func.5knot(t,flexmod.params",scen,"[sim,5,1],flexmod.params",scen,"[sim,5,2],flexmod.params",scen,"[sim,5,3],flexmod.params",scen,"[sim,5,4],
                         flexmod.params",scen,"[sim,5,5],flexmod.params",scen,"[sim,5,6],flexmod.params",scen,"[sim,5,7],
                         flexmod.knots",scen,"[sim,5,1],flexmod.knots",scen,"[sim,5,7],flexmod.knots",scen,"[sim,5,2],flexmod.knots",scen,"[sim,5,3],
                         flexmod.knots",scen,"[sim,5,4],flexmod.knots",scen,"[sim,5,5],flexmod.knots",scen,"[sim,5,6])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.6knot.scenario",scen,"=function(t,sim){
                         haz.func.6knot(t,flexmod.params",scen,"[sim,6,1],flexmod.params",scen,"[sim,6,2],flexmod.params",scen,"[sim,6,3],flexmod.params",scen,"[sim,6,4],
                         flexmod.params",scen,"[sim,6,5],flexmod.params",scen,"[sim,6,6],flexmod.params",scen,"[sim,6,7],flexmod.params",scen,"[sim,6,8],
                         flexmod.knots",scen,"[sim,6,1],flexmod.knots",scen,"[sim,6,8],flexmod.knots",scen,"[sim,6,2],flexmod.knots",scen,"[sim,6,3],
                         flexmod.knots",scen,"[sim,6,4],flexmod.knots",scen,"[sim,6,5],flexmod.knots",scen,"[sim,6,6],flexmod.knots",scen,"[sim,6,7])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.7knot.scenario",scen,"=function(t,sim){
                         haz.func.7knot(t,flexmod.params",scen,"[sim,7,1],flexmod.params",scen,"[sim,7,2],flexmod.params",scen,"[sim,7,3],flexmod.params",scen,"[sim,7,4],
                         flexmod.params",scen,"[sim,7,5],flexmod.params",scen,"[sim,7,6],flexmod.params",scen,"[sim,7,7],flexmod.params",scen,"[sim,7,8],
                         flexmod.params",scen,"[sim,7,9],
                         flexmod.knots",scen,"[sim,7,1],flexmod.knots",scen,"[sim,7,9],flexmod.knots",scen,"[sim,7,2],flexmod.knots",scen,"[sim,7,3],
                         flexmod.knots",scen,"[sim,7,4],flexmod.knots",scen,"[sim,7,5],flexmod.knots",scen,"[sim,7,6],flexmod.knots",scen,"[sim,7,7],
                         flexmod.knots",scen,"[sim,7,8])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.8knot.scenario",scen,"=function(t,sim){
                         haz.func.8knot(t,flexmod.params",scen,"[sim,8,1],flexmod.params",scen,"[sim,8,2],flexmod.params",scen,"[sim,8,3],flexmod.params",scen,"[sim,8,4],
                         flexmod.params",scen,"[sim,8,5],flexmod.params",scen,"[sim,8,6],flexmod.params",scen,"[sim,8,7],flexmod.params",scen,"[sim,8,8],
                         flexmod.params",scen,"[sim,8,9],flexmod.params",scen,"[sim,8,10],
                         flexmod.knots",scen,"[sim,8,1],flexmod.knots",scen,"[sim,8,10],flexmod.knots",scen,"[sim,8,2],flexmod.knots",scen,"[sim,8,3],
                         flexmod.knots",scen,"[sim,8,4],flexmod.knots",scen,"[sim,8,5],flexmod.knots",scen,"[sim,8,6],flexmod.knots",scen,"[sim,8,7],
                         flexmod.knots",scen,"[sim,8,8],flexmod.knots",scen,"[sim,8,9])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.9knot.scenario",scen,"=function(t,sim){
                         haz.func.9knot(t,flexmod.params",scen,"[sim,9,1],flexmod.params",scen,"[sim,9,2],flexmod.params",scen,"[sim,9,3],flexmod.params",scen,"[sim,9,4],
                         flexmod.params",scen,"[sim,9,5],flexmod.params",scen,"[sim,9,6],flexmod.params",scen,"[sim,9,7],flexmod.params",scen,"[sim,9,8],
                         flexmod.params",scen,"[sim,9,9],flexmod.params",scen,"[sim,9,10],flexmod.params",scen,"[sim,9,11],
                         flexmod.knots",scen,"[sim,9,1],flexmod.knots",scen,"[sim,9,11],flexmod.knots",scen,"[sim,9,2],flexmod.knots",scen,"[sim,9,3],
                         flexmod.knots",scen,"[sim,9,4],flexmod.knots",scen,"[sim,9,5],flexmod.knots",scen,"[sim,9,6],flexmod.knots",scen,"[sim,9,7],
                         flexmod.knots",scen,"[sim,9,8],flexmod.knots",scen,"[sim,9,9],flexmod.knots",scen,"[sim,9,10])
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("haz.func.10knot.scenario",scen,"=function(t,sim){
                         haz.func.10knot(t,flexmod.params",scen,"[sim,10,1],flexmod.params",scen,"[sim,10,2],flexmod.params",scen,"[sim,10,3],flexmod.params",scen,"[sim,10,4],
                         flexmod.params",scen,"[sim,10,5],flexmod.params",scen,"[sim,10,6],flexmod.params",scen,"[sim,10,7],flexmod.params",scen,"[sim,10,8],
                         flexmod.params",scen,"[sim,10,9],flexmod.params",scen,"[sim,10,10],flexmod.params",scen,"[sim,10,11],flexmod.params",scen,"[sim,10,12],
                         flexmod.knots",scen,"[sim,10,1],flexmod.knots",scen,"[sim,10,12],flexmod.knots",scen,"[sim,10,2],flexmod.knots",scen,"[sim,10,3],
                         flexmod.knots",scen,"[sim,10,4],flexmod.knots",scen,"[sim,10,5],flexmod.knots",scen,"[sim,10,6],flexmod.knots",scen,"[sim,10,7],
                         flexmod.knots",scen,"[sim,10,8],flexmod.knots",scen,"[sim,10,9],flexmod.knots",scen,"[sim,10,10],flexmod.knots",scen,"[sim,10,11])
}")))
}

#-------------------------
#absolute and squared difference between flex survival curves and true survival curves
#for each flexsurv model 
#-------------------------

#grid of time points in 0:10
t.grid=seq(0.1,10,0.1)

#for 0 knots (weibull): we use 'knot11' instead of 'knot0' just to make it easier a bit frther below this to select the best fitting model
for(scen in 1:4){
  eval(parse(text=paste0("diffsurv.abs.knot11.scenario",scen,"=rep(0,1000)")))
  eval(parse(text=paste0("diffsurv.sq.knot11.scenario",scen,"=rep(0,1000)")))
  for(sim in 1:1000){
    eval(parse(text=paste0("diffsurv.abs.knot11.scenario",scen,"[sim]=sum(0.1*abs(sapply(t.grid,function(t){surv.func.0knot.scenario",scen,"(t,sim)})-sapply(t.grid,function(t){true.surv.func",scen,"(t)})))")))
    eval(parse(text=paste0("diffsurv.sq.knot11.scenario",scen,"[sim]=sum((0.1*sapply(t.grid,function(t){surv.func.0knot.scenario",scen,"(t,sim)})-0.1*sapply(t.grid,function(t){true.surv.func",scen,"(t)}))^2)")))
  }
}

#for 1:10 knots
for(k in 1:10){  
  for(scen in 1:4){
    eval(parse(text=paste0("diffsurv.abs.knot",k,".scenario",scen,"=rep(0,1000)")))
    eval(parse(text=paste0("diffsurv.sq.knot",k,".scenario",scen,"=rep(0,1000)")))
    for(sim in 1:1000){
      eval(parse(text=paste0("diffsurv.abs.knot",k,".scenario",scen,"[sim]=sum(0.1*abs(sapply(t.grid,function(t){surv.func.",k,"knot.scenario",scen,"(t,sim)})-sapply(t.grid,function(t){true.surv.func",scen,"(t)})))")))
      eval(parse(text=paste0("diffsurv.sq.knot",k,".scenario",scen,"[sim]=sum((0.1*sapply(t.grid,function(t){surv.func.",k,"knot.scenario",scen,"(t,sim)})-0.1*sapply(t.grid,function(t){true.surv.func",scen,"(t)}))^2)")))
    }
  }
}

#-------------------------
#absolute and squared difference between flex hazard curves and true hazard curves
#for each flexsurv model 
#-------------------------

#grid of time points in 0:10
t.grid=seq(0.1,10,0.1)

#for 0 knots (weibull)
for(k in 1:10){  
  for(scen in 1:4){
    eval(parse(text=paste0("diffhaz.abs.knot11.scenario",scen,"=rep(0,1000)")))
    eval(parse(text=paste0("diffhaz.sq.knot11.scenario",scen,"=rep(0,1000)")))
    for(sim in 1:1000){
      eval(parse(text=paste0("diffhaz.abs.knot11.scenario",scen,"[sim]=sum(0.1*abs(sapply(t.grid,function(t){haz.func.0knot.scenario",scen,"(t,sim)})-sapply(t.grid,function(t){true.haz.func",scen,"(t)})))")))
      eval(parse(text=paste0("diffhaz.sq.knot11.scenario",scen,"[sim]=sum((0.1*sapply(t.grid,function(t){haz.func.0knot.scenario",scen,"(t,sim)})-0.1*sapply(t.grid,function(t){true.haz.func",scen,"(t)}))^2)")))
    }
  }
}

#for 1:10 knots
for(k in 1:10){  
  for(scen in 1:4){
    eval(parse(text=paste0("diffhaz.abs.knot",k,".scenario",scen,"=rep(0,1000)")))
    eval(parse(text=paste0("diffhaz.sq.knot",k,".scenario",scen,"=rep(0,1000)")))
    for(sim in 1:1000){
      eval(parse(text=paste0("diffhaz.abs.knot",k,".scenario",scen,"[sim]=sum(0.1*abs(sapply(t.grid,function(t){haz.func.",k,"knot.scenario",scen,"(t,sim)})-sapply(t.grid,function(t){true.haz.func",scen,"(t)})))")))
      eval(parse(text=paste0("diffhaz.sq.knot",k,".scenario",scen,"[sim]=sum((0.1*sapply(t.grid,function(t){haz.func.",k,"knot.scenario",scen,"(t,sim)})-0.1*sapply(t.grid,function(t){true.haz.func",scen,"(t)}))^2)")))
    }
  }
}


#-------------------------
#select best fitting flexsurv model in each simulated data set, using AIC
#save the difference between the curves from the best fitting model and the true model
#-------------------------

for(scen in 1:4){
  eval(parse(text=paste0("bestfit.scenario",scen,"=rep(0,1000)")))
  eval(parse(text=paste0("diffsurv.abs.bestfit.scenario",scen,"=rep(0,1000)")))  
  eval(parse(text=paste0("diffsurv.sq.bestfit.scenario",scen,"=rep(0,1000)")))
  eval(parse(text=paste0("diffhaz.abs.bestfit.scenario",scen,"=rep(0,1000)")))  
  eval(parse(text=paste0("diffhaz.sq.bestfit.scenario",scen,"=rep(0,1000)")))
  for(sim in 1:1000){
    eval(parse(text=paste0("bestfit.scenario",scen,"[sim]=which.min(flexmod.aic",scen,"[sim,])")))
    eval(parse(text=paste0("bestfit=bestfit.scenario",scen,"[sim]")))
    eval(parse(text=paste0("diffsurv.abs.bestfit.scenario",scen,"[sim]=diffsurv.abs.knot",bestfit,".scenario",scen,"[sim]")))
    eval(parse(text=paste0("diffsurv.sq.bestfit.scenario",scen,"[sim]=diffsurv.sq.knot",bestfit,".scenario",scen,"[sim]")))
    eval(parse(text=paste0("diffhaz.abs.bestfit.scenario",scen,"[sim]=diffhaz.abs.knot",bestfit,".scenario",scen,"[sim]")))
    eval(parse(text=paste0("diffhaz.sq.bestfit.scenario",scen,"[sim]=diffhaz.sq.knot",bestfit,".scenario",scen,"[sim]")))
  }
}

#-------------------------
#table of differences between survival curves
#-------------------------

survdiff.sq=rbind(c(mean(diffsurv.sq.knot11.scenario1),mean(diffsurv.sq.knot11.scenario2),mean(diffsurv.sq.knot11.scenario3),mean(diffsurv.sq.knot11.scenario4)),
c(mean(diffsurv.sq.knot1.scenario1),mean(diffsurv.sq.knot1.scenario2),mean(diffsurv.sq.knot1.scenario3),mean(diffsurv.sq.knot1.scenario4)),
c(mean(diffsurv.sq.knot2.scenario1),mean(diffsurv.sq.knot2.scenario2),mean(diffsurv.sq.knot2.scenario3),mean(diffsurv.sq.knot2.scenario4)),
c(mean(diffsurv.sq.knot3.scenario1),mean(diffsurv.sq.knot3.scenario2),mean(diffsurv.sq.knot3.scenario3),mean(diffsurv.sq.knot3.scenario4)),
c(mean(diffsurv.sq.knot4.scenario1),mean(diffsurv.sq.knot4.scenario2),mean(diffsurv.sq.knot4.scenario3),mean(diffsurv.sq.knot4.scenario4)),
c(mean(diffsurv.sq.knot5.scenario1),mean(diffsurv.sq.knot5.scenario2),mean(diffsurv.sq.knot5.scenario3),mean(diffsurv.sq.knot5.scenario4)),
c(mean(diffsurv.sq.knot6.scenario1),mean(diffsurv.sq.knot6.scenario2),mean(diffsurv.sq.knot6.scenario3),mean(diffsurv.sq.knot6.scenario4)),
c(mean(diffsurv.sq.knot7.scenario1),mean(diffsurv.sq.knot7.scenario2),mean(diffsurv.sq.knot7.scenario3),mean(diffsurv.sq.knot7.scenario4)),
c(mean(diffsurv.sq.knot8.scenario1),mean(diffsurv.sq.knot8.scenario2),mean(diffsurv.sq.knot8.scenario3),mean(diffsurv.sq.knot8.scenario4)),
c(mean(diffsurv.sq.knot9.scenario1),mean(diffsurv.sq.knot9.scenario2),mean(diffsurv.sq.knot9.scenario3),mean(diffsurv.sq.knot9.scenario4)),
c(mean(diffsurv.sq.knot10.scenario1),mean(diffsurv.sq.knot10.scenario2),mean(diffsurv.sq.knot10.scenario3),mean(diffsurv.sq.knot10.scenario4)),
c(mean(diffsurv.sq.bestfit.scenario1),mean(diffsurv.sq.bestfit.scenario2),mean(diffsurv.sq.bestfit.scenario3),mean(diffsurv.sq.bestfit.scenario4)))

xtable(survdiff.sq*10000)

survdiff.abs=rbind(c(mean(diffsurv.abs.knot11.scenario1),mean(diffsurv.abs.knot11.scenario2),mean(diffsurv.abs.knot11.scenario3),mean(diffsurv.abs.knot11.scenario4)),
       c(mean(diffsurv.abs.knot1.scenario1),mean(diffsurv.abs.knot1.scenario2),mean(diffsurv.abs.knot1.scenario3),mean(diffsurv.abs.knot1.scenario4)),
      c(mean(diffsurv.abs.knot2.scenario1),mean(diffsurv.abs.knot2.scenario2),mean(diffsurv.abs.knot2.scenario3),mean(diffsurv.abs.knot2.scenario4)),
      c(mean(diffsurv.abs.knot3.scenario1),mean(diffsurv.abs.knot3.scenario2),mean(diffsurv.abs.knot3.scenario3),mean(diffsurv.abs.knot3.scenario4)),
      c(mean(diffsurv.abs.knot4.scenario1),mean(diffsurv.abs.knot4.scenario2),mean(diffsurv.abs.knot4.scenario3),mean(diffsurv.abs.knot4.scenario4)),
      c(mean(diffsurv.abs.knot5.scenario1),mean(diffsurv.abs.knot5.scenario2),mean(diffsurv.abs.knot5.scenario3),mean(diffsurv.abs.knot5.scenario4)),
      c(mean(diffsurv.abs.knot6.scenario1),mean(diffsurv.abs.knot6.scenario2),mean(diffsurv.abs.knot6.scenario3),mean(diffsurv.abs.knot6.scenario4)),
      c(mean(diffsurv.abs.knot7.scenario1),mean(diffsurv.abs.knot7.scenario2),mean(diffsurv.abs.knot7.scenario3),mean(diffsurv.abs.knot7.scenario4)),
      c(mean(diffsurv.abs.knot8.scenario1),mean(diffsurv.abs.knot8.scenario2),mean(diffsurv.abs.knot8.scenario3),mean(diffsurv.abs.knot8.scenario4)),
      c(mean(diffsurv.abs.knot9.scenario1),mean(diffsurv.abs.knot9.scenario2),mean(diffsurv.abs.knot9.scenario3),mean(diffsurv.abs.knot9.scenario4)),
      c(mean(diffsurv.abs.knot10.scenario1),mean(diffsurv.abs.knot10.scenario2),mean(diffsurv.abs.knot10.scenario3),mean(diffsurv.abs.knot10.scenario4)),
      c(mean(diffsurv.abs.bestfit.scenario1),mean(diffsurv.abs.bestfit.scenario2),mean(diffsurv.abs.bestfit.scenario3),mean(diffsurv.abs.bestfit.scenario4)))

xtable(survdiff.abs*10,digits=3)

#-------------------------
#table of differences between hazard curves
#-------------------------

hazdiff.sq=rbind(c(mean(diffhaz.sq.knot11.scenario1),mean(diffhaz.sq.knot11.scenario2),mean(diffhaz.sq.knot11.scenario3),mean(diffhaz.sq.knot11.scenario4)),
  c(mean(diffhaz.sq.knot1.scenario1),mean(diffhaz.sq.knot1.scenario2),mean(diffhaz.sq.knot1.scenario3),mean(diffhaz.sq.knot1.scenario4)),
      c(mean(diffhaz.sq.knot2.scenario1),mean(diffhaz.sq.knot2.scenario2),mean(diffhaz.sq.knot2.scenario3),mean(diffhaz.sq.knot2.scenario4)),
      c(mean(diffhaz.sq.knot3.scenario1),mean(diffhaz.sq.knot3.scenario2),mean(diffhaz.sq.knot3.scenario3),mean(diffhaz.sq.knot3.scenario4)),
      c(mean(diffhaz.sq.knot4.scenario1),mean(diffhaz.sq.knot4.scenario2),mean(diffhaz.sq.knot4.scenario3),mean(diffhaz.sq.knot4.scenario4)),
      c(mean(diffhaz.sq.knot5.scenario1),mean(diffhaz.sq.knot5.scenario2),mean(diffhaz.sq.knot5.scenario3),mean(diffhaz.sq.knot5.scenario4)),
      c(mean(diffhaz.sq.knot6.scenario1),mean(diffhaz.sq.knot6.scenario2),mean(diffhaz.sq.knot6.scenario3),mean(diffhaz.sq.knot6.scenario4)),
      c(mean(diffhaz.sq.knot7.scenario1),mean(diffhaz.sq.knot7.scenario2),mean(diffhaz.sq.knot7.scenario3),mean(diffhaz.sq.knot7.scenario4)),
      c(mean(diffhaz.sq.knot8.scenario1),mean(diffhaz.sq.knot8.scenario2),mean(diffhaz.sq.knot8.scenario3),mean(diffhaz.sq.knot8.scenario4)),
      c(mean(diffhaz.sq.knot9.scenario1),mean(diffhaz.sq.knot9.scenario2),mean(diffhaz.sq.knot9.scenario3),mean(diffhaz.sq.knot9.scenario4)),
      c(mean(diffhaz.sq.knot10.scenario1),mean(diffhaz.sq.knot10.scenario2),mean(diffhaz.sq.knot10.scenario3),mean(diffhaz.sq.knot10.scenario4)),
      c(mean(diffhaz.sq.bestfit.scenario1),mean(diffhaz.sq.bestfit.scenario2),mean(diffhaz.sq.bestfit.scenario3),mean(diffhaz.sq.bestfit.scenario4)))

xtable(hazdiff.sq*10000)

hazdiff.abs=rbind(c(mean(diffhaz.abs.knot11.scenario1),mean(diffhaz.abs.knot11.scenario2),mean(diffhaz.abs.knot11.scenario3),mean(diffhaz.abs.knot11.scenario4)),
  c(mean(diffhaz.abs.knot1.scenario1),mean(diffhaz.abs.knot1.scenario2),mean(diffhaz.abs.knot1.scenario3),mean(diffhaz.abs.knot1.scenario4)),
      c(mean(diffhaz.abs.knot2.scenario1),mean(diffhaz.abs.knot2.scenario2),mean(diffhaz.abs.knot2.scenario3),mean(diffhaz.abs.knot2.scenario4)),
      c(mean(diffhaz.abs.knot3.scenario1),mean(diffhaz.abs.knot3.scenario2),mean(diffhaz.abs.knot3.scenario3),mean(diffhaz.abs.knot3.scenario4)),
      c(mean(diffhaz.abs.knot4.scenario1),mean(diffhaz.abs.knot4.scenario2),mean(diffhaz.abs.knot4.scenario3),mean(diffhaz.abs.knot4.scenario4)),
      c(mean(diffhaz.abs.knot5.scenario1),mean(diffhaz.abs.knot5.scenario2),mean(diffhaz.abs.knot5.scenario3),mean(diffhaz.abs.knot5.scenario4)),
      c(mean(diffhaz.abs.knot6.scenario1),mean(diffhaz.abs.knot6.scenario2),mean(diffhaz.abs.knot6.scenario3),mean(diffhaz.abs.knot6.scenario4)),
      c(mean(diffhaz.abs.knot7.scenario1),mean(diffhaz.abs.knot7.scenario2),mean(diffhaz.abs.knot7.scenario3),mean(diffhaz.abs.knot7.scenario4)),
      c(mean(diffhaz.abs.knot8.scenario1),mean(diffhaz.abs.knot8.scenario2),mean(diffhaz.abs.knot8.scenario3),mean(diffhaz.abs.knot8.scenario4)),
      c(mean(diffhaz.abs.knot9.scenario1),mean(diffhaz.abs.knot9.scenario2),mean(diffhaz.abs.knot9.scenario3),mean(diffhaz.abs.knot9.scenario4)),
      c(mean(diffhaz.abs.knot10.scenario1),mean(diffhaz.abs.knot10.scenario2),mean(diffhaz.abs.knot10.scenario3),mean(diffhaz.abs.knot10.scenario4)),
      c(mean(diffhaz.abs.bestfit.scenario1),mean(diffhaz.abs.bestfit.scenario2),mean(diffhaz.abs.bestfit.scenario3),mean(diffhaz.abs.bestfit.scenario4)))

xtable(hazdiff.abs*10,digits=2)

#-------------------------
#table of best fitting models
#-------------------------

xtable(cbind(table(factor(bestfit.scenario1,levels=c(1:11))),
             table(factor(bestfit.scenario2,levels=c(1:11))),
             table(factor(bestfit.scenario3,levels=c(1:11))),
             table(factor(bestfit.scenario4,levels=c(1:11))))/10)

#-----------------------------------
#-----------------------------------
#plots of estimated survival and hazard curves (under the model with lowest AIC)
#-----------------------------------
#-----------------------------------

options(expressions = 10000)

#--------
#scenario 1
#--------

#survival 

sim1.0=which(bestfit.scenario1==11)
coeflines.0 <-
  alply(as.matrix(sim1.0), 1, function(sim) {
    stat_function(fun=function(x){surv.func.0knot.scenario1(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim1.",k,"=which(bestfit.scenario1==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim1.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){surv.func.",k,"knot.scenario1(x,sim)}, colour=\"grey\")
                         })")))
}

survplot.1=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.surv.func1(x)) + 
  labs(x="Follow-up time",y="Survival function")+ggtitle("Scenario 1")

#hazard

sim1.0=which(bestfit.scenario1==11)
coeflines.0 <-
  alply(as.matrix(sim1.0), 1, function(sim) {
    stat_function(fun=function(x){haz.func.0knot.scenario1(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim1.",k,"=which(bestfit.scenario1==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim1.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){haz.func.",k,"knot.scenario1(x,sim)}, colour=\"grey\")
                         })")))
}

hazplot.1=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.haz.func1(x)) + 
  labs(x="Follow-up time",y="Hazard function")+ggtitle("Scenario 1")

#--------
#scenario 2
#--------

#survival 

sim2.0=which(bestfit.scenario2==11)
coeflines.0 <-
  alply(as.matrix(sim2.0), 1, function(sim) {
    stat_function(fun=function(x){surv.func.0knot.scenario2(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim2.",k,"=which(bestfit.scenario2==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim2.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){surv.func.",k,"knot.scenario2(x,sim)}, colour=\"grey\")
                         })")))
}

survplot.2=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.surv.func2(x)) + 
  labs(x="Follow-up time",y="Survival function")+ggtitle("Scenario 2")

#hazard

sim2.0=which(bestfit.scenario2==11)
coeflines.0 <-
  alply(as.matrix(sim2.0), 1, function(sim) {
    stat_function(fun=function(x){haz.func.0knot.scenario2(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim2.",k,"=which(bestfit.scenario2==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim2.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){haz.func.",k,"knot.scenario2(x,sim)}, colour=\"grey\")
                         })")))
}

hazplot.2=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.haz.func2(x)) + 
  labs(x="Follow-up time",y="Hazard function")+ggtitle("Scenario 2")

#--------
#scenario 3
#--------

#survival 

sim3.0=which(bestfit.scenario3==11)
coeflines.0 <-
  alply(as.matrix(sim3.0), 1, function(sim) {
    stat_function(fun=function(x){surv.func.0knot.scenario3(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim3.",k,"=which(bestfit.scenario3==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim2.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){surv.func.",k,"knot.scenario3(x,sim)}, colour=\"grey\")
                         })")))
}

survplot.3=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.surv.func3(x)) + 
  labs(x="Follow-up time",y="Survival function")+ggtitle("Scenario 3")

#hazard

sim3.0=which(bestfit.scenario3==11)
coeflines.0 <-
  alply(as.matrix(sim3.0), 1, function(sim) {
    stat_function(fun=function(x){haz.func.0knot.scenario3(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim3.",k,"=which(bestfit.scenario3==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim3.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){haz.func.",k,"knot.scenario3(x,sim)}, colour=\"grey\")
                         })")))
}

hazplot.3=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.haz.func3(x)) + 
  labs(x="Follow-up time",y="Hazard function")+ggtitle("Scenario 3")

#--------
#scenario 4
#--------

#survival 

sim4.0=which(bestfit.scenario4==11)
coeflines.0 <-
  alply(as.matrix(sim4.0), 1, function(sim) {
    stat_function(fun=function(x){surv.func.0knot.scenario4(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim4.",k,"=which(bestfit.scenario4==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim2.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){surv.func.",k,"knot.scenario4(x,sim)}, colour=\"grey\")
                         })")))
}

#coeflines.0knots =alply(as.matrix(1:1000), 1, function(sim) {
#                         stat_function(fun=function(x){surv.func.0knot.scenario4(x,sim)}, colour="blue")})

#coeflines.1knots =alply(as.matrix(1:1000), 1, function(sim) {
#  stat_function(fun=function(x){surv.func.1knot.scenario4(x,sim)}, colour="blue")})

#coeflines.10knots =alply(as.matrix(1:1000), 1, function(sim) {
#  stat_function(fun=function(x){surv.func.10knot.scenario4(x,sim)}, colour="green")})

survplot.4=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10  + stat_function(fun=function(x) true.surv.func4(x)) + 
  labs(x="Follow-up time",y="Survival function")+ggtitle("Scenario 4")

#hazard

sim4.0=which(bestfit.scenario4==11)
coeflines.0 <-
  alply(as.matrix(sim4.0), 1, function(sim) {
    stat_function(fun=function(x){haz.func.0knot.scenario4(x,sim)}, colour="grey")
  })

for(k in 1:10){
  eval(parse(text=paste0("sim4.",k,"=which(bestfit.scenario4==",k,")")))
  eval(parse(text=paste0("coeflines.",k," =
                         alply(as.matrix(sim4.",k,"), 1, function(sim) {
                         stat_function(fun=function(x){haz.func.",k,"knot.scenario4(x,sim)}, colour=\"grey\")
                         })")))
}

hazplot.4=ggplot(data.frame(x = c(0, 10)), aes(x)) + coeflines.0 + coeflines.1 + coeflines.2 + coeflines.3 + coeflines.4 + coeflines.5 +
  coeflines.6 + coeflines.7 + coeflines.8 + coeflines.9 + coeflines.10 + stat_function(fun=function(x) true.haz.func4(x)) + 
  labs(x="Follow-up time",y="Hazard function")+ggtitle("Scenario 4")

#--------
#combine plots
#--------

windows(10,18)
grid.arrange(survplot.1,hazplot.1,survplot.2,hazplot.2,survplot.3,hazplot.3,survplot.4,hazplot.4,nrow = 4)

ggsave("H:/CF_UK_lifeexp/BARS_sim/plots/sim_results_plots.jpeg")

#-----------------------------------
#-----------------------------------
#coverage
#-----------------------------------
#-----------------------------------

#0 knots - use 11 in place of zero for later convenience below

surv.varcov.11knot=function(t,flexmod.params.1,flexmod.params.2,cov.mat){
  vmat=as.matrix(c(1,log(t)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.11knot.scenario",scen,"=function(t,sim){log(-log(surv.func.0knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.11knot(t,flexmod.params",scen,"[sim,11,1],flexmod.params",scen,"[sim,11,2],
                        flexmod.varcov",scen,"[sim,11,1:2,1:2]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.11knot.scenario",scen,"=function(t,sim){log(-log(surv.func.0knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.11knot(t,flexmod.params",scen,"[sim,11,1],flexmod.params",scen,"[sim,11,2],
                        flexmod.varcov",scen,"[sim,11,1:2,1:2]))
}")))
}

#1 knots

surv.varcov.1knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.1knot.scenario",scen,"=function(t,sim){log(-log(surv.func.1knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.1knot(t,flexmod.params",scen,"[sim,1,1],flexmod.params",scen,"[sim,1,2],flexmod.params",scen,"[sim,1,3],
                         flexmod.knots",scen,"[sim,1,1],flexmod.knots",scen,"[sim,1,3],flexmod.knots",scen,"[sim,1,2],
                        flexmod.varcov",scen,"[sim,1,1:3,1:3]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.1knot.scenario",scen,"=function(t,sim){log(-log(surv.func.1knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.1knot(t,flexmod.params",scen,"[sim,1,1],flexmod.params",scen,"[sim,1,2],flexmod.params",scen,"[sim,1,3],
                         flexmod.knots",scen,"[sim,1,1],flexmod.knots",scen,"[sim,1,3],flexmod.knots",scen,"[sim,1,2],
                        flexmod.varcov",scen,"[sim,1,1:3,1:3]))
}")))
}

#2 knots

surv.varcov.2knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.2knot.scenario",scen,"=function(t,sim){log(-log(surv.func.2knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.2knot(t,flexmod.params",scen,"[sim,2,1],flexmod.params",scen,"[sim,2,2],flexmod.params",scen,"[sim,2,3],flexmod.params",scen,"[sim,2,4],
                         flexmod.knots",scen,"[sim,2,1],flexmod.knots",scen,"[sim,2,4],flexmod.knots",scen,"[sim,2,2],flexmod.knots",scen,"[sim,2,3],
                        flexmod.varcov",scen,"[sim,2,1:4,1:4]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.2knot.scenario",scen,"=function(t,sim){log(-log(surv.func.2knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.2knot(t,flexmod.params",scen,"[sim,2,1],flexmod.params",scen,"[sim,2,2],flexmod.params",scen,"[sim,2,3],flexmod.params",scen,"[sim,2,4],
                         flexmod.knots",scen,"[sim,2,1],flexmod.knots",scen,"[sim,2,4],flexmod.knots",scen,"[sim,2,2],flexmod.knots",scen,"[sim,2,3],
                        flexmod.varcov",scen,"[sim,2,1:4,1:4]))
}")))
}

#3 knots

surv.varcov.3knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.3knot.scenario",scen,"=function(t,sim){log(-log(surv.func.3knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.3knot(t,flexmod.params",scen,"[sim,3,1],flexmod.params",scen,"[sim,3,2],flexmod.params",scen,"[sim,3,3],flexmod.params",scen,"[sim,3,4],
                         flexmod.params",scen,"[sim,3,5],
                         flexmod.knots",scen,"[sim,3,1],flexmod.knots",scen,"[sim,3,5],flexmod.knots",scen,"[sim,3,2],flexmod.knots",scen,"[sim,3,3],
                         flexmod.knots",scen,"[sim,3,4],
                        flexmod.varcov",scen,"[sim,3,1:5,1:5]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.3knot.scenario",scen,"=function(t,sim){log(-log(surv.func.3knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.3knot(t,flexmod.params",scen,"[sim,3,1],flexmod.params",scen,"[sim,3,2],flexmod.params",scen,"[sim,3,3],flexmod.params",scen,"[sim,3,4],
                         flexmod.params",scen,"[sim,3,5],
                         flexmod.knots",scen,"[sim,3,1],flexmod.knots",scen,"[sim,3,5],flexmod.knots",scen,"[sim,3,2],flexmod.knots",scen,"[sim,3,3],
                         flexmod.knots",scen,"[sim,3,4],
                        flexmod.varcov",scen,"[sim,3,1:5,1:5]))
}")))
}


#4 knots

surv.varcov.4knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.4knot.scenario",scen,"=function(t,sim){log(-log(surv.func.4knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.4knot(t,flexmod.params",scen,"[sim,4,1],flexmod.params",scen,"[sim,4,2],flexmod.params",scen,"[sim,4,3],flexmod.params",scen,"[sim,4,4],
                         flexmod.params",scen,"[sim,4,5],flexmod.params",scen,"[sim,4,6],
                         flexmod.knots",scen,"[sim,4,1],flexmod.knots",scen,"[sim,4,6],flexmod.knots",scen,"[sim,4,2],flexmod.knots",scen,"[sim,4,3],
                         flexmod.knots",scen,"[sim,4,4],flexmod.knots",scen,"[sim,4,5],
                        flexmod.varcov",scen,"[sim,4,1:6,1:6]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.4knot.scenario",scen,"=function(t,sim){log(-log(surv.func.4knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.4knot(t,flexmod.params",scen,"[sim,4,1],flexmod.params",scen,"[sim,4,2],flexmod.params",scen,"[sim,4,3],flexmod.params",scen,"[sim,4,4],
                         flexmod.params",scen,"[sim,4,5],flexmod.params",scen,"[sim,4,6],
                         flexmod.knots",scen,"[sim,4,1],flexmod.knots",scen,"[sim,4,6],flexmod.knots",scen,"[sim,4,2],flexmod.knots",scen,"[sim,4,3],
                         flexmod.knots",scen,"[sim,4,4],flexmod.knots",scen,"[sim,4,5],
                        flexmod.varcov",scen,"[sim,4,1:6,1:6]))
}")))
}


#5 knots

surv.varcov.5knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                           flexmod.params.7,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.5knot.scenario",scen,"=function(t,sim){log(-log(surv.func.5knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.5knot(t,flexmod.params",scen,"[sim,5,1],flexmod.params",scen,"[sim,5,2],flexmod.params",scen,"[sim,5,3],flexmod.params",scen,"[sim,5,4],
                         flexmod.params",scen,"[sim,5,5],flexmod.params",scen,"[sim,5,6],flexmod.params",scen,"[sim,5,7],
                         flexmod.knots",scen,"[sim,5,1],flexmod.knots",scen,"[sim,5,7],flexmod.knots",scen,"[sim,5,2],flexmod.knots",scen,"[sim,5,3],
                         flexmod.knots",scen,"[sim,5,4],flexmod.knots",scen,"[sim,5,5],flexmod.knots",scen,"[sim,5,6],
                        flexmod.varcov",scen,"[sim,5,1:7,1:7]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.5knot.scenario",scen,"=function(t,sim){log(-log(surv.func.5knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.5knot(t,flexmod.params",scen,"[sim,5,1],flexmod.params",scen,"[sim,5,2],flexmod.params",scen,"[sim,5,3],flexmod.params",scen,"[sim,5,4],
                         flexmod.params",scen,"[sim,5,5],flexmod.params",scen,"[sim,5,6],flexmod.params",scen,"[sim,5,7],
                         flexmod.knots",scen,"[sim,5,1],flexmod.knots",scen,"[sim,5,7],flexmod.knots",scen,"[sim,5,2],flexmod.knots",scen,"[sim,5,3],
                         flexmod.knots",scen,"[sim,5,4],flexmod.knots",scen,"[sim,5,5],flexmod.knots",scen,"[sim,5,6],
                        flexmod.varcov",scen,"[sim,5,1:7,1:7]))
}")))
}


#6 knots

surv.varcov.6knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                           flexmod.params.7,flexmod.params.8,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                           flexmod.knots.7,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.6knot.scenario",scen,"=function(t,sim){log(-log(surv.func.6knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.6knot(t,flexmod.params",scen,"[sim,6,1],flexmod.params",scen,"[sim,6,2],flexmod.params",scen,"[sim,6,3],flexmod.params",scen,"[sim,6,4],
                         flexmod.params",scen,"[sim,6,5],flexmod.params",scen,"[sim,6,6],flexmod.params",scen,"[sim,6,7],flexmod.params",scen,"[sim,6,8],
                         flexmod.knots",scen,"[sim,6,1],flexmod.knots",scen,"[sim,6,8],flexmod.knots",scen,"[sim,6,2],flexmod.knots",scen,"[sim,6,3],
                         flexmod.knots",scen,"[sim,6,4],flexmod.knots",scen,"[sim,6,5],flexmod.knots",scen,"[sim,6,6],flexmod.knots",scen,"[sim,6,7],
                        flexmod.varcov",scen,"[sim,6,1:8,1:8]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.6knot.scenario",scen,"=function(t,sim){log(-log(surv.func.6knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.6knot(t,flexmod.params",scen,"[sim,6,1],flexmod.params",scen,"[sim,6,2],flexmod.params",scen,"[sim,6,3],flexmod.params",scen,"[sim,6,4],
                         flexmod.params",scen,"[sim,6,5],flexmod.params",scen,"[sim,6,6],flexmod.params",scen,"[sim,6,7],flexmod.params",scen,"[sim,6,8],
                         flexmod.knots",scen,"[sim,6,1],flexmod.knots",scen,"[sim,6,8],flexmod.knots",scen,"[sim,6,2],flexmod.knots",scen,"[sim,6,3],
                         flexmod.knots",scen,"[sim,6,4],flexmod.knots",scen,"[sim,6,5],flexmod.knots",scen,"[sim,6,6],flexmod.knots",scen,"[sim,6,7],
                        flexmod.varcov",scen,"[sim,6,1:8,1:8]))
}")))
}


#7 knots

surv.varcov.7knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                           flexmod.params.7,flexmod.params.8,flexmod.params.9,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                           flexmod.knots.7,flexmod.knots.8,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.7knot.scenario",scen,"=function(t,sim){log(-log(surv.func.7knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.7knot(t,flexmod.params",scen,"[sim,7,1],flexmod.params",scen,"[sim,7,2],flexmod.params",scen,"[sim,7,3],flexmod.params",scen,"[sim,7,4],
                         flexmod.params",scen,"[sim,7,5],flexmod.params",scen,"[sim,7,6],flexmod.params",scen,"[sim,7,7],flexmod.params",scen,"[sim,7,8],
                         flexmod.params",scen,"[sim,7,9],
                         flexmod.knots",scen,"[sim,7,1],flexmod.knots",scen,"[sim,7,9],flexmod.knots",scen,"[sim,7,2],flexmod.knots",scen,"[sim,7,3],
                         flexmod.knots",scen,"[sim,7,4],flexmod.knots",scen,"[sim,7,5],flexmod.knots",scen,"[sim,7,6],flexmod.knots",scen,"[sim,7,7],
                         flexmod.knots",scen,"[sim,7,8],
                        flexmod.varcov",scen,"[sim,7,1:9,1:9]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.7knot.scenario",scen,"=function(t,sim){log(-log(surv.func.7knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.7knot(t,flexmod.params",scen,"[sim,7,1],flexmod.params",scen,"[sim,7,2],flexmod.params",scen,"[sim,7,3],flexmod.params",scen,"[sim,7,4],
                         flexmod.params",scen,"[sim,7,5],flexmod.params",scen,"[sim,7,6],flexmod.params",scen,"[sim,7,7],flexmod.params",scen,"[sim,7,8],
                         flexmod.params",scen,"[sim,7,9],
                         flexmod.knots",scen,"[sim,7,1],flexmod.knots",scen,"[sim,7,9],flexmod.knots",scen,"[sim,7,2],flexmod.knots",scen,"[sim,7,3],
                         flexmod.knots",scen,"[sim,7,4],flexmod.knots",scen,"[sim,7,5],flexmod.knots",scen,"[sim,7,6],flexmod.knots",scen,"[sim,7,7],
                         flexmod.knots",scen,"[sim,7,8],
                        flexmod.varcov",scen,"[sim,7,1:9,1:9]))
}")))
}


#8 knots

surv.varcov.8knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                           flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,
                           flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                           flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.8knot.scenario",scen,"=function(t,sim){log(-log(surv.func.8knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.8knot(t,flexmod.params",scen,"[sim,8,1],flexmod.params",scen,"[sim,8,2],flexmod.params",scen,"[sim,8,3],flexmod.params",scen,"[sim,8,4],
                         flexmod.params",scen,"[sim,8,5],flexmod.params",scen,"[sim,8,6],flexmod.params",scen,"[sim,8,7],flexmod.params",scen,"[sim,8,8],
                         flexmod.params",scen,"[sim,8,9],flexmod.params",scen,"[sim,8,10],
                         flexmod.knots",scen,"[sim,8,1],flexmod.knots",scen,"[sim,8,10],flexmod.knots",scen,"[sim,8,2],flexmod.knots",scen,"[sim,8,3],
                         flexmod.knots",scen,"[sim,8,4],flexmod.knots",scen,"[sim,8,5],flexmod.knots",scen,"[sim,8,6],flexmod.knots",scen,"[sim,8,7],
                         flexmod.knots",scen,"[sim,8,8],flexmod.knots",scen,"[sim,8,9],
                        flexmod.varcov",scen,"[sim,8,1:10,1:10]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.8knot.scenario",scen,"=function(t,sim){log(-log(surv.func.8knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.8knot(t,flexmod.params",scen,"[sim,8,1],flexmod.params",scen,"[sim,8,2],flexmod.params",scen,"[sim,8,3],flexmod.params",scen,"[sim,8,4],
                         flexmod.params",scen,"[sim,8,5],flexmod.params",scen,"[sim,8,6],flexmod.params",scen,"[sim,8,7],flexmod.params",scen,"[sim,8,8],
                         flexmod.params",scen,"[sim,8,9],flexmod.params",scen,"[sim,8,10],
                         flexmod.knots",scen,"[sim,8,1],flexmod.knots",scen,"[sim,8,10],flexmod.knots",scen,"[sim,8,2],flexmod.knots",scen,"[sim,8,3],
                         flexmod.knots",scen,"[sim,8,4],flexmod.knots",scen,"[sim,8,5],flexmod.knots",scen,"[sim,8,6],flexmod.knots",scen,"[sim,8,7],
                         flexmod.knots",scen,"[sim,8,8],flexmod.knots",scen,"[sim,8,9],
                        flexmod.varcov",scen,"[sim,8,1:10,1:10]))
}")))
}

#9 knots

surv.varcov.9knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                            flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,flexmod.params.11,
                            flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                            flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,flexmod.knots.10,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.9knot.scenario",scen,"=function(t,sim){log(-log(surv.func.9knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.9knot(t,flexmod.params",scen,"[sim,9,1],flexmod.params",scen,"[sim,9,2],flexmod.params",scen,"[sim,9,3],flexmod.params",scen,"[sim,9,4],
                         flexmod.params",scen,"[sim,9,5],flexmod.params",scen,"[sim,9,6],flexmod.params",scen,"[sim,9,7],flexmod.params",scen,"[sim,9,8],
                         flexmod.params",scen,"[sim,9,9],flexmod.params",scen,"[sim,9,10],flexmod.params",scen,"[sim,9,11],
                         flexmod.knots",scen,"[sim,9,1],flexmod.knots",scen,"[sim,9,11],flexmod.knots",scen,"[sim,9,2],flexmod.knots",scen,"[sim,9,3],
                         flexmod.knots",scen,"[sim,9,4],flexmod.knots",scen,"[sim,9,5],flexmod.knots",scen,"[sim,9,6],flexmod.knots",scen,"[sim,9,7],
                         flexmod.knots",scen,"[sim,9,8],flexmod.knots",scen,"[sim,9,9],flexmod.knots",scen,"[sim,9,10],
                        flexmod.varcov",scen,"[sim,9,1:11,1:11]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.9knot.scenario",scen,"=function(t,sim){log(-log(surv.func.9knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.9knot(t,flexmod.params",scen,"[sim,9,1],flexmod.params",scen,"[sim,9,2],flexmod.params",scen,"[sim,9,3],flexmod.params",scen,"[sim,9,4],
                         flexmod.params",scen,"[sim,9,5],flexmod.params",scen,"[sim,9,6],flexmod.params",scen,"[sim,9,7],flexmod.params",scen,"[sim,9,8],
                         flexmod.params",scen,"[sim,9,9],flexmod.params",scen,"[sim,9,10],flexmod.params",scen,"[sim,9,11],
                         flexmod.knots",scen,"[sim,9,1],flexmod.knots",scen,"[sim,9,11],flexmod.knots",scen,"[sim,9,2],flexmod.knots",scen,"[sim,9,3],
                         flexmod.knots",scen,"[sim,9,4],flexmod.knots",scen,"[sim,9,5],flexmod.knots",scen,"[sim,9,6],flexmod.knots",scen,"[sim,9,7],
                         flexmod.knots",scen,"[sim,9,8],flexmod.knots",scen,"[sim,9,9],flexmod.knots",scen,"[sim,9,10],
                        flexmod.varcov",scen,"[sim,9,1:11,1:11]))
}")))
}


#10 knots

surv.varcov.10knot=function(t,flexmod.params.1,flexmod.params.2,flexmod.params.3,flexmod.params.4,flexmod.params.5,flexmod.params.6,
                          flexmod.params.7,flexmod.params.8,flexmod.params.9,flexmod.params.10,flexmod.params.11,flexmod.params.12,
                          flexmod.knots.1,flexmod.knots.max,flexmod.knots.2,flexmod.knots.3,flexmod.knots.4,flexmod.knots.5,flexmod.knots.6,
                          flexmod.knots.7,flexmod.knots.8,flexmod.knots.9,flexmod.knots.10,flexmod.knots.11,cov.mat){
  vmat=as.matrix(c(1,log(t),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.2),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.3),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.4),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.5),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.6),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.7),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.8),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.9),
                   v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.10),v(t,flexmod.knots.1,flexmod.knots.max,flexmod.knots.11)),ncol=1)
  t(vmat)%*%cov.mat%*%vmat
}


for(scen in 1:4){
  eval(parse(text=paste0("surv.upper.10knot.scenario",scen,"=function(t,sim){log(-log(surv.func.10knot.scenario",scen,"(t,sim)))+1.96*sqrt(
                         surv.varcov.10knot(t,flexmod.params",scen,"[sim,10,1],flexmod.params",scen,"[sim,10,2],flexmod.params",scen,"[sim,10,3],flexmod.params",scen,"[sim,10,4],
                         flexmod.params",scen,"[sim,10,5],flexmod.params",scen,"[sim,10,6],flexmod.params",scen,"[sim,10,7],flexmod.params",scen,"[sim,10,8],
                         flexmod.params",scen,"[sim,10,9],flexmod.params",scen,"[sim,10,10],flexmod.params",scen,"[sim,10,11],flexmod.params",scen,"[sim,10,12],
                         flexmod.knots",scen,"[sim,10,1],flexmod.knots",scen,"[sim,10,12],flexmod.knots",scen,"[sim,10,2],flexmod.knots",scen,"[sim,10,3],
                         flexmod.knots",scen,"[sim,10,4],flexmod.knots",scen,"[sim,10,5],flexmod.knots",scen,"[sim,10,6],flexmod.knots",scen,"[sim,10,7],
                         flexmod.knots",scen,"[sim,10,8],flexmod.knots",scen,"[sim,10,9],flexmod.knots",scen,"[sim,10,10],flexmod.knots",scen,"[sim,10,11],
                        flexmod.varcov",scen,"[sim,10,1:12,1:12]))
}")))
}

for(scen in 1:4){
  eval(parse(text=paste0("surv.lower.10knot.scenario",scen,"=function(t,sim){log(-log(surv.func.10knot.scenario",scen,"(t,sim)))-1.96*sqrt(
                         surv.varcov.10knot(t,flexmod.params",scen,"[sim,10,1],flexmod.params",scen,"[sim,10,2],flexmod.params",scen,"[sim,10,3],flexmod.params",scen,"[sim,10,4],
                         flexmod.params",scen,"[sim,10,5],flexmod.params",scen,"[sim,10,6],flexmod.params",scen,"[sim,10,7],flexmod.params",scen,"[sim,10,8],
                         flexmod.params",scen,"[sim,10,9],flexmod.params",scen,"[sim,10,10],flexmod.params",scen,"[sim,10,11],flexmod.params",scen,"[sim,10,12],
                         flexmod.knots",scen,"[sim,10,1],flexmod.knots",scen,"[sim,10,12],flexmod.knots",scen,"[sim,10,2],flexmod.knots",scen,"[sim,10,3],
                         flexmod.knots",scen,"[sim,10,4],flexmod.knots",scen,"[sim,10,5],flexmod.knots",scen,"[sim,10,6],flexmod.knots",scen,"[sim,10,7],
                         flexmod.knots",scen,"[sim,10,8],flexmod.knots",scen,"[sim,10,9],flexmod.knots",scen,"[sim,10,10],flexmod.knots",scen,"[sim,10,11],
                        flexmod.varcov",scen,"[sim,10,1:12,1:12]))
}")))
}

t=1
sim=1

t.grid=seq(0.1,10,0.1)

sapply(t.grid,function(t){log(-log(true.surv.func1(t)))>surv.lower.11knot.scenario1(t,sim) & log(-log(true.surv.func1(t)))<surv.upper.11knot.scenario1(t,sim)})



for(scen in 1:4){
  eval(parse(text=paste0("coverage.bestfit.scenario",scen,"=rep(0,1000)")))  
  for(sim in 1:1000){
    eval(parse(text=paste0("bestfit=bestfit.scenario",scen,"[sim]")))
    eval(parse(text=paste0("coverage.bestfit.scenario",scen,"[sim]=mean(sapply(t.grid,function(t){log(-log(true.surv.func",1,"(t)))>surv.lower.",bestfit,"knot.scenario",1,"(t,sim) & 
                       log(-log(true.surv.func",1,"(t)))<surv.upper.",bestfit,"knot.scenario",1,"(t,sim)}))")))
  }
}

sum(1-is.na(coverage.bestfit.scenario1))/1000
sum(1-is.na(coverage.bestfit.scenario2))/1000
sum(1-is.na(coverage.bestfit.scenario3))/1000
sum(1-is.na(coverage.bestfit.scenario4))/1000


mean(coverage.bestfit.scenario1,na.rm=T)
mean(coverage.bestfit.scenario2,na.rm=T)
mean(coverage.bestfit.scenario3,na.rm=T)
mean(coverage.bestfit.scenario4,na.rm=T)


