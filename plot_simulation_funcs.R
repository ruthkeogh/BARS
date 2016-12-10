
#------------------------------------------
#libraries
#------------------------------------------

library(ggplot2)
library(gridExtra)

#------------------------------------------
#define functions needed for plots
#------------------------------------------

#survivor function: scenario 1

surv.func1=function(lambda1,gamma1,t){
  exp(-lambda1*t^(gamma1))
}

#hazard function: scenario 1

hazard.func1=function(lambda1,gamma1,t){
  (lambda1*gamma1*(t^(gamma1-1))*exp(-lambda1*t^(gamma1)))/
    (exp(-lambda1*t^(gamma1)))
}


#survivor function: scenarios 2-4

surv.func=function(p,lambda1,lambda2,gamma1,gamma2,t){
  p*exp(-lambda1*t^(gamma1))+(1-p)*exp(-lambda2*t^(gamma2))
}

#hazard function: scenarios 2-4

hazard.func=function(p,lambda1,lambda2,gamma1,gamma2,t){
  (p*lambda1*gamma1*(t^(gamma1-1))*exp(-lambda1*t^(gamma1))+(1-p)*lambda2*gamma2*(t^(gamma2-1))*exp(-lambda2*t^(gamma2)))/
  (p*exp(-lambda1*t^(gamma1))+(1-p)*exp(-lambda2*t^(gamma2)))
}

#------------------------------------------
#plot survivor function
#------------------------------------------

windows(10,7)
cols <- c("Scenario 1"="black", "Scenario 2"="black", "Scenario 3"="grey", "Scenario 4"="grey")
ltys=c("Scenario 1"=1,"Scenario 2"=2,"Scenario 3"=1,"Scenario 4"=2)
myplot.surv=ggplot(data.frame(x=seq(0,10,1),y=seq(0,10,1)),aes(x))+
  stat_function(fun=function(x) surv.func1(0.6,0.8,x), aes(colour = "Scenario 1",linetype="Scenario 1"))+
  stat_function(fun=function(x) surv.func(0.2,0.2,1.6,0.8,1,x), aes(colour = "Scenario 2",linetype="Scenario 2"))+
  stat_function(fun=function(x) surv.func(0.5,1,1,1.5,0.5,x), aes(colour = "Scenario 3",linetype="Scenario 3"))+
  stat_function(fun=function(x) surv.func(0.7,0.03,0.3,1.9,2.5,x), aes(colour = "Scenario 4",linetype="Scenario 4"))+
  labs(x="Follow-up time",y="Survival probability")+theme(axis.title=element_text(size = 15))+
  theme(axis.text.x=element_text(size=12,margin=margin(20,20,10,0)),
      axis.text.y=element_text(size=12,margin=margin(20,20,30,10)))+
  scale_colour_manual(name=element_blank(),values=cols,guide = guide_legend(override.aes = list(
    linetype = c(1,2,1,2),
    size = c(0, 0, 0, 0))))+
  scale_linetype_manual(values = ltys,guide = FALSE)+theme(legend.text=element_text(size=rel(1.2)))

#ggsave("H:/CF_UK_lifeexp/BARS_sim/plots/surv_func_plot.jpeg")

#------------------------------------------
#plot hazard function
#------------------------------------------

windows(10,7)
cols <- c("Scenario 1"="black", "Scenario 2"="black", "Scenario 3"="grey", "Scenario 4"="grey")
ltys=c("Scenario 1"=1,"Scenario 2"=2,"Scenario 3"=1,"Scenario 4"=2)
myplot.haz=ggplot(data.frame(x=seq(0,10,1),y=seq(0,10,1)),aes(x))+
  stat_function(fun=function(x) hazard.func1(0.6,0.8,x), aes(colour = "Scenario 1",linetype="Scenario 1"))+
  stat_function(fun=function(x) hazard.func(0.2,0.2,1.6,0.8,1,x), aes(colour = "Scenario 2",linetype="Scenario 2"))+
  stat_function(fun=function(x) hazard.func(0.5,1,1,1.5,0.5,x), aes(colour = "Scenario 3",linetype="Scenario 3"))+
  stat_function(fun=function(x) hazard.func(0.7,0.03,0.3,1.9,2.5,x), aes(colour = "Scenario 4",linetype="Scenario 4"))+
  labs(x="Follow-up time",y="Hazard function")+theme(axis.title=element_text(size = 15))+
  theme(axis.text.x=element_text(size=12,margin=margin(20,20,10,0)),
        axis.text.y=element_text(size=12,margin=margin(20,20,30,10)))+
  scale_colour_manual(name=element_blank(),values=cols,guide = guide_legend(override.aes = list(
    linetype = c(1,2,1,2),
    size = c(0, 0, 0, 0))))+
  scale_linetype_manual(values = ltys,guide = FALSE)+theme(legend.text=element_text(size=rel(1.2)))


windows(12,6)
grid.arrange(myplot.surv, myplot.haz,nrow = 1)

ggsave("H:/CF_UK_lifeexp/BARS_sim/plots/hazard_func_plot.jpeg")
