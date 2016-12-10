
#-----------------
#sample size (i.e. number of individuals in each simulated data set)
#-----------------

n=5000

#-----------------
#values for parameters used to simulate the data
#-----------------

#note that 'scenario' is set in the file 'standard_analysis.R'
#scenario is 1,2,3 or 4

#scenario 1
if(scenario==1){
lambda1=0.6
gamma1=0.8}

#scenario 2
if(scenario==2){
p=0.2
lambda1=0.2
lambda2=1.6
gamma1=0.8
gamma2=1}

#scenario 3
if(scenario==3){
p=0.5
lambda1=1
lambda2=1
gamma1=1.5
gamma2=0.5}

#scenario 4
if(scenario==4){
p=0.7
lambda1=0.03
lambda2=0.3
gamma1=1.9
gamma2=2.5}

#-----------------
#simulate from two weibull distributions used in the mixture of weibulls
#-----------------

#weibull distribution 1

unif1=runif(n)
t.weib1=((-1/lambda1)*log(unif1))^(1/gamma1)

#weibull distribution 2

if(scenario>1){
  unif2=runif(n)
  t.weib2=((-1/lambda2)*log(unif2))^(1/gamma2)
  }

#-----------------
#generate the mixture distribution
#-----------------

if(scenario==1){
  t=t.weib1
  }

if(scenario>1){
  bern=rbinom(n,1,p)
  t=bern*t.weib1+(1-bern)*t.weib2
  }

#-----------------
#censor event times at time 10
#and generate event indicator
#-----------------

#event indicator
d=ifelse(t>10,0,1)

#censor event times
t=ifelse(t>10,10,t)








  
  
  

