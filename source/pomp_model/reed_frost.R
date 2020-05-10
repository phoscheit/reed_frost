library(tidyverse)
library(pomp)

# First, specify the simulation function

reed_frost.proc.sim <- function(gamma,h0,beta,A,pI,S,E,DI,H,...){
  infection_prob <- 1-A*exp(-beta*H)
  new_infected <- rbinom(size = E,prob = pI,n=1)
  new_exposed <- rbinom(size = S,prob = infection_prob,n=1)
  c(E=E+new_exposed-new_infected,
    DI=new_infected,
    H=gamma*H+h0*DI,
    S=S-new_exposed)
}

reed_frost.init.sim <- function(s0,...){
  c(E=0,DI=0,H=0,S=s0)
}

reed_frost.meas.sim <- function(DI,sym_prob,...){
  c(Y=rbinom(size = DI,prob = sym_prob,n = 1))
}

reed_frost.meas.dens <- function(Y,DI,sym_prob,log,t,...){
  dbinom(x = Y,size = DI,prob = sym_prob,log=log)
}

reed_frost.param.transform <- parameter_trans(log=c("h0","beta"),
                                              logit=c("gamma","A","sym_prob","pI"))

reedfrost <- pomp(t0=0,  
                  data=data.frame(time=1:100,Y=NA),
                  time="time",
                  rinit=reed_frost.init.sim,
                  rprocess=discrete_time(reed_frost.proc.sim,delta.t=1),
                  rmeasure = reed_frost.meas.sim,
                  dmeasure = reed_frost.meas.dens,
                  obsnames = c("Y"),
                  statenames = c("E","DI","H","S"),
                  paramnames = c("gamma","h0","beta","A","pI","s0","sym_prob"),
                  partrans = reed_frost.param.transform
)

true_params <- c(gamma=0.5,
  h0=1,
  beta=0.1,
  A=0.99,
  pI=1/3,
  s0=5,
  sym_prob=0.9)

sim <- simulate(reedfrost,params=true_params,nsim = 10000)

test_data <- simulate(reedfrost,params=true_params,nsim = 100)

# Implementing the basic particle filter

pf <- pfilter(sim,Np=1000)
logLik(pf)

test_params <- c(gamma=0.4,
                 h0=1,
                 beta=0.1,
                 A=0.99,
                 pI=1/3,
                 s0=5,
                 sym_prob=0.5)

test_pf <- pfilter(sim,params=test_params,Np=1000)
logLik(test_pf)

# Using Iterated Filtering to estimate parameters

res_filtrage <- mif2(data=sim,
     Nmif=100,
     params=test_params,
     Np=1000,
     cooling.fraction=0.7,
     rw.sd=rw.sd(gamma=0.02,h0=0.02,beta=0.02,A=0.02,pI=0.02,sym_prob=0.02))


