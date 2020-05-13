library(tidyverse)
library(pomp)
library(panelPomp)

# First, specify the simulation function

reed_frost.proc.sim <- function(gamma,h0,beta,A,pI,pS,S,E,IP,IS,H,...){
  infection_prob <- 1-A*exp(-beta*(H+IP*h0))
  new_infected_pre <- rbinom(size = E,prob = pI,n=1)
  new_infected <- rbinom(size=IP,prob=pS,n=1)
  new_exposed <- rbinom(size = S,prob = infection_prob,n=1)
  c(E=E+new_exposed-new_infected_pre,
    IP=IP+new_infected_pre-new_infected,
    IS=IS+new_infected,
    H=gamma*H+h0*new_infected,
    S=S-new_exposed)
}

reed_frost.init.sim <- function(s0,...){
  c(S=s0,E=0,IP=0,IS=0,H=0)
}

reed_frost.meas.sim <- function(IS,sym_prob,...){
  c(Y=rbinom(size = IS,prob = sym_prob,n = 1))
}

reed_frost.meas.dens <- function(Y,IS,sym_prob,log,t,...){
  dbinom(x = Y,size = IS,prob = sym_prob,log=log)
}

reed_frost.param.transform <- parameter_trans(log=c("h0","beta"),
                                              logit=c("gamma","A","sym_prob","pI","pS"))

reedfrost <- pomp(t0=0,  
                  data=data.frame(time=1:100,Y=NA),
                  time="time",
                  rinit=reed_frost.init.sim,
                  rprocess=discrete_time(reed_frost.proc.sim,delta.t=1),
                  rmeasure = reed_frost.meas.sim,
                  dmeasure = reed_frost.meas.dens,
                  obsnames = c("Y"),
                  statenames = c("E","IP","IS","H","S"),
                  paramnames = c("gamma","h0","beta","A","pI","pS","s0","sym_prob"),
                  partrans = reed_frost.param.transform,
                  accumvars = "IS"
)

true_params <- c(gamma=0.5,
  h0=1,
  beta=0.1,
  A=0.99,
  pI=1/3,
  pS=1/2,
  s0=5,
  sym_prob=0.9)

sim <- simulate(reedfrost,params=true_params,nsim = 10000)

test_data_1 <- simulate(reedfrost,params=true_params,nsim = 1)
test_data_2 <- simulate(reedfrost,params=true_params,nsim = 1)


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

res_filtrage <- mif2(data=test_data_1,
     Nmif=100,
     Np=1000,
     cooling.fraction=0.7,
     rw.sd=rw.sd(gamma=0.02,h0=0.02,beta=0.02,A=0.02,pI=0.02,sym_prob=0.02))

# What about panelPOMP?

model_list <- list(u1 = test_data_1,u2=test_data_2)
panel_reed_frost <- panelPomp(object = model_list,shared=true_params)

test_ppf <- pfilter(panel_reed_frost,Np=1000)

res_filtrage <- mif2(data=panel_reed_frost,
                     shared.start=true_params,
                     Nmif=100,
                     Np=1000,
                     cooling.fraction.50=0.7,
                     cooling.type="geometric",
                     rw.sd=rw.sd(gamma=0.02,h0=0.02,beta=0.02,A=0.02,pI=0.02,sym_prob=0.02))
