library(tidyverse)
library(pomp)

simulate(t0=0, times=1:100,
    params = c(gamma=0.8,
              h0=1,
              beta=1,
              A=0.95,
              pI=1/3,
              s0=5,
              sym_prob=0.5),
    rinit=function(s0,...){
      c(E=0,DI=0,H=0,S=5)
    },
    rprocess=discrete_time(
      function(gamma,h0,beta,A,pI,S,E,DI,H,...){
        infection_prob <- 1-A*exp(-beta*H)
        new_infected <- rbinom(size = E,prob = pI,n=1)
        new_exposed <- rbinom(size = S,prob = infection_prob,n=1)
        c(E=E+new_exposed-new_infected,
          DI=new_infected,
          H=gamma*H+h0*DI,
          S=S-new_exposed)
      },
      delta.t=1
    ),
    rmeasure = function(DI,sym_prob,...){
      c(Y=rbinom(size = DI,prob = sym_prob,n = 1))
    }
) -> sim1