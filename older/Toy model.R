# Dynamic optimisation model with look-ahead policy

# Variables and parameters:

# s = symbiont population
# nh = host/diet-derived nutrient reserves
# ns = symbiont-derived nutrient reserves
s.max <- 3 # maximum symbiont population
nh.max <- 3 # maximum reserves of host-dervied nutrients
ns.max <- 3 # maximum reserves of symbiont-derived nutrients
N <- 2 # nutritional resources available for allocation, per time step
mh <- 1 # metabolic expenditure from host-derived reserves, per time step
ms <- 1 # metabolic expenditure from symbiont-derived reserves, per time step
T <- 8 # time horizon
H <- 2 # planing horizon

# Ecological functions
source("C:/Users/apico/Documents/Mathilda/PhD/Notes/functions.R")
source("C:/Users/apico/Documents/Mathilda/PhD/Notes/runge-kutta.R")
# functions: chop(x, x.min, x.max); provisioning(s, a); population.change(a); fitness(nh, ns); rk4(derivs, sym.pop, a, n)

########################################## Look-ahead function ##########################################

look.ahead <- function(s, nh, ns, N, t, H){ # The argument t is the same value as H
  if (t == -1){
    return(fitness(nh, ns)) # if every step of the planning period has occured, return the fitness
  } else if (t > -1){
    W <- array(data = NA, dim = (N+1)) # create an array to be populated by the fitness values
    for (a in 0:N){
      nh. = chop(nh-mh+N-a, 0, nh.max)
      ns. = chop(ns-ms+provisioning(s, a), 0, ns.max)
      # s. = chop(rk4(derivs, s, a, 100), 0, s.max)
      s. = chop(population.change(a)*s, 0, s.max)
      W[(a+1)] = look.ahead(s., nh., ns., N, (t-1), H) # recursive function, the next time step involves the same process
    }
    if (t < H){
      return(max(W)) # if the time step is not the first, return the maximum fitness for all allocation decisions
    } else if (t == H){
      return(which.max(W)-1) # return the allocation decision of the first time step
    }
  }
}

########################################## Optimisation array ##########################################

library(reshape2)

a.opt <- array(data = NA, dim = c((T-1), s.max, nh.max, ns.max)) # array to be populated with best decisions
for (t in 1:(T-1)){
  for (s in 1:s.max){
    for (nh in 1:nh.max){
      for (ns in 1:ns.max){
        if (t < (T-H)){
          a.opt[t, s, nh, ns] <- look.ahead(s, nh, ns, N, H, H)} # look-ahead model with the full planning period
        else if (t >= (T-H)){
          a.opt[t, s, nh, ns] <- look.ahead(s, nh, ns, N, (T-t-1), (T-t-1))} # look-ahead model with planning period of the remaining time steps
      }
    }
  }
}

a.opt.figure <- melt(a.opt)
names(a.opt.figure) <- c("time.step", "sym.pop", "Nh", "Ns", "allocation")
a.opt.plot <- ggplot(data = a.opt.figure, aes(time.step, sym.pop, fill = allocation)) +
  geom_tile()+
  facet_grid(Nh ~ Ns)+
  scale_fill_gradient(low = "#80CBC4", high = "#00695C")

a.opt.plot

############################################ Simulations #############################################

library(tidyverse)
library(oce)
library(ggplot2)

sim.plot <- ggplot(NULL, aes(x = time.step, y = sym.pop, group = 1)) # empty plot to be populated with simulation runs

for (i in 1:20){
  sim <- data.frame(time.step = integer(), sym.pop = integer(), nh.reserves = integer(),
                    ns.reserves = integer(), allocation = integer()) %>% # empty data frame to be populated with state values and allocations
    add_row(time.step = 1, sym.pop = round(rnorm(1, mean = 1.5, sd = 0.5)), nh.reserves = round(rnorm(1, mean = 1.5, sd = 0.5)),
            ns.reserves = round(rnorm(1, mean = 1.5, sd = 0.5)), allocation = NA) # initial state values chosen at random
  for (t in 1:(T-1)){
    if (fitness(sim$nh.reserves[t], sim$ns.reserves[t]) >= 1 && sim$sym.pop[t] >= 1){
      a = a.opt[t, sim$sym.pop[t], sim$nh.reserves[t], sim$ns.reserves[t]]
      nh. = chop(sim$nh.reserves[t]-mh+N-a, 0, nh.max)
      ns. = chop(sim$ns.reserves[t]-ms+provisioning(s, a), 0, ns.max)
      s. = chop(rk4(derivs, sim$sym.pop[t], a, 100), 0, s.max)  # rk4: runge-kutta function numerically solves differential equation
      sim <- add_row(sim, time.step = (t+1), sym.pop = s., nh.reserves = nh., ns.reserves = ns.,
                     allocation = NA) # populate data frame with new state values
    } else {
      break # if fitness or symbiont population becomes zero, stop the simulation
    }
  }
  print(i)
  print(sim)
  sim.plot <- sim.plot+
    geom_line(data = sim, aes(x = time.step, y = sym.pop)) # add data frame to plot
}

sim.plot