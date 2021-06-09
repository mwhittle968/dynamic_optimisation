# Dynamic optimisation model with look-ahead policy

# Packages
library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
library(parallel)

# Variables and parameters:
# s = symbiont population
# nh = host/diet-derived nutrient reserves
# ns = symbiont-derived nutrient reserves
s.max <- 10 # maximum symbiont population
nh.max <- 20 # maximum reserves of host-dervied nutrients
ns.max <- 5 # maximum reserves of symbiont-derived nutrients
nh.critical <- 6 # minimum reserves needed to reproduce
ns.critical <- 3 # minimum reserves needed to reproduce
N <- 10 # nutritional resources available for allocation, per time step
mh <- 3 # metabolic expenditure from host-derived reserves, per time step
ms <- 2 # metabolic expenditure from symbiont-derived reserves, per time step
T <- 10 # time horizon
H <- 4 # planning horizon
n <- 20 # time steps for runge-kutta

# Functions
source("C:/Users/apico/Documents/Mathilda/PhD/Notes/R_code/functions.R")
source("C:/Users/apico/Documents/Mathilda/PhD/Notes/R_code/interpolate.R")

########################################## Look-ahead function ##########################################

# look-ahead function returns the decision which maximises fitness for the end of the planning period

look.ahead <- function(s, nh, ns, t, H, h){ # The argument h has to be the same value as H
  if (h == 0){
    return(fitness(nh, ns, t)) # if every step of the planning period has occured, return the immediate contribution to fitness
  } else if (h > 0){
    W. = fitness(nh, ns, t) # Immediate contribution to fitness
    S = ex.survival*in.survival(nh)
    W <- array(data = NA, dim = (N+1)) # create an array to be populated by the fitness values
    for (a in 0:N){
      nh. = as.numeric(new.state(t, a, s, ns, nh)[1])
      ns. = as.numeric(new.state(t, a, s, ns, nh)[2])
      s. = as.numeric(new.state(t, a, s, ns, nh)[3])
      W[(a+1)] = W. + S*look.ahead(s., nh., ns., (t+1), H, (h-1)) # recursive function, the next time step involves the same process
    }
    if (h < H){
      return(max(W)) # if the time step is not the first, return the maximum fitness for all allocation decisions + the immediate contribution to fitness
    } else if (h == H){
      return(which.max(W)-1) # return the allocation decision of the first time step
    }
  }
}

########################################## Optimisation array ##########################################

a.opt.fun <- function(){
  a.opt <- array(data = NA, dim = c((T-1), s.max, nh.max, ns.max)) # array to be populated with best decisions
for (t in 1:(T-1)){
  for (s in 1:s.max){
    for (nh in 1:nh.max){
      for (ns in 1:ns.max){
        if (t <= (T-H)){
          a.opt[t, s, nh, ns] <- look.ahead(s, nh, ns, t, H, H) # look-ahead model with the full planning period
          print("t, s, nh, ns =")
          print(t)
          print(s)
          print(nh)
          print(ns)}
        else if (t > (T-H)){
          a.opt[t, s, nh, ns] <- look.ahead(s, nh, ns, t, (T-t), (T-t)) # look-ahead model with planning period of the remaining time steps
        print("t, s, nh, ns =")
        print(t)
        print(s)
        print(nh)
        print(ns)
        }
      }
      }
    }
}
  return(a.opt)
}

system.time({a.opt <- parLapply(a.opt.fun(), mc.cores = numcores)})

a.opt.figure <- melt(a.opt)
names(a.opt.figure) <- c("time.step", "sym.pop", "Nh", "Ns", "allocation")
a.opt.plot <- ggplot(data = a.opt.figure, aes(time.step, sym.pop, fill = allocation)) +
  geom_tile()+
  facet_grid(Nh ~ Ns)+
  scale_fill_gradient(low = "#80CBC4", high = "#00695C")
a.opt.plot

############################################ Simulations #############################################

sim.plot <- ggplot(NULL, aes(x = time.step, y = sym.pop)) # empty plot to be populated with simulation runs

for (i in 1:500){
sim <- data.frame(time.step = integer(), sym.pop = integer(), nh.reserves = integer(),
                  ns.reserves = integer(), allocation = integer())%>% # empty data frame to be populated with state values and allocations
add_row(time.step = 1, sym.pop = round(rnorm(1, mean = 5, sd = 1)), nh.reserves = round(rnorm(1, mean = 10, sd = 1)),
        ns.reserves = round(rnorm(1, mean = 2, sd = 1)), allocation = NA) # initial state values chosen at random
    for (t in 1:(T-1)){
      if (sim$nh.reserves[t] >= 1 & sim$ns.reserves[t] >= 1){
        a = interpolate(t, sim$sym.pop[t], sim$nh.reserves[t], sim$ns.reserves[t]) # extract optimum allocation amount based on state values
        sim$allocation[t] = a
        nh. = as.numeric(new.state(t, a, sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t])[1])
        ns. = as.numeric(new.state(t, a, sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t])[2])
        s. = as.numeric(new.state(t, a, sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t])[3])
        sim <- add_row(sim, time.step = (t+1), sym.pop = s., nh.reserves = nh., ns.reserves = ns.,
                       allocation = NA) # populate data frame with new state values
        } else {
          break # if reserves get to zero, stop the simulation
        }
    }
print(i)
print(sim)
sim.plot <- sim.plot+
geom_line(data = sim, aes(x = time.step, y = sym.pop)) # add data frame to plot
}

sim.plot

