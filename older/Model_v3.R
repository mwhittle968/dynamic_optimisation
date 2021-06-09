# Dynamic optimisation model with look-ahead policy

# Packages
library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
#library(parallel)
library(utils)

# Variables and parameters:
# s = symbiont population
# nh = host/diet-derived nutrient reserves
# ns = symbiont-derived nutrient reserves
s.max <- 5 # maximum symbiont population
nh.max <- 5 # maximum reserves of host-dervied nutrients
ns.max <- 5 # maximum reserves of symbiont-derived nutrients
nh.critical <- 2 # minimum reserves needed to reproduce
ns.critical <- 1 # minimum reserves needed to reproduce
N <- 4 # nutritional resources available for allocation, per time step
mh <- 1 # metabolic expenditure from host-derived reserves, per time step
ms <- 2 # metabolic expenditure from symbiont-derived reserves, per time step
T <- 6 # time horizon
H <- 3 # planning horizon

# Functions
source("C:/Users/apico/Documents/Mathilda/PhD/Notes/R_code/functions_v3.R")
source("C:/Users/apico/Documents/Mathilda/PhD/Notes/R_code/interpolate_v3.R")

################################################ Decision ###############################################

d <- c(0:(2*N)) # vector containing possible decisions, index 1 is no investment or regulation, indices 2:N+1 are investment, N+2:2N+1 are regulation

########################################## Look-ahead function ##########################################

# look-ahead function returns the decision which maximises fitness for the end of the planning period

look.ahead <- function(s, nh, ns, t, H, h){ # The argument h has to be the same value as H
  if (h == 0){
    return(fitness(nh, ns, t)) # if every step of the planning period has occured, return the immediate contribution to fitness
  } else if (h > 0){
    W. = fitness(nh, ns, t) # Immediate contribution to fitness
    W <- array(data = NA, dim = ((2*N)+1)) # create an array to be populated by the fitness values
    for (d in 1:((2*N)+1)){
      nh. = as.numeric(new.state(t, d, s, ns, nh)[1])
      ns. = as.numeric(new.state(t, d, s, ns, nh)[2])
      s. = as.numeric(new.state(t, d, s, ns, nh)[3])
      S = ex.survival*in.survival(nh.)
      W[(d)] = W. + S*look.ahead(s., nh., ns., (t+1), H, (h-1)) # recursive function, the next time step involves the same process
      }
    if (h < H){
      return(max(W)) # if the time step is not the first, return the maximum fitness for all allocation decisions + the immediate contribution to fitness
    } else if (h == H){
      return(which.max(W)) # return the index for the decision of the first time step
    }
  }
}

########################################## Optimisation array ##########################################

#d.opt.fun <- function(){
pb <- txtProgressBar(min = 0, max = 100, style = 3)
d.opt <- array(data = NA, dim = c((T-1), s.max, nh.max, ns.max)) # array to be populated with best decisions
for (t in 1:(T-1)){
  for (s in 1:s.max){
    for (nh in 1:nh.max){
      for (ns in 1:ns.max){
        if (t <= (T-H)){
          setTxtProgressBar(pb, i)
          d.opt[t, s, nh, ns] <- look.ahead(s, nh, ns, t, H, H) # look-ahead model with the full planning period
          }
        else if (t > (T-H)){
          setTxtProgressBar(pb, i)
          d.opt[t, s, nh, ns] <- look.ahead(s, nh, ns, t, (T-t), (T-t)) # look-ahead model with planning period of the remaining time steps
        }
      
      }
      }
    }
}

close(pb)
#  return(d.opt)
#}

#system.time({a.opt <- parLapply(a.opt.fun(), mc.cores = numcores)})

d.opt.figure <- melt(d.opt)
names(d.opt.figure) <- c("time.step", "sym.pop", "Nh", "Ns", "decision")
d.opt.plot <- ggplot(data = d.opt.figure, aes(time.step, sym.pop, fill = decision)) +
  geom_tile()+
  facet_grid(Nh ~ Ns)+
  scale_fill_gradient(low = "#80CBC4", high = "#00695C")+
  theme(axis.title.y.right = element_text("Nh"),
        axis.title.x.top = element_text("Ns"))
d.opt.plot

############################################ Simulations #############################################

sim.plot <- ggplot(NULL, aes(x = time.step, y = sym.pop)) # empty plot to be populated with simulation runs

for (i in 1:1000){
sim <- data.frame(time.step = integer(), sym.pop = integer(), nh.reserves = integer(),
                  ns.reserves = integer(), allocation = integer())%>% # empty data frame to be populated with state values and allocations
add_row(time.step = 1, sym.pop = chop(round(rnorm(1, mean = (s.max/3), sd = 0.5)), 0, s.max), nh.reserves = chop(round(rnorm(1, mean = (nh.max/3), sd = 0.5)), 0, nh.max),
        ns.reserves = chop(round(rnorm(1, mean = (ns.max/2), sd = 0.5)), 0, ns.max), allocation = NA) # initial state values chosen at random
    for (t in 1:(T-1)){
      if (sim$nh.reserves[t] >= 1 & sim$ns.reserves[t] >= 1 & sim$sym.pop[t] >= 1){
        d = interpolate(t, sim$sym.pop[t], sim$nh.reserves[t], sim$ns.reserves[t]) # extract optimum allocation amount based on state values
        sim$allocation[t] = d
        nh. = as.numeric(new.state(t, d, sim$sym.pop[t], sim$ns.reserves[t], sim$nh.reserves[t])[1])
        ns. = as.numeric(new.state(t, d, sim$sym.pop[t], sim$ns.reserves[t], sim$nh.reserves[t])[2])
        s. = as.numeric(new.state(t, d, sim$sym.pop[t], sim$ns.reserves[t], sim$nh.reserves[t])[3])
        sim <- add_row(sim, time.step = (t+1), sym.pop = s., nh.reserves = nh., ns.reserves = ns.,
                       allocation = NA) # populate data frame with new state values
        } else {
          break # if reserves get to zero, stop the simulation
        }
    }
print(i)
print(sim)
sim.plot <- sim.plot+
geom_line(data = sim, aes(x = time.step, y = sym.pop, group = 1)) # add data frame to plot
}

sim.plot <- sim.plot+
  labs(x = "time step", y = "symbiont density")+
  ylim(0, 5)

sim.plot

