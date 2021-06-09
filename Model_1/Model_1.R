# Dynamic optimisation model with look-ahead policy

# Packages
library(tidyverse)
library(reshape2)
library(utils)

# Functions
source("Model_1/functions_1.R")
source("interpolate.R")

########################################## Look-ahead function ##########################################

# look-ahead function returns the decision which maximises fitness for the end of the planning period

look.ahead <- function(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, nh.critical, mh, ms, t, H, h, N, m, i, j, k, l){ # The argument h has to be the same value as H
  if (h == 0){
    return(fitness(nh, ns, nh.reproduction, ns.reproduction, t, m)) # if every step of the planning period has occurred, return the immediate contribution to fitness
    }
  else {
    W. = fitness(nh, ns, nh.reproduction, ns.reproduction, t, m) # Immediate contribution to fitness
    W <- array(data = NA, dim = ((2*N)+1)) # create an array to be populated by reproductive values
    for (d in 1:((2*N)+1)){
      nh. = as.numeric(new.state(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d)[1])
      ns. = as.numeric(new.state(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d)[2])
      s. = as.numeric(new.state(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d)[3])
      W[d] = W. + survival(nh., nh.critical)*look.ahead(nh., ns., s., nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, nh.critical, mh, ms,
                                                        (t+1), H, (h-1), N, m, i, j, k, l)
      }# recursive function, the next time step involves the same process
    if (h < H){
      return(max(W)) # if the time step is not the first, return the maximum reproductive value for all allocation decisions
      }
    else {
      return(which.max(W)) # return the index for the decision of the first time step
    }
  }
}

########################################## Optimisation array ##########################################

opt <- function(nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, nh.critical, mh, ms, T, H, h, N, m, i, j, k, l) {
pb <- txtProgressBar(min = 0, max = ns.max*nh.max*s.max*(T-1), style = 3)
d.opt <- array(data = NA, dim = c((T-1), nh.max, ns.max, s.max)) # array to be populated with best decisions
for (t in 1:(T-1)){
  for (nh in 1:nh.max){
    for (ns in 1:ns.max){
      for (s in 1:s.max){
        if (t <= (T-H)){
          d.opt[t, nh, ns, s] <- look.ahead(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, nh.critical, mh, ms, t, H, h, N, m, i, j, k, l) # look-ahead model with the full planning period
          setTxtProgressBar(pb, (s-1) + (ns-1)*s.max + (nh-1)*ns.max*s.max + (t-1)*nh.max*ns.max*s.max)
          }
        else {
          d.opt[t, nh, ns, s] <- look.ahead(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, nh.critical, mh, ms, t, (T-t), (T-t), N, m, i, j, k, l) # look-ahead model with planning period of the remaining time steps
          setTxtProgressBar(pb, (s-1) + (ns-1)*s.max + (nh-1)*ns.max*s.max + (t-1)*nh.max*ns.max*s.max)
          }
        }
      }
    }
  }
close(pb)
return(d.opt)
}

####################################### Optimisation array plot #######################################

opt.plot <- function(d.opt){
  d.opt = melt(d.opt)
  names(d.opt) <- c("time.step", "Nh", "Ns", "sym.pop", "decision")
  plot <- ggplot(data = d.opt, aes(time.step, sym.pop, fill = decision)) +
    geom_tile()+
    facet_grid(Nh ~ Ns)+
    scale_fill_gradient(low = "#80CBC4", high = "#00695C")+
    theme(axis.title.y.right = element_text("Nh"),
          axis.title.x.top = element_text("Ns"))
return(plot)
}

############################################ Simulations #############################################

sim <- function(d.opt, nh.max, ns.max, sm.max, nh.reproduction, ns.reproduction, nh.critical, mh, ms, T, N, i, j, k, l, I){
sim.plot <- ggplot(NULL, aes(x = time.step, y = sym.pop)) # empty plot to be populated with simulation runs
for (i in 1:I){
sim <- data.frame(time.step = integer(), nh.reserves = integer(), ns.reserves = integer(), sym.pop = integer(), decision = integer())%>% # empty data frame to be populated with state values and allocations
add_row(time.step = 1,  nh.reserves = chop(round(rnorm(1, mean = (nh.max/3), sd = 0.5)), 0, nh.max), ns.reserves = chop(round(rnorm(1, mean = (ns.max/2), sd = 0.5)), 0, ns.max),
        sym.pop = chop(round(rnorm(1, mean = (s.max/3), sd = 0.5)), 0, s.max), decision = NA) # initial state values chosen at random
for (t in 1:(T-1)){
  if (sim$nh.reserves[t] >= nh.critical & sim$ns.reserves[t] >= 1 & sim$sym.pop[t] >= 1){
    d = interpolate(d.opt, sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], t) # extract optimum allocation amount based on state values
    sim$decision[t] = d
    nh. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d)[1])
    ns. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d)[2])
    s. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d)[3])
    sim <- add_row(sim, time.step = (t+1), sym.pop = s., nh.reserves = nh., ns.reserves = ns., decision = NA) # populate data frame with new state values
    }
  else {
    break # if reserves get to zero, stop the simulation
  }
  }
sim.plot <- sim.plot+
geom_line(data = sim, aes(x = time.step, y = sym.pop, group = 1)) # add data frame to plot
}
sim.plot <- sim.plot+
  labs(x = "time step", y = "symbiont density")
return(sim.plot)
}

########################################### Model runs ###############################################

d.opt.1a <- opt(10, 5, 10, 4, 2, 1, 1, 1, 15, 4, 4, 5, 1, 0.2, 2, 2, 0.5)
opt.plot.1a <- opt.plot(d.opt.1a)
sim.1a <- sim(10, 5, 10, 4, 2, 1, 1, 1, 15, 5, 0.2, 2, 2, 0.5, 100)

save(d.opt.1a, file = "Model_1/d.opt.1a.RData")
save(opt.plot.1a, file = "Model_1/opt.plot.1a.RData")
save(sim.1a, file = "Model_1/sim.1a.RData")
