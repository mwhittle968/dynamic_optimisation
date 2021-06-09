# Dynamic optimisation model with look-ahead policy

# Packages
library(tidyverse)
library(reshape2)
library(utils)

# Functions
source("Model_3/functions_3.R")
source("interpolate.R")

########################################## Data frame for events and parameters ##########################################

time.step <- c(1:61)
female <- as.data.frame(time.step)%>%
  mutate(stage = ifelse(time.step <= 13, "pupae", "adult"))%>%
  mutate(reproduction = ifelse(time.step %in% seq(from = 21, to = 61, by = 4), "yes", "no"))%>%
  mutate(n = ifelse(stage == "pupae", 0, 5))%>%
  mutate(mh = ifelse(stage == "pupae", 0.2, 1), ms = ifelse(stage == "pupae", 0.2, 1))%>%
  mutate(ex.survival = ifelse(stage == "pupae", 1, (1 - 1/61)))

########################################## Dynamic Programming with for loops ##########################################

opt <- function(nh.max = 10, ns.max = 5, s.max = 10, nh.repro = 4, ns.repro = 2, nh.crit = 1, T = 61, N = 5, m = 1, i, j, k, l){
  pb <- txtProgressBar(min = 0, max = ns.max*nh.max*s.max*(T-1), style = 3)
  V <- array(data = NA, dim = c((T), nh.max, ns.max, s.max)) # Empty fitness array, to be populated with values from T back to 1
  d.opt <- array(data = NA, dim = c((T-1), nh.max, ns.max, s.max)) # Empty array for best decisions at each state and time
  V[T,,,] <- 0 # terminal fitness function
  for (t in (T-1):1){ # Iterates backwards in time
    for (nh in 1:nh.max){
      for (ns in 1:ns.max){
        for (s in 1:s.max){
          H = array(data = NA, dim = (2*N+1)) # Empty reproductive value array (for one state)
          for (d in 1:(2*N+1)){
            nh. = as.numeric(new.state(nh, ns, s, nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d)[1]) # Calculate new values for states
            ns. = as.numeric(new.state(nh, ns, s, nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d)[2])
            s. = as.numeric(new.state(nh, ns, s, nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d)[3])
            H[d] = B(nh, ns, nh.repro, ns.repro, t, m) + S(nh., nh.crit, t)*interpolate(V, chop(nh., 1, nh.max), chop(ns., 1, ns.max),
            chop(s., 1, s.max), (t+1)) # Calculate reproductive value given decision
            }
          V[t, nh, ns, s] <- max(H) # define fitness as maximum of the reproductive values
          d.opt[t, nh, ns, s] <- which.max(H) # Set best decision as the one which maximises fitness
          setTxtProgressBar(pb, (s-1) + (ns-1)*s.max + (nh-1)*ns.max*s.max + (T-t-1)*nh.max*ns.max*s.max)
        }
      }
    }
  }
  return(d.opt)
}

####################################### Optimisation array plot ######################################

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

sim <- function(data, nh.max = 10, ns.max = 5, s.max = 10, nh.repro = 4, ns.repro = 2, nh.crit = 1, mh = 1, ms = 1, T = 61, N = 5, i, j, k, l, I){
sim.plot <- ggplot(NULL, aes(x = time.step, y = sym.pop)) # empty plot to be populated with simulation runs
for (i in 1:I){
sim <- data.frame(time.step = integer(), nh.reserves = integer(), ns.reserves = integer(), sym.pop = integer(), decision = integer())%>% # empty data frame to be populated with state values and allocations
add_row(time.step = 1,  nh.reserves = chop(round(rnorm(1, mean = (nh.max), sd = 0.5)), 0, nh.max), ns.reserves = chop(round(rnorm(1, mean = (ns.max/2), sd = 0.5)), 0, ns.max),
        sym.pop = chop(round(rnorm(1, mean = (s.max/3), sd = 0.5)), 0, s.max), decision = NA) # initial state values chosen at random
for (t in 1:(T-1)){
  if (sim$nh.reserves[t] >= nh.crit & sim$ns.reserves[t] >= 1 & sim$sym.pop[t] >= 1){
    d = interpolate(data, sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], t) # extract optimum allocation amount based on state values
    sim$decision[t] = d
    nh. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d)[1])
    ns. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d)[2])
    s. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d)[3])
    sim <- add_row(sim, time.step = (t+1), nh.reserves = nh., ns.reserves = ns., sym.pop = s., decision = NA) # populate data frame with new state values
    print(t)}
  else {
    break # if reserves get to zero, stop the simulation
  }
}
print(sim)
sim.plot <- sim.plot+
geom_line(data = sim, aes(x = time.step, y = sym.pop, group = 1)) # add data frame to plot
}
sim.plot <- sim.plot+
  labs(x = "time step", y = "symbiont density")
return(sim.plot)
}

########################################### Model runs ###############################################

d.opt.3a <- opt(i = 2, j = 3, k = 3, l = 2)
opt.plot.3a <- opt.plot(d.opt.3a)
opt.plot.3a
sim.3a <- sim(d.opt.3a, i = 2, j = 3, k = 3, l = 2, I = 100)
sim.3a

save(d.opt.3a, file = "Model_3/d.opt.3a.RData")
save(opt.plot.3a, file = "Model_3/opt.plot.3a.RData")
save(sim.3a, file = "Model_3/sim.3a.RData")
