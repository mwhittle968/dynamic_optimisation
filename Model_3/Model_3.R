# Dynamic optimisation model with look-ahead policy

# Packages
library(tidyverse)
library(reshape2)
library(utils)

# Functions
source("Model_3/functions_3.R")
source("interpolate.R")
source("Model_3/test_parameters.R")

###################################################################################################################################################

dyanamic.optimisation <- function(model.name, nh.max, ns.max, s.max, nh.repro, ns.repro, s.repro, nh.crit.pupae, nh.crit.adult, nh.larva, ns.larva,
                                  mh.pupae, mh.adult, ms.pupae, ms.adult, ex.surv.pupae, ex.surv.adult, T, n, m, i, j, k, l, Q){
  
###################################################### Data frames for events and parameters ######################################################

parameters <- c(model.name, nh.max, ns.max, s.max, nh.repro, ns.repro, s.repro, nh.crit.pupae, nh.crit.adult, nh.larva, ns.larva,
                mh.pupae, mh.adult, ms.pupae, ms.adult, ex.surv.pupae, ex.surv.adult, T, n, m, i, j, k, l, Q)
parameters <- t(as.data.frame(parameters))
colnames(parameters) <- c("model.name", "nh.max", "ns.max", "s.max", "nh.repro", "ns.repro", "s.repro", "nh.crit.pupae", "nh.crit.adult", "nh.larva", "ns.larva",
                     "mh.pupae", "mh.adult", "ms.pupae", "ms.adult", "ex.surv.pupae", "ex.surv.adult", "T", "n", "m", "i", "j", "k", "l", "Q")
save(parameters, file = paste("Model_3/parameters_", model.name, ".RData", sep = ""))

time.step <- c(1:61)
F <- as.data.frame(time.step)%>%
  mutate(stage = ifelse(time.step <= 13, "pupae", "adult"))%>%
  mutate(repro = ifelse(time.step %in% seq(from = 21, to = 61, by = 4), "yes", "no"))%>%
  mutate(N = ifelse(stage == "pupae", 0, n))%>%
  mutate(mh = ifelse(stage == "pupae", mh.pupae, mh.adult), ms = ifelse(stage == "pupae", ms.pupae, ms.adult))%>%
  mutate(ex.surv = ifelse(stage == "pupae", ex.surv.pupae, ex.surv.adult))%>%
  mutate(nh.crit = ifelse(stage == "pupae", nh.crit.pupae, nh.crit.adult))

######################################################## Dynamic Programming with for loops #######################################################

pb <- txtProgressBar(min = 0, max = ns.max*nh.max*s.max*(T-1), style = 3)
V <- array(data = NA, dim = c((T), nh.max, ns.max, s.max)) # Empty fitness array, to be populated with values from T back to 1
d.opt <- array(data = NA, dim = c((T-1), nh.max, ns.max, s.max)) # Empty array for best decisions at each state and time
V[T,,,] <- 0 # terminal fitness function
for (t in (T-1):1){ # Iterates backwards in time
  for (nh in 1:nh.max){
    for (ns in 1:ns.max){
      for (s in 1:s.max){
        H = array(data = NA, dim = (2*n+1)) # Empty reproductive value array (for one state)
        for (d in 1:(2*n+1)){
          nh. = as.numeric(new.state(nh, ns, s, t, d)[1]) # Calculate new values for states
          ns. = as.numeric(new.state(nh, ns, s, t, d)[2])
          s. = as.numeric(new.state(nh, ns, s, t, d)[3])
          H[d] = B(nh, ns, s, t) + S(nh., t)*interpolate(V, chop(nh., 1, nh.max), chop(ns., 1, ns.max), chop(s., 1, s.max), (t+1))
          }
        V[t, nh, ns, s] <- max(H) # define fitness as maximum of the reproductive values
        d.opt[t, nh, ns, s] <- which.max(H) # Set best decision as the one which maximises fitness
        setTxtProgressBar(pb, (s-1) + (ns-1)*s.max + (nh-1)*ns.max*s.max + (T-t-1)*nh.max*ns.max*s.max)
      }
    }
  }
}
save(d.opt, file = paste("Model_3/d.opt_", model.name, ".RData", sep = ""))

############################################################# Optimisation array plot ############################################################

d.opt. <- melt(d.opt)
  names(d.opt.) <- c("time.step", "Nh", "Ns", "sym.pop", "decision")
  d.opt.plot <- ggplot(data = d.opt., aes(time.step, sym.pop, fill = decision)) +
    geom_tile()+
    facet_grid(Nh ~ Ns)+
    scale_fill_gradient(low = "#80CBC4", high = "#00695C")+
    theme(axis.title.y.right = element_text("Nh"),
          axis.title.x.top = element_text("Ns"))
save(d.opt.plot, file = paste("Model_3/d.opt.plot_", model.name, ".RData", sep = ""))

################################################################### Simulations ##################################################################

sim.plot <- ggplot(NULL, aes(x = time.step, y = sym.pop)) # empty plot to be populated with simulation runs
for (q in 1:Q){
sim <- data.frame(time.step = integer(), nh.reserves = integer(), ns.reserves = integer(), sym.pop = integer(), decision = integer())%>%
# empty data frame to be populated with state values and allocations
add_row(time.step = 1,  nh.reserves = chop(round(rnorm(1, mean = (nh.max), sd = 1)), 0, nh.max),
        ns.reserves = chop(round(rnorm(1, mean = (ns.max), sd = 1)), 0, ns.max),
        sym.pop = chop(round(rnorm(1, mean = (s.max/5), sd = 1)), 0, s.max), decision = NA) # initial state values chosen at random
for (t in 1:(T-1)){
  if (sim$nh.reserves[t] >= F$nh.crit[t]){
    d = interpolate(d.opt, chop(sim$nh.reserves[t], 1, nh.max), chop(sim$ns.reserves[t], 1, ns.max), chop(sim$sym.pop[t], 1, s.max), t)
    sim$decision[t] = d
    nh. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], t, d)[1])
    ns. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], t, d)[2])
    s. = as.numeric(new.state(sim$nh.reserves[t], sim$ns.reserves[t], sim$sym.pop[t], t, d)[3])
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
save(sim.plot, file = paste("Model_3/sim.plot_", model.name, ".RData", sep = ""))

################################################################################################################################################
}