
############################################### Dynamic optimisation model with dynamic programming ###############################################

# Required packages
library(dplyr)
library(reshape2)
library(utils)
library(pracma)
library(ggplot2)

# Set working directory
# setwd("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/dynamic_optimisation")

# Source code
source("bilinear_interpolation.R")
load("rio.female.RData")

###################################################################################################################################################

dynamic.optimisation <- function(model.version, model.name, nh.max = 20, s.max = 20, nh.repro = 7, s.repro = 6, nh.crit = 4, nh.transfer = 0.6,
                                 ex.surv = (1 - 1/48), T = 48, N = 100, n = 10, i = 0.1, j = 2, k = 2, Q = 100){

# model.version : character string for version of model e.g "6"
# model.name : character string for name of model
# nh.max : maximum level of nh / number of states of nh
# s.max : maximum level of s / number of states of s
# nh.repro : minimum level of nh for reproduction to occur
# s.repro : minimum level of s for an immediate contribution to fitness
# nh.crit : minimum level of nh for survival
# nh.transfer : proportion of nh transferred to larva / spent on reproduction
# nh.exp[t] : decrease in nh due to background metabolic expenditure, time dependant
# ex.surv : extrinsic (state-independent) probability of survival
# T : time horizon / host life expectancy
# N : number of resources gained during feeding  
# n : number of resources available for allocation by host
# m : parameter for fitness as a linear function of nh.larva
# i : cost of maintaining symbiont population (resources per s)
# j : amount of additional symbiont density (s per resource invested)
# k : amount of symbiont density removed (s per resource invested)
# Q : number of iterations for simulation
  
ptm <- proc.time()
  
########################################################## Folder to save model outputs #########################################################

dir.create(paste("Model_", model.version, sep = ""))

####################################################### Vector of reproductive events ############################################################

# Vector of time steps, metabolic expenditure, and if the female reproduces or not; yes = reproductive time step
# 1st reproductive event at t = 8, and every 4 time steps thereafter
repro <- c(1:T)
for (t in 1:T){
  repro[t] <- ifelse(t != 4 & t%%4 == 0, "yes", "no")
}
nh.exp <- c(1:T)
for (t in 1:T){
  nh.exp[t] <- ifelse(t <= 5, 2, ifelse(t <= 10, 1.5, 1))
}

#################################################################### Functions ###################################################################

# Chop function limits the range of values the state variables can take:
chop <- function(x, x.min, x.max){
  x = ifelse(x >= x.max, x.max, ifelse(x <= x.min, x.min, x))
  return(x)
}

# Probability of survival to next time step
survival <- function(nh){
  in.surv = ifelse(nh < nh.crit, 0, 1) # Below the threshold level of nh, probability of survival is zero
  S = ex.surv*in.surv
  return(S)
}

# Immediate contribution to fitness (B)
immediate.fitness <- function(t, nh, s){
  if (repro[t] == "no" | nh < nh.repro | s < s.repro){
    B = 0 # If it's a non-reproductive time step, or states are insufficient, the contribution to fitness is zero
  } else {# If it's a reproductive time step and the states are sufficient, the fitness is calculated by reserves transferred to larva
    nh.larva = (nh-nh.exp[t])*nh.transfer # After metabolic expenditure, female transfers some of her reserves and symbiont to larva
    nh.larva.min = (nh.repro-nh.exp[t])*nh.transfer
    Bh = 1/(1+exp(nh.larva.min-nh.larva)) # Fitness is logistic by nh.larva and s
    Bs = 1/(1+exp(s.repro-s))
    B = Bh*Bs
  }
  return(B)
}

# Maintenance, investment, regulation and production
maintenance <- function(s){
  m = s*i
  return(m)
}

increase <- function(s, d){
  s. = ifelse(s == 0, 0, s + (d-1)*j)
  return(s.)
}

decrease <- function(s, d){
  s. = ifelse(s == 0, 0, s - (d-n-1)*k)
  return(s.)
}

# New state variables:
new.state <- function(nh, s, t, d){
  # If it's a reproductive time step and states are sufficient, reserves are transferred to larva
  nh.larva = ifelse(repro[t] == "no" | nh < nh.repro, 0, (nh - nh.exp[t])*nh.transfer)
  if (d == 1){ 
    new.nh = chop(nh - nh.larva - nh.exp[t] + (N - maintenance(s))/20, 0, nh.max)
    new.s = s # If d = 1, there is no investment into changing s
  } else if (d > 1 & d <= (n+1)){ 
    new.nh = chop(nh - nh.exp[t] - nh.larva + (N - maintenance(s) - (d-1))/20, 0, nh.max)
    new.s = chop(increase(s, d), 0, s.max) # If d is between 2 and n+1, there is investment into increasing s
  } else { 
    new.nh = chop(nh - nh.exp[t] - nh.larva + (N - maintenance(s) - (d-n-1))/20, 0 , nh.max)
    new.s = chop(decrease(s, d), 0, s.max) # If d is greater than n+1, there is  investment into decreasing s
  }
  output = data.frame(new.nh, new.s)
  return(output)
}

######################################################## Dynamic Programming with for loops #######################################################

pb <- txtProgressBar(min = 0, max = nh.max*s.max*(T-1), style = 3)
reproductive.value <- array(data = NA, dim = c(T, nh.max, s.max)) # Empty array for reproductive values (V), to be populated from T back to 1
d.opt <- array(data = NA, dim = c((T-1), nh.max, s.max)) # Empty array for optimum decisions at each state and time
for (nh in 1:nh.max){
    for (s in 1:s.max){
      # Terminal fitness function is just the immediate contribution to fitness at final time step
      reproductive.value[T, nh, s] <- immediate.fitness(T, nh, s)
    }
  }
for (t in (T-1):1){ # Iterates backwards in time
  for (nh in 1:nh.max){
      for (s in 1:s.max){
        H = array(data = NA, dim = (2*n+1)) # Empty array for current reproductive value (H) (for one state)
        for (d in 1:(2*n+1)){
          # New states calculated for each decision
          nh. = as.numeric(new.state(nh, s, t, d)[1])
          s. = as.numeric(new.state(nh, s, t, d)[2])
          # Calculate current reproductive value (H) given each decision
          # Bilinear interpolation of reproductive.value array to find reproductive value (V) for the new states
          H[d] = immediate.fitness(t, nh, s) + survival(nh.)*interpolate(reproductive.value, chop(nh., 0, nh.max), chop(s., 0, s.max), (t+1))
        }
        reproductive.value[t, nh, s] <- max(H) # Reproductive value (V) at time t is maximum value of H
        d.opt[t, nh, s] <- which.max(H) # optimal decision is that which maximises current reproductive value (H)
       setTxtProgressBar(pb, (s-1) + (nh-1)*s.max + (T-t-1)*nh.max*s.max)
      }
  }
}

# Assign object name to array d.opt based on model name, and save
assign(paste("d.opt.", model.name, sep = ""), d.opt)
do.call(save, list(paste("d.opt.", model.name, sep = ""), file = paste("Model_", model.version, "/d.opt.", model.name, ".RData", sep = "")))

############################################################# Optimisation array plot ############################################################

# Visualisation of optimisation array
d.opt. <- melt(d.opt)
  names(d.opt.) <- c("time.step", "Nh", "sym.pop", "decision")
  d.opt.plot <- ggplot(data = d.opt., aes(time.step, sym.pop, fill = decision)) +
    geom_tile()+
    facet_wrap(~Nh)+
    scale_fill_gradient(low = "#80CBC4", high = "#00695C")#+
#    theme(axis.title.y.right = element_text("Nh"))
  
# Assign object name to plot based on model name, and save
assign(paste("d.opt.plot.", model.name, sep = ""), d.opt.plot)
do.call(save, list(paste("d.opt.plot.", model.name, sep = ""), file = paste("Model_", model.version, "/d.opt.plot.", model.name, ".RData", sep = "")))

################################################################## Simulations ##################################################################

# Run simulations of hosts using the optimal strategy to regulate symbiont density
sim.data <- data.frame(sim = integer(), time.step = integer(), nh.reserves = integer(), sym.pop = integer(), decision = integer())
for (q in 1:Q){
  sim.data <- add_row(sim.data, sim = q, time.step = 1, nh.reserves = chop(round(rnorm(1, mean = (nh.max), sd = 1)), 0, nh.max),
                      sym.pop = chop(round(rnorm(1, mean = (s.max/5), sd = 1)), 0, s.max), decision = NA)
  for (t in 1:(T-1)){
    sim.q <- filter(sim.data, sim == q)
    if (sim.q$nh.reserves[t] >= nh.crit){
    d = interpolate(d.opt, chop(sim.q$nh.reserves[t], 1, nh.max), chop(sim.q$sym.pop[t], 1, s.max), t)
    sim.data$decision[sim.data$sim == q & sim.data$time.step == t] <- d
    nh. = as.numeric(new.state(sim.q$nh.reserves[t], sim.q$sym.pop[t], t, d)[1])
    s. = as.numeric(new.state(sim.q$nh.reserves[t], sim.q$sym.pop[t], t, d)[2])
    sim.data <- add_row(sim.data, sim = q, time.step = (t+1), nh.reserves = nh., sym.pop = s., decision = NA)
    } else {
    break
  }
  }
}

# Assign object name to data frame based on model name, and save
assign(paste("sim.data.", model.name, sep = ""), sim.data)
do.call(save, list(paste("sim.data.", model.name, sep = ""), file = paste("Model_", model.version, "/sim.data.", model.name, ".RData", sep = "")))

# Plot of simulations
sim.plot <- ggplot(data = sim.data, aes(x = time.step, y = sym.pop, group = sim))+
  geom_line()+
  labs(x = "time step", y = "symbiont density")

# Assign object name to plot based on model name, and save
assign(paste("sim.plot.", model.name, sep = ""), sim.plot)
do.call(save, list(paste("sim.plot.", model.name, sep = ""), file = paste("Model_", model.version, "/sim.plot.", model.name, ".RData", sep = "")))

################################################################### Model fit ##################################################################

# Calculate sum-of-squares to estimate goodness-of-fit of simulations to empirical data
squares <- vector()
for (q in 1:Q){
  sim.q <- filter(sim.data, sim == q)
  # rio.female are the symbiont densities in females from Rio et al., 2006
  Obs <- filter(rio.female, time.step <= max(sim.q$time.step))
  Exp.density <- interp1(sim.q$time.step, sim.q$sym.pop, Obs$time.step)
  values <- (Exp.density-Obs$W.density)^2
  squares <- c(squares, values)
}
rss <- mean(squares)

# save sum-of-squares
save(rss, file = paste("Model_", model.version, "/rss.", model.name, ".RData", sep = ""))
print(proc.time() - ptm)
}
