
############################################### Dynamic optimisation model with dynamic programming ###############################################

###################################################################################################################################################

dynamic.optimisation.4 <- function(model.name, nh.max = 20, ns.max = 5, s.max = 20, nh.repro = 8, ns.repro = 3, s.repro = 7, nh.crit.pupae = 4, nh.crit.adult = 4,
                                   nh.larva = 0.6, ns.larva = 0.6, mh.pupae = 1, mh.adult = 2, ms.pupae = 2, ms.adult = 2, ex.surv.pupae = 1, ex.surv.adult = (1 - 1/61),
                                   T = 61, n = 5, m = 1, i, j, k, l, Q = 100){
ptm <- proc.time()

####################################################### Data frames for events and parameters ######################################################

parameters <- c(model.name, nh.max, ns.max, s.max, nh.repro, ns.repro, s.repro, nh.crit.pupae, nh.crit.adult, nh.larva, ns.larva,
                mh.pupae, mh.adult, ms.pupae, ms.adult, ex.surv.pupae, ex.surv.adult, T, n, m, i, j, k, l, Q)
parameters <- t(as.data.frame(parameters))
colnames(parameters) <- c("model.name", "nh.max", "ns.max", "s.max", "nh.repro", "ns.repro", "s.repro", "nh.crit.pupae", "nh.crit.adult", "nh.larva", "ns.larva",
                     "mh.pupae", "mh.adult", "ms.pupae", "ms.adult", "ex.surv.pupae", "ex.surv.adult", "T", "n", "m", "i", "j", "k", "l", "Q")
#save(parameters, file = paste("Model_4/parameters.", model.name, ".RData", sep = ""))

time.step <- c(1:61)
F <- as.data.frame(time.step)%>%
  mutate(stage = ifelse(time.step <= 13, "pupae", "adult"))%>%
  mutate(repro = ifelse(time.step %in% seq(from = 21, to = 61, by = 4), "yes", "no"))%>%
  mutate(N = ifelse(stage == "pupae", 0, n))%>%
  mutate(mh = ifelse(stage == "pupae", mh.pupae, mh.adult), ms = ifelse(stage == "pupae", ms.pupae, ms.adult))%>%
  mutate(ex.surv = ifelse(stage == "pupae", ex.surv.pupae, ex.surv.adult))%>%
  mutate(nh.crit = ifelse(stage == "pupae", nh.crit.pupae, nh.crit.adult))

#################################################################### Functions ######################################################

# Chop function limits the range of values the state variables can take:
chop <- function(x, x.min, x.max){
  x = ifelse(x >= x.max, x.max, ifelse(x <= x.min, x.min, x))
  return(x)
}

# Extrinsic survival and state-dependent survival function based on maternal reserves
S <- function(nh, t){
  ex.survival = F$ex.surv[t]
  in.survival = ifelse(nh < F$nh.crit[t], 0, 1)
  S = ex.survival*in.survival
  return(S)
}

# Fitness function based on maternal reserves and event step:
B <- function(nh, ns, s, t){
  if (F$repro[t] == "no" | nh < nh.repro | ns < ns.repro | s < s.repro){
    B = 0
  } else {# If it's a reproductive time step  and the reserves are sufficient the fitness is calculated
    nh.l = (nh-F$mh[t])*nh.larva # Mother gives some of her reserves to larva
    ns.l = (ns-F$ms[t])*ns.larva
    B.nh = (m*nh.l)-m*(nh.repro-F$mh[t])*nh.larva # linear by nh
    B.ns = ifelse(ns < (ns.repro-F$ms[t])*ns.larva, 0, 1) # step function by ns
    B = B.nh*B.ns
  }
  return(B)
}

# Maintenance, investment, regulation and production
maintenance <- function(s){
  m = s*i
  return(m)
}

investment <- function(s, d){
  s. = s + (d-1)*j
  return(s.)
}

regulation <- function(s, d){
  s. = s - (d-n-1)*k
  return(s.)
}

production <- function(s){
  p = s*l
  return(p)
}

# New state variables:
new.state <- function(nh, ns, s, t, d){
  nh.l = ifelse(F$repro[t] == "no" | nh < nh.repro | ns < ns.repro, 0, (nh - F$mh[t])*nh.larva)
  ns.l = ifelse(F$repro[t] == "no" | nh < nh.repro | ns < ns.repro, 0, (ns - F$ms[t])*ns.larva)
  new.ns = chop((ns - F$ms[t] - ns.larva + production(s)), 0, ns.max)
  if (d == 1){ # If d = 1, there is no investment or regulation
    new.nh = chop((nh - nh.l - F$mh[t] + F$N[t] - maintenance(s)), 0, nh.max)
    new.s = s
  } else if (d > 1 & d <= (n+1)){ # If d is between 2 and N+1, there is investment
    new.nh = chop((nh - F$mh[t] - nh.l + F$N[t] - maintenance(s) - (d-1)), 0, nh.max)
    new.s = chop((investment(s, d)), 0, s.max)
  } else { # If d is greater than N+1, there is regulation
    new.nh = chop((nh - F$mh[t] - nh.l + F$N[t] - maintenance(s) - (d-n-1)), 0 , nh.max)
    new.s = chop((regulation(s, d)), 0, s.max)
  }
  output = data.frame(new.nh, new.ns, new.s)
  return(output)
}

######################################################## Dynamic Programming with for loops #######################################################

pb <- txtProgressBar(min = 0, max = ns.max*nh.max*s.max*(T-1), style = 3)
V <- array(data = NA, dim = c(T, nh.max, ns.max, s.max)) # Empty fitness array, to be populated with values from T back to 1
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
assign(paste("d.opt.", model.name, sep = ""), d.opt)
do.call(save, list(paste("d.opt.", model.name, sep = ""), file = paste("Model_4/d.opt.", model.name, ".RData", sep = "")))

############################################################# Optimisation array plot ############################################################

d.opt. <- melt(d.opt)
  names(d.opt.) <- c("time.step", "Nh", "Ns", "sym.pop", "decision")
  d.opt.plot <- ggplot(data = d.opt., aes(time.step, sym.pop, fill = decision)) +
    geom_tile()+
    facet_grid(Nh ~ Ns)+
    scale_fill_gradient(low = "#80CBC4", high = "#00695C")+
    theme(axis.title.y.right = element_text("Nh"),
          axis.title.x.top = element_text("Ns"))
assign(paste("d.opt.plot.", model.name, sep = ""), d.opt.plot)
do.call(save, list(paste("d.opt.plot.", model.name, sep = ""), file = paste("Model_4/d.opt.plot.", model.name, ".RData", sep = "")))

################################################################## Simulations ##################################################################

sim.data <- data.frame(sim = integer(), time.step = integer(), nh.reserves = integer(), ns.reserves = integer(), sym.pop = integer(), decision = integer())
for (q in 1:Q){
  sim.data <- add_row(sim.data, sim = q, time.step = 1, nh.reserves = chop(round(rnorm(1, mean = (nh.max), sd = 1)), 0, nh.max),
                      ns.reserves = chop(round(rnorm(1, mean = (ns.max), sd = 1)), 0, ns.max),
                      sym.pop = chop(round(rnorm(1, mean = (s.max/5), sd = 1)), 0, s.max), decision = NA)
  for (t in 1:(T-1)){
    sim.q <- filter(sim.data, sim == q)
    if (sim.q$nh.reserves[t] >= F$nh.crit[t]){
    d = interpolate(d.opt, chop(sim.q$nh.reserves[t], 1, nh.max), chop(sim.q$ns.reserves[t], 1, ns.max), chop(sim.q$sym.pop[t], 1, s.max), t)
    sim.data$decision[sim.data$sim == q & sim.data$time.step == t] <- d
    nh. = as.numeric(new.state(sim.q$nh.reserves[t], sim.q$ns.reserves[t], sim.q$sym.pop[t], t, d)[1])
    ns. = as.numeric(new.state(sim.q$nh.reserves[t], sim.q$ns.reserves[t], sim.q$sym.pop[t], t, d)[2])
    s. = as.numeric(new.state(sim.q$nh.reserves[t], sim.q$ns.reserves[t], sim.q$sym.pop[t], t, d)[3])
    sim.data <- add_row(sim.data, sim = q, time.step = (t+1), nh.reserves = nh., ns.reserves = ns., sym.pop = s., decision = NA)
    } else {
    break
  }
  }
}
assign(paste("sim.data.", model.name, sep = ""), sim.data)
do.call(save, list(paste("sim.data.", model.name, sep = ""), file = paste("Model_4/sim.data.", model.name, ".RData", sep = "")))

sim.plot <- ggplot(data = sim.data, aes(x = time.step, y = sym.pop, group = sim))+
  geom_line()+
  labs(x = "time step", y = "symbiont density")
assign(paste("sim.plot.", model.name, sep = ""), sim.plot)
do.call(save, list(paste("sim.plot.", model.name, sep = ""), file = paste("Model_4/sim.plot.", model.name, ".RData", sep = "")))

################################################################### Model fit ##################################################################

squares <- vector()
for (q in 1:Q){
  sim.q <- filter(sim.data, sim == q)
  Obs <- filter(Obs.F.data, time.step <= max(sim.q$time.step))
  Exp.density <- interp1(sim.q$time.step, sim.q$sym.pop, Obs$time.step)
  values <- (Exp.density-Obs$W.density)^2
  squares <- c(squares, values)
}
rss <- mean(squares)
save(rss, file = paste("Model_4/rss.", model.name, ".RData", sep = ""))
print(proc.time() - ptm)
}