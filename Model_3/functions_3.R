# Chop function limits the range of values the state variables can take:

chop <- function(x, x.min, x.max){
  x = ifelse(x >= x.max, x.max, ifelse(x <= x.min, x.min, x))
  return(x)
}

# Extrinsic survival and state-dependent survival function based on maternal reserves

S <- function(nh, nh.crit, t){
  ex.survival = female$ex.survival[t] # Approx 36 time steps so 1/36 chance of dying
  in.survival = ifelse(nh < nh.crit, 0, 1)
  S = ex.survival*in.survival
  return(S)
}

# Fitness function based on maternal reserves and event step:

B <- function(nh, ns, nh.repro, ns.repro, t, m){
  if (female$reproduction == "no" | nh < nh.repro | ns < ns.repro){
      B = 0
      } else {# If it's a reproductive time step  and the reserves are sufficient the fitness is calculated
    nh.larva = (nh-female$mh[t])*0.6 # Mother gives half her reserves to larva
    ns.larva = (ns-female$ms[t])*0.6
    B.nh = (m*nh.larva)-m*(nh.repro-female$mh[t])*0.6 # linear by nh
    B.ns = ifelse(ns < (ns.repro-female$ms[t])*0.6, 0, 1) # step function by ns
    B = B.nh*B.ns
    }
  return(B)
}

# Maintenance, investment, regulation and production

maintenance <- function(s, i){
  m = s*i
  return(m)
}

investment <- function(s, j, d){
  s. = s + (d-1)*j
  return(s.)
}

regulation <- function(s, N, k, d){
  s. = s - (d-N-1)*k
  return(s.)
}

production <- function(s, l){
  p = s*l
  return(p)
}

# New state variables:

new.state <- function(nh, ns, s, nh.max, ns.max, s.max, nh.repro, ns.repro, t, N, i, j, k, l, d){
  nh.larva = ifelse(female$reproduction == "no" | nh < nh.repro | ns < ns.repro, 0, (nh - female$mh[t])*0.6)
  ns.larva = ifelse(female$reproduction == "no" | nh < nh.repro | ns < ns.repro, 0, (ns - female$ms[t])*0.6)
  new.ns = chop((ns - female$ms[t] - ns.larva + production(s, l)), 0, ns.max)
    if (d == 1){ # If d = 1, there is no investment or regulation
      new.nh = chop((nh - nh.larva - female$mh[t] + female$n[t] - maintenance(s, i)), 0, nh.max)
      new.s = s
      } else if (d > 1 & d <= (N+1)){ # If d is between 2 and N+1, there is investment
      new.nh = chop((nh - female$mh[t] - nh.larva + female$n[t] - maintenance(s, i) - (d-1)), 0, nh.max)
      new.s = chop((investment(s, j, d)), 0, s.max)
      } else { # If d is greater than N+1, there is regulation
      new.nh = chop((nh - female$mh[t] - nh.larva + female$n[t] - maintenance(s, i) - (d-N-1)), 0 , nh.max)
      new.s = chop((regulation(s, N, k, d)), 0, s.max)
      }
  output = data.frame(new.nh, new.ns, new.s)
  return(output)
}
