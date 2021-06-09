# Chop function limits the range of values the state variables can take:

chop <- function(x, x.min, x.max){
  x = ifelse(x >= x.max, x.max, ifelse(x <= x.min, x.min, x))
  return(x)
}

# Extrinsic survival and state-dependent survival function based on maternal reserves

ex.survival <- 1 - 1/36 # Approx 36 time steps so 1/36 chance of dying
in.survival <- function(nh){
  S = 1/(1+exp((2-nh)))
  return(S)
}
# Fitness function based on maternal reserves and time step:

fitness <- function(nh, ns, t){
  if (t%%4 != 0 | nh < nh.critical | ns < ns.critical){
      W = 0
      }
  else {# If it's a reproductive time step  and the reserves are sufficient the fitness is calculated
    nh. = (nh-mh)/2 # Mother gives half her reserves to larva
    ns. = (ns-ms)/2
    W.nh = nh. # linear by nh
    W.ns = 1/(1+exp(ns.critical-ns.)) # logistic by ns
    W = W.nh*W.ns
    }
  return(W)
}

# Maintenance, investment, regulation and production

maintenance <- function(s, i){
  m = s/i
  return(m)
}

investment <- function(d, s, j){
  s. = s + (d-1)/(j*s)
  return(s.)
}

regulation <- function(d, s, k){
  s. = s - (d-N-1)/k
  return(s.)
}

production <- function(s, l){
  p = s/l
  return(p)
}

# New state variables:

new.state <- function(t, d, s, ns, nh, i, j, k, l){
  if (t%%4 != 0 | nh < nh.critical | ns < ns.critical){ # If it's a non-reproductive time step, or the reserves are below critical levels,
    new.ns = chop((ns - ms + production(s, l)), 0, ns.max) # no resources are passed to larva
    if (d == 1){ # If d = 1, there is no investment or regulation
    new.nh = chop((nh - mh + N - maintenance(s, i)), 0, nh.max)
    new.s = s
    } else if (d > 1 & d <= (N+1)){ # If d is between 2 and N+1, there is investment
      new.nh = chop((nh - mh + N - maintenance(s, i) - (d-1)), 0, nh.max)
      new.s = chop((investment(d, s, j)), 0, s.max)
    } else { # If d is greater than N+1, there is regulation
        new.nh = chop((nh - mh + N - maintenance(s, i) - (d-N-1)), 0 , nh.max)
        new.s = chop((regulation(d, s, k)), 0, s.max)
    }
    } else { # If it's a reproducing time step and the reserves are sufficient, half the reserves are passed to larva
      new.ns = chop((((ns - ms)/2) + production(s, l)), 0, ns.max)
     if (d == 1){
       new.nh = chop(((nh - mh)/2 + N - maintenance(s, i)), 0, nh.max)
       new.s = s
     } else if (d > 1 & d <= (N+1)){
       new.nh = chop(((nh - mh)/2 + N - maintenance(s, i) - (d-1)), 0, nh.max)
       new.s = chop((investment(d, s, j)), 0, s.max)
     } else {
       new.nh = chop(((nh - mh)/2 + N - maintenance(s, i) - (d-N-1)), 0 , nh.max)
       new.s = chop((regulation(d, s, k)), 0, s.max)
     }
    }
  output = data.frame(new.nh, new.ns, new.s)
  return(output)
}
