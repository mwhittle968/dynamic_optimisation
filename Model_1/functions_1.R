# Chop function limits the range of values the state variables can take:

chop <- function(x, x.min, x.max){
  x = ifelse(x >= x.max, x.max, ifelse(x <= x.min, x.min, x))
  return(x)
}

# Extrinsic survival and state-dependent survival function based on maternal reserves

survival <- function(nh, nh.critical){
ex.survival <- 1 - 1/36 # Approx 36 time steps so 1/36 chance of dying
in.survival <- ifelse(nh < nh.critical, 0, 1)
S = ex.survival*in.survival
  return(S)
}

# Fitness function based on maternal reserves and time step:

fitness <- function(nh, ns, nh.reproduction, ns.reproduction, t, m){
  if (t%%4 != 0 | nh < nh.reproduction | ns < ns.reproduction){
      B = 0
      }
  else {# If it's a reproductive time step  and the reserves are sufficient the fitness is calculated
    nh.larva = (nh-mh)*0.6 # Mother gives half her reserves to larva
    ns.larva = (ns-ms)*0.6
    B.nh = (m*nh.larva)-m*(nh.reproduction-mh)*0.6 # linear by nh
    B.ns = ifelse(ns < (ns.reproduction-ms)*0.6, 0, 1) # step function by ns
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

regulation <- function(s, k, d){
  s. = s - (d-N-1)*k
  return(s.)
}

production <- function(s, l){
  p = s*l
  return(p)
}

# New state variables:

new.state <- function(nh, ns, s, nh.max, ns.max, s.max, nh.reproduction, ns.reproduction, mh, ms, t, N, i, j, k, l, d){
  nh.larva = (nh - mh)*0.6
  ns.larva = (ns - ms)*0.6
  if (t%%4 != 0 | nh < nh.reproduction | ns < ns.reproduction){ # If it's a non-reproductive time step, or the reserves are below critical levels,
    new.ns = chop((ns - ms + production(s, l)), 0, ns.max) # no resources are passed to larva
    if (d == 1){ # If d = 1, there is no investment or regulation
    new.nh = chop((nh - mh + N - maintenance(s, i)), 0, nh.max)
    new.s = s
    } else if (d > 1 & d <= (N+1)){ # If d is between 2 and N+1, there is investment
      new.nh = chop((nh - mh + N - maintenance(s, i) - (d-1)), 0, nh.max)
      new.s = chop((investment(s, j, d)), 0, s.max)
    } else { # If d is greater than N+1, there is regulation
        new.nh = chop((nh - mh + N - maintenance(s, i) - (d-N-1)), 0 , nh.max)
        new.s = chop((regulation(s, k, d)), 0, s.max)
    }
    } else { # If it's a reproducing time step and the reserves are sufficient, half the reserves are passed to larva
      new.ns = chop((ns - ns.larva) + production(s, l), 0, ns.max)
     if (d == 1){
       new.nh = chop(((nh - nh.larva) + N - maintenance(s, i)), 0, nh.max)
       new.s = s
     } else if (d > 1 & d <= (N+1)){
       new.nh = chop(((nh - nh.larva) + N - maintenance(s, i) - (d-1)), 0, nh.max)
       new.s = chop((investment(s, j, d)), 0, s.max)
     } else {
       new.nh = chop(((nh - nh.larva) + N - maintenance(s, i) - (d-N-1)), 0 , nh.max)
       new.s = chop((regulation(s, k, d)), 0, s.max)
     }
    }
  output = data.frame(new.nh, new.ns, new.s)
  return(output)
}
