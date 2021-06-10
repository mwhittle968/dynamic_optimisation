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
