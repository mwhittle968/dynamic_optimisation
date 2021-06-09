source("C:/Users/apico/Documents/Mathilda/PhD/Notes/R_code/symbiont_growth_productivity.R")
# Chop function limits the range of values the state variables can take:

chop <- function(x, x.min, x.max){
  x = ifelse(x >= x.max, x.max, ifelse(x <= x.min, x.min, x))
  return(x)
}

# Extrinsic survival and state-dependent survival function based on maternal reserves

ex.survival <- 1 - 1/36
in.survival <- function(nh){
  S = 1/(1+0.35*exp(0.74*(9-nh)))
}
# Fitness function based on maternal reserves and time step:

fitness <- function(nh, ns, t){
  if (t%%4 != 0 | nh <= nh.critical | ns <= ns.critical){
      W = 0
      }
  else {# If it's a reproductive time step  and the reserves are sufficient the fitness is calulated
    nh. = (nh-mh)/2 # Mother gives half her reserves to larva
    ns. = (ns-ms)/2
    W.nh = nh. # linear by nh
    W.ns = 1/(1+exp(ns.critical-ns.)) # logistic by ns
    W = W.nh*W.ns
    }
  return(W)
}

# New state variables:

new.state <- function(t, a, s, ns, nh){
  if (t%%4 != 0 | nh <= nh.critical | ns <= ns.critical){
    ns. = ns - ms
    new.nh = chop((nh - mh + N - a + symbiont_growth_productivity(a, s, ns.)[1]), 0, nh.max)
  }
   else {
     ns. = (ns - ms)/2 # mother gives half her reserves to larva
    new.nh = chop(((nh - mh)/2 + N - a + symbiont_growth_productivity(a, s, ns.)[1]), 0, nh.max)   # mother gives half her reserves to larva
   }
  new.ns = chop((symbiont_growth_productivity(a, s, ns.)[3]), 0, ns.max)
  new.s = chop((symbiont_growth_productivity(a, s, ns.)[2]), 0, s.max)
  output = data.frame(new.nh, new.ns, new.s)
  return(output)
}
