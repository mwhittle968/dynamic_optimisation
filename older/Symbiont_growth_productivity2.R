# Dynamics of symbiont growth and nutrient production

# Packages
library(ggplot2)
library(cowplot)

# Variables:
# S = Symbiont density in bacteriocyte
# A = Allocated nutrient concentration in bacteriocyte
# P = Symbiont-produced nutrient concentration in bacteriocyte

symbiont_growth_productivity <- function(t, a, s, ns){

k1 = 0.05 # Rate of symbiont growth
k2 = 0.1 # Rate of nutrient production by symbiont
K = 10 # Symbiont population carrying capacity

eqn.1 <- function(t, a, s, ns){
  dadt = -(k1*a*s*(1-s/K))-(k2*a*s)
  return(dadt)
}  
eqn.2 <- function(t, a, s, ns){
  dsdt = k1*a*s*(1-s/K)
  return(dsdt)
}  
eqn.3 <- function(t, a, s, ns){
  dnsdt = k2*a*s
  return(dnsdt)
}

fns <- c(eqn.1, eqn.2, eqn.3)

rk4 <- function(a, s, ns){
  h <- 1/n
  vbs <- array(data = NA, dim = c((n+1), 3))
  t <- array(data = NA, dim = n+1)
  vbs[1,1] = a
  vbs[1,2] = s
  vbs[1,3] = ns
  t[1] = 0
  for (j in 1:n){
    k <- array(data = NA, dim = c(4,3))
    for (i in 1:3){
      k[1,i] = fns[[i]](t[j], vbs[j,1], vbs[j,2], vbs[j,3])
      k[2,i] = fns[[i]](t[j]+(h/2), vbs[j,1]+((h*k[1,i])/2), vbs[j,2]+((h*k[1,i])/2), vbs[j,3]+((h*k[1,i])/2))
      k[3,i] = fns[[i]](t[j]+(h/2), vbs[j,1]+((h*k[2,i])/2), vbs[j,2]+((h*k[2,i])/2), vbs[j,3]+((h*k[2,i])/2))
      k[4,i] = fns[[i]](t[j]+h, vbs[j,1]+(h*k[3,i]), vbs[j,2]+(h*k[3,i]), vbs[j,3]+(h*k[3,i]))
      vbs[(j+1), i] = vbs[j, i] + (1/6)*h*(k[1,i]+(2*k[2,i])+(2*k[3,i])+k[4,i])
      t[j+1] = t[j] + h
    }
  }
  return(vbs)
}
return(rk4(a,s,ns))
}

output = symbiont_growth_productivity(10, 4, 6, 3)
output <- as.data.frame(output)

A_plot<- ggplot(data = output)+
  geom_line(aes(x = c(1:21), y = V1))+
  xlab("time")+
  ylab("Nutrients available in bacteriocyte (A)")
A_plot 

Ns_plot <- ggplot(data = output)+
  geom_line(aes(x = c(1:21), y = V3))+
  xlab("time")+
  ylab("Nutrients provisioned back to host (Ns)")
Ns_plot

S_plot<- ggplot(data = output)+
  geom_line(aes(x = c(1:21), y = V2))+
  xlab("time")+
  ylab("Symbiont density")
S_plot

plots <- plot_grid(A_plot, Ns_plot, S_plot, align = "h")
plots
