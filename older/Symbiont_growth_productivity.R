# Dynamics of symbiont growth and nutrient production

# Packages
library(ggplot2)
library(cowplot)

# Variables:
# S = Symbiont density
# A = Allocated nutrients in host cell
# A. = Nutrients taken up by symbiont and available for use
# P = Nutrients produced by symbiont
# P. = Nutrients provisioned back to host and available for use

k1 <- 0.08 # Rate of acquisition of nutrients by the symbiont
k2 <- 0.1 # Rate of nutrient production by symbiont
k3 <- 0.005 # Rate of symbiont growth
k4 <- 10 # Rate of provisioning of nutrients back to host

a <- 0.7 # Proportion of nutrients used for conversion to other nutrients
b <- 0.3 # Proportion of nutrients used for symbiont growth

K <- 130 # Symbiont population carrying capacity

eqn.1 <- function(t, A, A., Ns., Ns, S){
  dAdt = -k1*A*S
  return(dAdt)
}  
eqn.2 <- function(t, A, A., Ns., Ns, S){
  dA.dt = (k1*A*S) - (a*k2*A.*S) - (b*k3*A.*S)#*(1-S/K)
  return(dA.dt)
}  
eqn.3 <- function(t, A, A., Ns., Ns, S){
  dNs.dt = (a*k2*A.*S) - (k4*Ns.)
  return(dNs.dt)
}
eqn.4 <- function(t, A, A., Ns., Ns, S){
  dNsdt = k4*Ns.
  return(dNsdt)
}
eqn.5 <- function(t, A, A., Ns., Ns, S){
  dSdt = (b*k3*A.*S)#*(1-S/K)
  return(dSdt)
}
fns <- c(eqn.1, eqn.2, eqn.3, eqn.4, eqn.5)

rk4 <- function(A, A., Ns., Ns, S, n){
  h <- 1/n
  vbs <- array(data = NA, dim = c((n+1), 5))
  t <- array(data = NA, dim = n+1)
  vbs[1,1] = A
  vbs[1,2] = A.
  vbs[1,3] = Ns.
  vbs[1,4] = Ns
  vbs[1,5] = S
  t[1] = 0
  for (j in 1:n){
    k <- array(data = NA, dim = c(4,5))
    for (i in 1:5){
      k[1,i] = fns[[i]](t[j], vbs[j,1], vbs[j,2], vbs[j,3], vbs[j,4], vbs[j,5])
      k[2,i] = fns[[i]](t[j]+(h/2), vbs[j,1]+((h*k[1,i])/2), vbs[j,2]+((h*k[1,i])/2), vbs[j,3]+((h*k[1,i])/2), vbs[j,4]+((h*k[1,i])/2), vbs[j,5]+((h*k[1,i])/2))
      k[3,i] = fns[[i]](t[j]+(h/2), vbs[j,1]+((h*k[2,i])/2), vbs[j,2]+((h*k[2,i])/2), vbs[j,3]+((h*k[2,i])/2), vbs[j,4]+((h*k[2,i])/2), vbs[j,5]+((h*k[2,i])/2))
      k[4,i] = fns[[i]](t[j]+h, vbs[j,1]+(h*k[3,i]), vbs[j,2]+(h*k[3,i]), vbs[j,3]+(h*k[3,i]), vbs[j,4]+(h*k[3,i]), vbs[j,5]+(h*k[3,i]))
      vbs[(j+1), i] = vbs[j, i] + (1/6)*h*(k[1,i]+(2*k[2,i])+(2*k[3,i])+k[4,i])
      t[j+1] = t[j] + h
    }
  }
  return(vbs)
}

output = rk4(1000, 0, 0, 0, 100, 100)
output <- as.data.frame(output)

A_plot<- ggplot(data = output)+
  geom_line(aes(x = c(1:101), y = V1))+
  xlab("time")+
  ylab("Nutrients available in bacteriocyte (A)")
A_plot 

A._plot<- ggplot(data = output)+
  geom_line(aes(x = c(1:101), y = V2))+
  xlab("time")+
  ylab("Nutrients taken up by symbiont (A.)")
A._plot

Ns._plot <- ggplot(data = output)+
  geom_line(aes(x = c(1:101), y = V3))+
  xlab("time")+
  ylab("Nutrients produced by symbiont (Ns.)")
Ns._plot

Ns_plot <- ggplot(data = output)+
  geom_line(aes(x = c(1:101), y = V4))+
  xlab("time")+
  ylab("Nutrients provisioned back to host (Ns)")
Ns_plot
  
S_plot<- ggplot(data = output)+
  geom_line(aes(x = c(1:101), y = V5))+
  xlab("time")+
  ylab("Symbiont density")
S_plot

plots <- plot_grid(A_plot, A._plot, Ns._plot, Ns_plot, S_plot)
plots
