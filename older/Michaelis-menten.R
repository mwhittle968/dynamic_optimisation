# For an enzyme reaction E + S <-> C -> E + P, with rate constants k.f (forwards), k.b (backwards), k.cat (product formation)

k.f <- 1
k.b <- 0.3
k.cat <- 1.3

enzyme <- function(t, E, S, C, P){
  dEdt = -k.f*E*S + k.b*C + k.cat*C
  return(dEdt)
}  
substrate <- function(t, E, S, C, P){
  dSdt = -k.f*E*S + k.b*C
  return(dSdt)
}  
complex <- function(t, E, S, C, P){
  dCdt = k.f*E*S - k.b*C - k.cat*C
  return(dCdt)
}
product <- function(t, E, S, C, P){
  dPdt = k.cat*C
  return(dPdt)
}
fns <- c(enzyme, substrate, complex, product)

rk4.enzyme <- function(E.1, S.1, C.1, P.1, n){
  h <- 1/n
  vbs <- array(data = NA, dim = c((n+1), 4))
  t <- array(data = NA, dim = n+1)
  vbs[1,1] = E.1 # E[1] is the starting conc. of enzyme/symbiont
  vbs[1,2] = S.1 # S[1] is the starting conc. of substrate/allocation
  vbs[1,3] = C.1
  vbs[1,4] = P.1 # P[1] is the starting conc. of product/provisioned nutrients
  t[1] = 0
  for (j in 1:n){
    k <- array(data = NA, dim = c(4,4))
    for (i in 1:4){
    k[1,i] = fns[[i]](t[j], vbs[j,1], vbs[j,2], vbs[j,3], vbs[j,4])
    k[2,i] = fns[[i]](t[j]+(h/2), vbs[j,1]+((h*k[1,i])/2), vbs[j,2]+((h*k[1,i])/2), vbs[j,3]+((h*k[1,i])/2), vbs[j,4]+((h*k[1,i])/2))
    k[3,i] = fns[[i]](t[j]+(h/2), vbs[j,1]+((h*k[2,i])/2), vbs[j,2]+((h*k[2,i])/2), vbs[j,3]+((h*k[2,i])/2), vbs[j,4]+((h*k[2,i])/2))
    k[4,i] = fns[[i]](t[j]+h, vbs[j,1]+(h*k[3,i]), vbs[j,2]+(h*k[3,i]), vbs[j,3]+(h*k[3,i]), vbs[j,4]+(h*k[3,i]))
    vbs[(j+1), i] = vbs[j, i] + (1/6)*h*(k[1,i]+(2*k[2,i])+(2*k[3,i])+k[4,i])
    t[j+1] = t[j] + h
    }
  }
  return(vbs[,4])
}

output = rk4.enzyme(10, 10, 0, 0, 5)
output <- as.data.frame(output)
ggplot(data = output, aes(x = c(1:6), y = output))+
         geom_line()

       