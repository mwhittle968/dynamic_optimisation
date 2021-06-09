# runge-kutta method for numerically solving differential equations

rk4 <- function(population.change, sym.pop, a, n){
h <- 1/n # from t to t+1 (of the "big" model) there are n steps of h length
s <- array(data = NA, dim = n+1) # array to be populated with s values for t = 0 to n
t <- array(data = NA, dim = n+1) # array to be populated with time steps
s[1] = sym.pop # initial condition given by previous time step
t[1] = 0
for (j in 1:n){
  k.1 = population.change(t[j],s[j], a)
  k.2 = population.change(t[j]+(h/2), s[j]+((h*k.1)/2), a)
  k.3 = population.change(t[j]+(h/2), s[j]+((h*k.2)/2), a)
  k.4 = population.change(t[j]+h, s[j]+(h*k.3), a)
  s[j+1] = s[j] + (1/6)*h*(k.1+(2*k.2)+(2*k.3)+k.4)
  t[j+1] = t[j] + h
}
return(s[n+1])
}
