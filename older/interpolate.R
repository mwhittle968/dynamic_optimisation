# function for trilinear interpolation of values from a.opt array, using the new values of state variables as input
interpolate <- function(t, x, y, z){
  x.i = floor(x)
  x.ii = ceiling(x)
  y.i = floor(y)
  y.ii = ceiling(y)
  z.i = floor(z)
  z.ii = ceiling(z)
  if(x.i == x.ii && y.i == y.ii && z.i == z.ii){
    a = a.opt[t, x, y, z]
    } else if (x.i != x.ii && y.i == y.ii && z.i == z.ii){
      a = a.opt[t, x.i, y, z] + (x - x.i)*(a.opt[t, x.ii, y, z] - a.opt[t, x.i, y, z])/(x.ii - x.i)
    } else if (x.i == x.ii && y.i != y.ii && z.i == z.ii){
      a = a.opt[t, x, y.i, z] + (y - y.i)*(a.opt[t, x, y.ii, z] - a.opt[t, x, y.i, z])/(y.ii - y.i)
    } else if (x.i == x.ii && y.i == y.ii && z.i != z.ii){
      a = a.opt[t, x, y, z.i] + (z - z.i)*(a.opt[t, x, y, z.ii] - a.opt[t, x, y, z.i])/(z.ii - z.i)
    } else if(x.i != x.ii && y.i != y.ii && z.i == z.ii){
      a.i = ((x.ii - x)/(x.ii - x.i))*a.opt[t, x.i, y.i, z] + ((x - x.i)/(x.ii - x.i))*a.opt[t, x.ii, y.i, z]
      a.ii = ((x.ii - x)/(x.ii - x.i))*a.opt[t, x.i, y.ii, z] + ((x - x.i)/(x.ii - x.i))*a.opt[t, x.ii, y.ii, z]
      a = ((y.ii - y)/(y.ii - y.i))*a.i + ((y - y.i)/(y.ii - y.i))*a.ii
    } else if(x.i != x.ii && y.i == y.ii && z.i != z.ii){
      a.i = ((x.ii - x)/(x.ii - x.i))*a.opt[t, x.i, y, z.i] + ((x - x.i)/(x.ii - x.i))*a.opt[t, x.ii, y, z.i]
      a.ii = ((x.ii - x)/(x.ii - x.i))*a.opt[t, x.i, y, z.ii] + ((x - x.i)/(x.ii - x.i))*a.opt[t, x.ii, y, z.ii]
      a = ((z.ii - z)/(z.ii - z.i))*a.i + ((z - z.i)/(z.ii - z.i))*a.ii
    } else if(x.i == x.ii && y.i != y.ii && z.i != z.ii){
      a.i = ((y.ii - y)/(y.ii - y.i))*a.opt[t, x, y.i, z.i] + ((y - y.i)/(y.ii - y.i))*a.opt[t, x, y.ii, z.i]
      a.ii = ((y.ii - y)/(y.ii - y.i))*a.opt[t, x, y.i, z.ii] + ((y - y.i)/(y.ii - y.i))*a.opt[t, x, y.ii, z.ii]
      a = ((z.ii - z)/(z.ii - z.i))*a.i + ((z - z.i)/(z.ii - z.i))*a.ii
    } else {
      x.d = (x - x.i)/(x.ii - x.i)
      y.d = (y - y.i)/(y.ii - y.i)
      z.d = (z - z.i)/(z.ii - z.i)
      a.i.i = a.opt[t, x.i, y.i, z.i]*(1-x.d) + a.opt[t, x.ii, y.i, z.i]*x.d
      a.i.ii = a.opt[t, x.i, y.i, z.ii]*(1-x.d) + a.opt[t, x.ii, y.i, z.ii]*x.d
      a.ii.i = a.opt[t, x.i, y.ii, z.i]*(1-x.d) + a.opt[t, x.ii, y.ii, z.i]*x.d
      a.ii.ii = a.opt[t, x.i, y.ii, z.ii]*(1-x.d) + a.opt[t, x.ii, y.ii, z.ii]*x.d
      a.i = a.i.i*(1-y.d) + a.ii.i*y.d
      a.ii = a.i.ii*(1-y.d) + a.ii.ii*y.d
      a = a.i*(1-z.d) + a.ii*z.d
     }
  return(a)
}