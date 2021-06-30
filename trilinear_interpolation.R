# function for trilinear interpolation of values from array, using the new values of state variables as input
interpolate <- function(array, x, y, z, t){
  if(x < 1 | y < 1 | z < 1){
    d = 0
  } else {
  x.i = floor(x)
  x.ii = ceiling(x)
  y.i = floor(y)
  y.ii = ceiling(y)
  z.i = floor(z)
  z.ii = ceiling(z)
  if(x.i == x.ii && y.i == y.ii && z.i == z.ii){
    d = array[t, x, y, z]
    } else if (x.i != x.ii && y.i == y.ii && z.i == z.ii){
      d = array[t, x.i, y, z] + (x - x.i)*(array[t, x.ii, y, z] - array[t, x.i, y, z])/(x.ii - x.i)
    } else if (x.i == x.ii && y.i != y.ii && z.i == z.ii){
      d = array[t, x, y.i, z] + (y - y.i)*(array[t, x, y.ii, z] - array[t, x, y.i, z])/(y.ii - y.i)
    } else if (x.i == x.ii && y.i == y.ii && z.i != z.ii){
      d = array[t, x, y, z.i] + (z - z.i)*(array[t, x, y, z.ii] - array[t, x, y, z.i])/(z.ii - z.i)
    } else if(x.i != x.ii && y.i != y.ii && z.i == z.ii){
      d.i = ((x.ii - x)/(x.ii - x.i))*array[t, x.i, y.i, z] + ((x - x.i)/(x.ii - x.i))*array[t, x.ii, y.i, z]
      d.ii = ((x.ii - x)/(x.ii - x.i))*array[t, x.i, y.ii, z] + ((x - x.i)/(x.ii - x.i))*array[t, x.ii, y.ii, z]
      d = ((y.ii - y)/(y.ii - y.i))*d.i + ((y - y.i)/(y.ii - y.i))*d.ii
    } else if(x.i != x.ii && y.i == y.ii && z.i != z.ii){
      d.i = ((x.ii - x)/(x.ii - x.i))*array[t, x.i, y, z.i] + ((x - x.i)/(x.ii - x.i))*array[t, x.ii, y, z.i]
      d.ii = ((x.ii - x)/(x.ii - x.i))*array[t, x.i, y, z.ii] + ((x - x.i)/(x.ii - x.i))*array[t, x.ii, y, z.ii]
      d = ((z.ii - z)/(z.ii - z.i))*d.i + ((z - z.i)/(z.ii - z.i))*d.ii
    } else if(x.i == x.ii && y.i != y.ii && z.i != z.ii){
      d.i = ((y.ii - y)/(y.ii - y.i))*array[t, x, y.i, z.i] + ((y - y.i)/(y.ii - y.i))*array[t, x, y.ii, z.i]
      d.ii = ((y.ii - y)/(y.ii - y.i))*array[t, x, y.i, z.ii] + ((y - y.i)/(y.ii - y.i))*array[t, x, y.ii, z.ii]
      d = ((z.ii - z)/(z.ii - z.i))*d.i + ((z - z.i)/(z.ii - z.i))*d.ii
    } else {
      x.d = (x - x.i)/(x.ii - x.i)
      y.d = (y - y.i)/(y.ii - y.i)
      z.d = (z - z.i)/(z.ii - z.i)
      d.i.i = array[t, x.i, y.i, z.i]*(1-x.d) + array[t, x.ii, y.i, z.i]*x.d
      d.i.ii = array[t, x.i, y.i, z.ii]*(1-x.d) + array[t, x.ii, y.i, z.ii]*x.d
      d.ii.i = array[t, x.i, y.ii, z.i]*(1-x.d) + array[t, x.ii, y.ii, z.i]*x.d
      d.ii.ii = array[t, x.i, y.ii, z.ii]*(1-x.d) + array[t, x.ii, y.ii, z.ii]*x.d
      d.i = d.i.i*(1-y.d) + d.ii.i*y.d
      d.ii = d.i.ii*(1-y.d) + d.ii.ii*y.d
      d = d.i*(1-z.d) + d.ii*z.d
    }
  }
  return(d)
}