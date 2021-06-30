# function for trilinear interpolation of values from array, using the new values of state variables as input
interpolate <- function(array, x, y, t){
  if(x < 1 | y < 1){
    d = 0
  } else {
  x.i = floor(x)
  x.ii = ceiling(x)
  y.i = floor(y)
  y.ii = ceiling(y)
  if(x.i == x.ii && y.i == y.ii){
    d = array[t, x, y]
    } else if (x.i != x.ii && y.i == y.ii){
      d = array[t, x.i, y] + (x - x.i)*(array[t, x.ii, y] - array[t, x.i, y])/(x.ii - x.i)
    } else if (x.i == x.ii && y.i != y.ii){
      d = array[t, x, y.i] + (y - y.i)*(array[t, x, y.ii] - array[t, x, y.i])/(y.ii - y.i)
    } else {
      d.i = ((x.ii - x)/(x.ii - x.i))*array[t, x.i, y.i] + ((x - x.i)/(x.ii - x.i))*array[t, x.ii, y.i]
      d.ii = ((x.ii - x)/(x.ii - x.i))*array[t, x.i, y.ii] + ((x - x.i)/(x.ii - x.i))*array[t, x.ii, y.ii]
      d = ((y.ii - y)/(y.ii - y.i))*d.i + ((y - y.i)/(y.ii - y.i))*d.ii
  return(d)
    }
  }
}