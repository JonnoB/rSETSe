velocity <- function(u, a, t0, t1){
  #calculates the velocity of an object given and initial velocity an acceleration and time difference
  #u initial velocity
  #a acceleration
  #t0 initial time
  #t1 next time period
  
  return(u + a*(t1-t0))
}
