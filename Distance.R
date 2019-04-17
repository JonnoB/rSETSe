Distance <- function(z, v, a, t0, t1 ){
  #Calculates the total distance a node will be from the origin for a given time step
  #z, distance at current time
  #v velcoity at current time
  #a acceleration
  #t0 current time
  #t1 time at and of time step
  
  return(v*(t1-t0) + 0.5*a*(t1-t0) + z)
  
}