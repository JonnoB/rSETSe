#' Distance
#' Calculates the distance a node will be from the origin for a given time-step
#' @param z A numeric value. The distance at current time
#' @param v A numeric value. Velcoity at current time
#' @param a A numeric value. Acelleration at current time
#' @param t0 A numeric value. Current time in seconds from simulation start
#' @param t1 A numeric value. Time at next time step aka t1 = t0 tdelta t
#' 
#' @export

Distance <- function(z, v, a, t0, t1 ){
  #Calculates the total distance a node will be from the origin for a given time step
  #z, distance at current time
  #v velcoity at current time
  #a acceleration
  #t0 current time
  #t1 time at and of time step
  
  return(v*(t1-t0) + 0.5*a*(t1-t0) + z)
  
}