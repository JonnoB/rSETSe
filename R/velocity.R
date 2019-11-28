#' Node velocity
#' 
#' Calculate the velocity of the nodes at time step t. This is a helper function.
#' @param u A numeric. The velocity at time t0
#' @param a A numeric. The acceleration at time t0
#' @param t0 A numeric. The time at beginning of the time step
#' @param t1 A numeric. The time at the end of the time step.
#' @export
#' @examples 
#' velocity(0, 1, 2,5)

velocity <- function(u, a, t0, t1){
  #calculates the velocity of an object given and initial velocity an acceleration and time difference
  #u initial velocity
  #a acceleration
  #t0 initial time
  #t1 next time period
  
  return(u + a*(t1-t0))
}
