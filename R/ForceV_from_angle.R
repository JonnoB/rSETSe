#' Calculate for from a given angle
#' 
#' This is mostly useful just for single edges. It may be removed from the package for being to specialised
#' 
#' @param target_angle A numeric. The angle for which the force is needed
#' @param k A numeric. The spring constant
#' @param d A numeric. The original distance between nodes
#' @export

ForceV_from_angle <- function(target_angle = 5*pi/12, k, d=temp$d){
  
  H = sqrt(d^2 * (1 + tan(target_angle)^2))
  ForceT = k*(H-d)
  ForceV = ForceT*sin(target_angle)
  #IT is reallt slow to use a tibble here keeping everything thing base R is much faster
  # tibble(H = sqrt(d^2 * (1 + tan(target_angle)^2)),
  #        ForceT = k*(H-d), 
  #        ForceV = ForceT*sin(target_angle)) %>%
  #   pull(ForceV)
  
  return(ForceV)
  
}