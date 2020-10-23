#' Calculate for from a given angle
#' 
#' This is an internal function
#' Calculates the vertical force for a given angle for a network or bi-connected sub graph of two nodes.
#' 
#' @param target_angle A numeric. The angle for which the force is needed
#' @param k A numeric. The spring constant
#' @param d A numeric. The original distance between nodes
#' 
#' @details This function is used by the two_node_solution function to calculate the correct vertical distance between two nodes.
#' It has very little functionality outside small networks or bi-connected components.
#' 
#' @export

ForceV_from_angle <- function(target_angle = 5*pi/12, k, d=temp$d){
  
  H = sqrt(d^2 * (1 + tan(target_angle)^2))
  ForceT = k*(H-d)
  ForceV = ForceT*sin(target_angle)
  #It is really slow to use a tibble here keeping everything thing base R is much faster
  # tibble(H = sqrt(d^2 * (1 + tan(target_angle)^2)),
  #        ForceT = k*(H-d), 
  #        ForceV = ForceT*sin(target_angle)) %>%
  #   pull(ForceV)
  
  return(ForceV)
  
}