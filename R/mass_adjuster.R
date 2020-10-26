#' Mass adjuster
#' 
#' This function adjusts the mass of the nodes so that the force in each direction over the mass for that direction 
#' produces an acceleration of 1.
#' 
#' @param g An igraph object. the network
#' @param force A character string. The name of the network attribute contain the network forces. Default is "force"
#' @param resolution_limit logical. If the forces in the network are smaller than the square root of the machine floating point limit
#' then the mass is set to one. default is true
#' 
#' @details This function can help stabilise the convergence of networks by preventing major imbalanes between the force in the network
#' and the mass of the nodes. in certain cases acceleration can become very large or very small in 
#' if force and mass are not well parametrised. 
#' 
#' This function means that if the network were reduced to two nodes where each node contained all the mass and all the force of
#' one of the two directions, then each node would have an acceleration of 1ms^-2
#' 
#' The function can become important when using setset_bicomp as the force mass ratio of biconnection components can vary widely from
#' the total force mass ratio of the network.
#' 
#' @return A numeric value giving the adjusted mass of the nodes in the network.
#' 
#' @examples 
#' set.seed(234) #set the random see for generating the network
#' 
#' g <- generate_peels_network(type = "E") %>%
#' prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                      force_var = "class", 
#'                      positive_value = "A")
#' 
#' mass_adjuster(g, force = "force", resolution_limit = TRUE)
#' 
#' @export
#' 
mass_adjuster <- function(g, force = "force", resolution_limit = TRUE){
  
  total_force <- sum(abs(igraph::vertex_attr(g, name = force)))
  
  if(resolution_limit &  total_force > .Machine$double.eps^0.5){
    
    mass <- 1
    
  } else{
    
    mass <- total_force /igraph::vcount(g)
    
  }
  
  return(mass)
  
}
