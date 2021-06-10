#' Remove small components
#' 
#' keep only the largest component of graph
#' 
#' @param g An igraph object of the graph to embed.
#' 
#' @details As setse only works on connected components this function removes all but the largest component. 
#' This is a helper function to quickly project a network with setse.
#' 
#' @return An igraph object.
#' @examples 
#' library(igraph)
#' set.seed(1284)
#' #generate a random erdos renyi graph with 100 nodes and 150 edges
#' g <- erdos.renyi.game(n=100, p.or.m = 150, type = "gnm" )
#' #count the number of components
#' components(g)$no
#' 
#' #remove all but the largest component
#' g2 <-remove_small_components(g)
#' 
#' #Now there is only 1 component
#' igraph::components(g2)$no
#' 
#' @export
#' 

remove_small_components <- function(g){
  
  components_info <- igraph::components(g)
  max_comp <- which.max(components_info$csize)
  
  g <- igraph::delete.vertices(g, components_info$membership != max_comp)
  
  return(g)
  
}

