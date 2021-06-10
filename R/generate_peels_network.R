#'Create a random Peel network
#'
#' Creates an example of a network from Peel's quintet of the specified type.
#' 
#' @param type A character which is any of the capital letters A-E
#' @param k_values An integer vector. The spring constant for the edge types within sub class, within class but not sub-class,
#' between classes. The default value is 1000, 500, 100. This means the strongest connection is for nodes in the same
#' sub-class and the weakest connection is for nodes in different classes
#' @param single_component Logical. Guarantees a single component network. Set to TRUE as default
#'  
#' @details  This function generates networks matching the 5 types described in Peel et al 2019(\url{doi.org/10.1073/pnas.1713019115}). All networks have 40 nodes,
#'  160 edges, two node classes and four node sub-classes. The connections between the are equal across all 5 types.
#'  As a result all networks generated have identical assortativity. However, as the sub-classes have different connection
#'  probability the structures produced by the networks are very different. When projected into setse space the network types
#'  occupy there own area, see Bourne 2020 for details
#'  
#' @return An igraph object that matches one of the 5 Peel's quintet types. The nodes are labeled with class and sub class.
#' The edges have attribute k which is the spring constant of the edge given relationship between the nodes the edge connects to
#' @examples
#' set.seed(234)
#' g <- generate_peels_network(type = "E")
#' plot(g)
#' @export
generate_peels_network <- function(type, k_values = c(1000, 500, 100), single_component = TRUE){

  if(single_component){
  
  #Common with other node embeddings methods setse only works on a single connected component
  #The below code ensures that the network created has a single component
  num_components <- 2
  
  while(num_components >1){
    g <- generate_peels_network_internal(type, k_values)
    num_components <- igraph::components(g)$no
  }
  
  } else {
    
    g <- generate_peels_network_internal(type, k_values)
  }
  return(g)
  
  
}
