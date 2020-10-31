#' Calculate the spring constant
#' 
#' This function adds the graph characteristic k which is the spring constant. When A and distance are both set to 1
#' \code{k=E} and the spring constant is equivalent to Young's modulus.
#' 
#' @param g an igraph object. The graph representing the network
#' @param youngs_mod a character string. The Young's modulus of the edge. The default is E
#' @param A a character string. The cross sectional area of the line. The default is A. see details on values of A
#' @param distance A character string. See details on values of distance
#' @return and edge attribute called k with value \code{EA/distance}
#' @seealso [calc_spring_youngs_modulus]
#' 
#' @details 
#' The values A and distance are edge attributes referring to the cross-sectional area of the edge and the horizontal distance of the edge,
#' in other words the distance between the two nodes at each end of the edge. These values can be set to anything the user wishes, they may be
#' constant or not. However, consider carefully setting the values to anything else other than 1. There needs to be a clear reasoning
#' or the results will be meaningless. 
#' 
#' For example setting the distance of an edge that represents an electrical cable to the distance
#' of the electrical cable will return very different results when compared to a constant of one. However, the physical distance between two points
#' does not necessarily have an impact on the loading of the line and so the results would not be interpretable. In contrast setting the distance
#' metric to be some function of the line resistance may have meaning and be appropriate. As a general rule distance and area should be set to 1.
#' 
#' @examples
#' \dontrun{
#' library(igraph)
#' set.seed(234)
#' g_prep <- generate_peels_network("A") %>%
#'  set.edge.attribute(., name = "edge_characteristic", value = rep(1:16, each = 10)) %>%
#'  #set some pretend Young's modulus value
#'  set.edge.attribute(., name = "E", value = rep(c(1e5, 5e5, 2e5, 3e5), each = 40)) %>%
#'  #calculate the spring area from another edge characteristic
#'  calc_spring_area(., value = "edge_characteristic", minimum_value = 10, range = 20) %>%
#'  prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                     force_var = "class", 
#'                     positive_value = "A")
#'
#' g <- calc_spring_constant(g_prep, youngs_mod = "E", A = "Area", distance = "distance")
#'
#' }
#' 
#' 
#' @export
#' 
calc_spring_constant <- function(g, youngs_mod = "E", A = "Area", distance = "distance"){
# 
  
 youngs_mod_vect <- igraph::get.edge.attribute(g, name = youngs_mod)*igraph::get.edge.attribute(g, name = A)/
   igraph::get.edge.attribute(g, name = distance)
   # temp <- igraph::as_data_frame(g) %>% tibble::as_tibble() %>%
   #   dplyr::mutate(k = rlang::.data[[E]]*rlang::.data[[A]]/rlang::.data[[distance]])
   
   print(temp)

  g2 <- igraph::set.edge.attribute(g, "k", value = youngs_mod_vect
                                   #temp$k
                                   )
  return(g2)
}
