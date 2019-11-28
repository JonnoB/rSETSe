#' Calculate the spring constant
#' 
#' This function adds the graph characteristic k which is the spring constant. When A and distance are both set to 1
#' \code{k=E} and the spring constant is equivalent to youngs modulus.
#' 
#' The values A and distance are edge atributes referring to the crossectional area of the edge and the horizontal distance of the edge,
#' in other words the distance between the two nodes at each end of the edge. These values can be set to anything the user wishes, they may be
#' constant or not. However, consider carefully setting the values to anything else other than 1. There needs to be a clear resoning
#' or the results will be meaningless. 
#' 
#' For example setting the distance of an edge that represents an electrical cable to the distance
#' of the electrical cable will return very different results when compared to a constant of one. However, the physical distance between two points
#' does not necessarily have an impact on the loading of the line and so the results would not be interpretable. In contrast setting the distance
#' metric to be some function of the line resistance may have meaning and be appropriate. As a general rule distance and area should be set to 1.
#'
#' @param g an igraph object. The graph representing the network
#' @param E a character string. The youngs modulus of the edge. The default is E
#' @param A a character string. The cross sectional area of the line. The defualt is A. see details on values of A
#' @param distance A character string. See details on values of distance
#' @return and edge atribute called k with value \code{EA/distance}
#' @seealso [calc_spring_youngs_modulus]
#' @export
#' 
calc_spring_constant <- function(g, E = "E", A = "A", distance = "distance"){

  temp <- as_data_frame(g) %>% as.tibble %>%
    mutate(k = .data[[E]]*.data[[A]]/.data[[distance]])

  
  g2 <- set.edge.attribute(g, "k", value = temp$k)
  return(g2)
}
