#' Calculate the cross sectional area of the edge
#' 
#' This function adds the graph characteristic A which is the cross sectional area of the edge.
#'
#' @param g an igraph object. The graph representing the network
#' @param value a character string. The name of the edge attribute that is used as value from which Area will be calculated
#' @param minimum_value a numeric value. Indicating the most thinnest edge
#' @param range a numeric value. This gives the range of A values above the minimum.
#' 
#' @details This function is pretty niche but calculates a cross sectional area of an edge.
#' This is useful when you wish to calculate the spring coefficient k using Young's modulus. 
#' The function coerces and edge characteristic to be within a certain range of values preventing
#' negative/zero/infinite values.
#' 
#' @return a igraph object with the new edge attribute "Area" for each edge
#' 
#' @examples  
#' \dontrun{
#' library(igraph)
#' set.seed(234)
#' g_prep <- generate_peels_network("A") %>%
#'  set.edge.attribute(., name = "edge_characteristic", value = rep(1:16, each = 10))
#'
#' g <- calc_spring_area(g_prep, value = "edge_characteristic", minimum_value = 10, range = 20)
#'
#' get.edge.attribute(g, "Area")
#' }
#' @export

calc_spring_area <- function(g, value, minimum_value, range){

  temp <- igraph::as_data_frame(g) %>% tibble::as_tibble(.) %>%
    dplyr::rename(value_2 = !!value) %>%
    dplyr::mutate(
      value_2 = abs(value_2),
           A = dplyr::case_when(
             is.finite(value_2) ~  range*(value_2 - min(value_2))/(max(value_2)-min(value_2)) + minimum_value,
             TRUE ~range + minimum_value
           ),
            
           A = ifelse(is.finite(A), A, minimum_value)) #prevents NaNs from 0 values or other such annoying stuff
  
  g2 <- igraph::set.edge.attribute(g, "Area", value = temp$A)
  return(g2)
}
