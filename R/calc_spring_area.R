#' Calculate the youngs modulus of the spring
#' 
#' This function adds the graph characteristic A which is the cross sectional area of the springy.
#'
#' @param g an igraph object. The graph representing the network
#' @param value a character string. The name of the edge atribute that is used as value
#' @param minimum_value a numeric value. Indicating the most stretchy value of youngs modulus
#' @param range a numeric value. This gives the range of A values above the minimum, the maximum value value.
#' @export
#' 
calc_spring_area <- function(g, value, minimum_value, range){

  temp <- igraph::as_data_frame(g) %>% tibble::as_tibble %>%
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
