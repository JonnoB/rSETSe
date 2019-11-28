#' Calculate the youngs modulus of the spring
#' 
#' This function adds the graph characteristic E which is the youngs modulus of each edge in the graph.
#'
#' @param g an igraph object. The graph representing the network
#' @param force a character string. The name of the edge atribute that is used as force
#' @param capacity a character string. The ame of the edge atribute for edge capacity
#' @param minimum_value a numeric value. Indicating the most stretchy value of youngs modulus
#' @param stretch_range a numeric value. This gives the range of k values above the minimum and the point when loading is 100% of capacity
#' @export
#' 
calc_spring_youngs_modulus <- function(g, force, capacity, minimum_value, stretch_range){
  #g  an igraph object
  #force a character string that is the name of the edge atribute that is used as force
  #capacity a character string that is the ame of the edge atribute for edge capacity
  #minimum_value a numeric value indicating the most stretchy value of youngs modulos
  #stretch_range a numeric value giving the range of k values above the minimum
  #and the point when loading is 100% of capacity.
  
  temp <- as_data_frame(g) %>% as.tibble %>%
    rename(force_2 = !!force,
           capacity_2 = !!capacity) %>%
    mutate(LL = abs(force_2)/capacity_2,
           E = stretch_range*(1-LL) + minimum_value)
  
  g2 <- set.edge.attribute(g, "E", value = temp$E)
  return(g2)
}
