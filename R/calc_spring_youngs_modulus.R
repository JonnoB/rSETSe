#' Calculate the youngs modulus of the spring
#' 
#' This function adds the graph characteristic E which is the youngs modulus of each edge in the graph.
#'
#' @param g an igraph object. The graph representing the network
#' @param flow a character string. The name of the edge atribute that is used as flow
#' @param capacity a character string. The ame of the edge atribute for edge capacity
#' @param minimum_value a numeric value. Indicating the most stretchy value of youngs modulus
#' @param stretch_range a numeric value. This gives the range of k values above the minimum and the point when loading is 100% of capacity
#' 
#' @export

calc_spring_youngs_modulus <- function(g, flow, capacity, minimum_value, stretch_range){

  #This can be replaced by just multiplying vectors taken straight from edge attributes, it would be more
  #memory efficient especially on large graphs
  temp <- as_data_frame(g) %>% as_tibble %>%
    rename(flow_2 = !!flow,
           capacity_2 = !!capacity) %>%
    mutate(LL = abs(flow_2)/capacity_2,
           E = stretch_range*(1-LL) + minimum_value,
           E = ifelse(is.finite(E), E, minimum_value)) #if capacity is 0 NaNs are produced this prevents that
  
  g2 <- set.edge.attribute(g, "E", value = temp$E)
  
  return(g2)
}
