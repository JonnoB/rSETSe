Calc_Spring_Youngs_Modulus <- function(g, Force, Capacity, minimum_value, stretch_range){
  #g  an igraph object
  #Force a character string that is the name of the edge atribute that is used as force
  #Capacity a character string that is the ame of the edge atribute for edge capacity
  #minimum_value a numeric value indicating the most stretchy value of youngs modulos
  #stretch_range a numeric value giving the range of k values above the minimum
  #and the point when loading is 100% of capacity.
  
  temp <- as_data_frame(g) %>% as.tibble %>%
    rename(Force_2 = !!Force,
           Capacity_2 = !!Capacity) %>%
    mutate(LL = abs(Force_2)/Capacity_2,
           E = stretch_range*(1-LL) + minimum_value )
  
  g2 <- set.edge.attribute(g, "E", value = temp$E)
  return(g2)
}
