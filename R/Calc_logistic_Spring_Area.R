#' Logistic sprin area
#' 
#' Calculate the spring area of an edge based on a logistic curve
#' 
#' This function allows for variable cross sectional area of an edge. This might be because
#' there is a complex relationship between the flow and it's physical interpretation.
#' The output of this function is a graph with an area attribute "Area".
#' 
#' @param g  an igraph object
#' @param Force a character string that is the name of the edge atribute that is used as force
#' @param Capacity a character string that is the ame of the edge atribute for edge capacity
#' @param eta the range of possible logistic values
#' @param gamma the minimum value the logistic function will produce
#' @param kappa the squashing parameter of the logistic function
#' @export

Calc_logistic_Spring_Area <- function(g, Force, Capacity, eta, gamma, kappa){
  #g  an igraph object
  #Force a character string that is the name of the edge atribute that is used as force
  #Capacity a character string that is the ame of the edge atribute for edge capacity
  #eta the range of possible logistic values
  #gamma the minimum value the logistic function will produce
  #kappa the squashing parameter of the logistic function
  
  temp <- as_data_frame(g) %>% as.tibble %>%
    rename(Force_2 = !!Force,
           Capacity_2 = !!Capacity) %>%
    mutate(supports = ifelse(Force_2 > 0, from, to), #each edge "supports" a node. This node is the node that recieves flow.
           #The supporter node is the node that supplies flow. nodes can be both supporters and supported. However, a node cannot support a
           #node that it is supported by as the network is acyclic.
           Free_Capacity = Capacity_2-abs(Force_2)) %>%
    group_by(supports) %>%
    mutate(Total_Free_Capacity = sum(Free_Capacity), #find the free capacity of each supported node
           total_supports = n()) %>%
    ungroup %>%
    mutate(x =  abs(Force_2)*(1- 1/(Total_Free_Capacity - Free_Capacity)),
           x0 = Free_Capacity/2,
           A = abs(Force_2)*(eta/(1+exp(-kappa*(x-x0)))+ gamma))
  
  
  g2 <- set.edge.attribute(g, "Area", value = temp$A)
  
  return(g2)
  
}