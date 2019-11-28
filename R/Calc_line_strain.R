#' Calculate line strain
#' 
#' This function calculates the line strain characteristics for a graph
#' 
#' @param g An igraph object of the network
#' @param solved_heigh_df A data frame. This is the results of Create_stabilised_blocks or Find_network_balance
#' @param distance A character string. The name of the edge attribute that contains the distance between two nodes. The default is "distance"
#' @param capcity A character string. The name of the edge attribute that contains the flow capacity of the edge between two nodes. 
#' @param flow A character string. The name of the edge attribute that contains the flow between the two nodes at the end of the edge. The default is "power_flow".
#' @param edge_name   A character string. The name of the edge attribute that contains the edge name. The default is "edge_name".
#' @export

Calc_line_strain <- function(g, solved_height_df, distance = "distance", capacity, flow = "power_flow", edge_name = "edge_name"){

  #convert the character strings to symbols
  #afterwords the symbols can be evaluated by curly curly {{}}
  #Really replacing the sym function with inline .data[[]] would be better but I'll have to leave that for another time
  distance <- sym(distance)
  capacity <- sym(capacity)
  flow <- sym(flow)
  edge_name <- sym(edge_name)
  
  Out <- as_data_frame(g) %>% as_tibble %>%
    left_join(., solved_height_df %>% select(node, z), by = c("from"= "node")) %>%
    left_join(., solved_height_df %>% select(node, z), by = c("to"= "node")) %>%
    mutate(dz = abs(z.x-z.y),
           mean_z = (z.x+z.y)/2,
           H = sqrt(dz^2 +({{distance}})^2),
           strain = (H-{{distance}})/{{distance}},
           alpha = {{capacity}}/abs({{flow}}),
           line_load = abs({{flow}})/{{capacity}},
           percentile_strain = percent_rank(strain)) %>%
    select({{edge_name}}, alpha, line_load, dz, H, strain, percentile_strain, mean_z, {{flow}})
  
  return(Out)
  
}