#' Calculate line tension and strain
#' 
#' This function calculates the line  tension and strain characteristics for the edges in a graph
#' 
#' Whilst the height_embeddings_df contains the elevation of the SETS alogirthm this function produces a data frame that contains the Tension
#' and Strain. The dataframe that is returned contains a substantial amount of line information so reducing the number of variables may be
#' necessary if the data frame will be merged with previously generated data as there could be multiple columns of the same value.
#' 
#' 
#' @param g An igraph object of the network.
#' @param height_embeddings_df A data frame. This is the results of Create_stabilised_blocks or Find_network_balance
#' @param distance A character string. The name of the edge attribute that contains the distance between two nodes. The default is "distance"
#' @param capcity A character string. The name of the edge attribute that contains the flow capacity of the edge between two nodes. 
#' @param flow A character string. The name of the edge attribute that contains the flow between the two nodes at the end of the edge. The default is "power_flow".
#' @param edge_name   A character string. The name of the edge attribute that contains the edge name. The default is "edge_name".
#' @param k A character string. The name of the edge attribute that contains the sping coefficient
#' @return The function returns a data frame of 11 columns. These columns are the edge name, the alpha level, the line load (aka the inverse of alpha), 
#' the change in elevation, The final distance between the two nodes (the hypotenuse of the orignal distance and the vertical distance), 
#' the spring constant k, the edge tension, the edge strain, the percentile edge strain, the mean elevation and the total flow.
#' @export

calc_tension_strain <- function(g, height_embeddings_df, distance = "distance", capacity, flow = "power_flow", edge_name = "edge_name", k = "k"){

  #convert the character strings to symbols
  #afterwords the symbols can be evaluated by curly curly {{}}
  #Really replacing the sym function with inline .data[[]] would be better but I'll have to leave that for another time
  distance <- sym(distance)
  capacity <- sym(capacity)
  flow <- sym(flow)
  edge_name <- sym(edge_name)
  k <- sym(k)
  
  Out <- as_data_frame(g, what = "edges") %>% as_tibble %>%
    left_join(., height_embeddings_df %>% select(node, z), by = c("from"= "node")) %>%
    left_join(., height_embeddings_df %>% select(node, z), by = c("to"= "node")) %>%
    mutate(de = abs(z.x-z.y), #change in elevation
           mean_e = (z.x+z.y)/2, #mean elevation for the edge. simply the mean elevation for the nodes at both ends
           H = sqrt(de^2 +({{distance}})^2), #The total length of the edge. This is the hypotenuse of the triangle
           tension = {{k}}*H, #The tension in the edge following Hooks law
           strain = (H-{{distance}})/{{distance}}, #The mechanical strain of the edge
           alpha = {{capacity}}/abs({{flow}}), #The loading on the edge
           line_load = abs({{flow}})/{{capacity}}, #The line load of the edge aka inverse alpha
           percentile_strain = percent_rank(strain)) %>% #The strain rank relative to all edges in the network
    select({{edge_name}}, {{flow}}, {{capacity}}, alpha, line_load, de, mean_e, H, k, tension, strain, percentile_strain)
  
  return(Out)
  
}