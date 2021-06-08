#' Calculate line tension and strain from the topology and node embeddings for high dimensional feature networks
#' 
#' This function calculates the line tension and strain characteristics for the edges in a graph.
#' It is called by default by all the embedding functions (SETSe_*) but is included here for completeness.
#' 
#' @param g An igraph object of the network.
#' @param height_embeddings_df A data frame. This is the results of Create_stabilised_blocks or Find_network_balance
#' @param distance A character string. The name of the edge attribute that contains the distance between two nodes. The default is "distance"
#' @param edge_name   A character string. The name of the edge attribute that contains the edge name. The default is "edge_name".
#' @param k A character string. The name of the edge attribute that contains the spring coefficient
#' 
#' @return The function returns a data frame of 7 columns. These columns are the edge name,
#' the change in elevation, The final distance between the two nodes (the hypotenuse of the original distance and the vertical distance), 
#' the spring constant k, the edge tension, the edge strain, and the mean elevation.
#
#' @details Whilst the node embeddings dataframe contains the elevation of the SETSe algorithm this function produces a data frame that contains the Tension
#' and Strain. The dataframe that is returned contains a substantial amount of line information so reducing the number of variables may be
#' necessary if the data frame will be merged with previously generated data as there could be multiple columns of the same value.
#' This function is called by default at the end of all SETSe functions
#' 
#' @examples 
#' 
#' set.seed(234) #set the random see for generating the network
#' g <- generate_peels_network(type = "E") %>%
#' #prepare the network for a binary embedding
#' prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                      force_var = "class", 
#'                      positive_value = "A") 
#' #embed the network using auto setse
#'  embeddings <- SETSe_auto(g)
#'   
#'  edge_embeddings_df <- calc_tension_strain(g, embeddings$node_embeddings)
#'  all.equal(embeddings$edge_embeddings, edge_embeddings_df) 
#'  
#'  
#' @export

calc_tension_strain_hd <- function(g, height_embeddings_df, distance = "distance", edge_name = "edge_name", k = "k"){
  
  #convert the character strings to symbols
  #afterwords the symbols can be evaluated by curly curly {{}}
  #Really replacing the sym function with inline .data[[]] would be better but I'll have to leave that for another time
  distance <- rlang::sym(distance)
  edge_name <- rlang::sym(edge_name)
  k <- rlang::sym(k)
  
  #get the edge list for the network
  temp <- igraph::as_data_frame(g, what = "edges") %>% 
    tibble::as_tibble(.) 
  
  #get the embedded node elevation across all dimensions
  elevation_df <- height_embeddings_df %>% dplyr::select(node, starts_with("elevation"))
  
  #merge the edge list tand the node elevations for both the from and to nodes
  #The columns are also re-named for clarity
  #This is a very slow way of doing this. Match would be much faster however I can't be bothered right now
  #and I will change if this becomes a speed choke point
  temp <- temp %>%
    dplyr::left_join(., 
                     elevation_df %>% dplyr::rename_with(., .fn = ~paste0("from_",.), .cols = -node)
                     , by = c("from"= "node")) %>%
    dplyr::left_join(., 
                     elevation_df %>% dplyr::rename_with(., .fn = ~paste0("to_",.), .cols = -node)
                     , by = c("to"= "node"))
  
  from_elevation_df <- temp %>% dplyr::select(dplyr::starts_with("from_elevation")) 
  to_elevation_df <-  temp  %>% dplyr::select(dplyr::starts_with("to_elevation")) 
 
  diff_df <- {from_elevation_df - to_elevation_df}
  names(diff_df) <- gsub("^from_", "elev_diff_", names(diff_df))
  
  
  Out <- temp  %>%
    dplyr::bind_cols(diff_df) %>% #change in elevation in each axis
    dplyr::mutate( 
      mean_e = sqrt(rowSums((diff_df/2)^2)), #mean elevation for the edge. simply the mean elevation for the nodes at both ends
      H = sqrt(rowSums(diff_df^2) +({{distance}})^2), #The total length of the edge. This is the hypotenuse of the triangle
      tension = {{k}}*(H-{{distance}}), #The tension in the edge following Hooks law
      strain = (H-{{distance}})/{{distance}} #The mechanical strain of the edge
    ) %>% 
    dplyr::select({{edge_name}}, dplyr::starts_with("elev_diff_"), mean_e, H, k, tension, strain)
  
  return(Out)
  
}
