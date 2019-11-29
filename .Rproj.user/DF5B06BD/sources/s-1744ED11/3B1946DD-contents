#' calculate SET
#' 
#' A wrapper function that takes a prepared graph and outputs a list of the embeddings
#' @param 
#'@export
calculate_SET <- function(g, 
                          force = "force", 
                          flow = "flow", 
                          distance = "distance", 
                          capacity = "capacity", 
                          edge_name = "edge_name",
                          k = "k",
                          tstep,
                          tol,
                          maxIter = 20000,
                          mass = 1,
                          verbose = FALSE){
  
  #seperate out the network into blocks
  List_of_BiConComps <- create_balanced_blocks(g, 
                                               force = force, 
                                               flow = flow)
  
  #find the largest component and use that as the origin block
  giant_componant <-List_of_BiConComps %>% map_dbl(~vcount(.x)) %>% which.max()
  
  print("Giant component found")
  
  #use the largest block to set the simulation parameters k and m.
  #k needs to be sufficiently stretched to allow enough topology variation. otherwise all that happens is a 
  #surface angled in the direct of net power flow. Which is interesting but not that interesting
  OriginBlock_complete <- Find_network_balance(g = List_of_BiConComps[[giant_componant]],
                                               force =force,
                                               flow = flow,
                                               distance = distance,
                                               capacity = capacity,
                                               edge_name = edge_name,
                                               tstep = tstep,
                                               tol = tol,
                                               maxIter = maxIter,
                                               mass = mass,
                                               verbose = verbose)
  
  print("Origin block complete")
  
  #Calculate the height embeddings using the Orgin block as a base
  height_embeddings_df <- Create_stabilised_blocks(g = g,
                                                   OriginBlock = OriginBlock_complete,
                                                   OriginBlock_number = giant_componant,
                                                   force = force,
                                                   flow = flow,
                                                   distance = distance,
                                                   capacity = capacity,
                                                   edge_name = edge_name,
                                                   tstep = tstep,
                                                   tol = tol,
                                                   maxIter = maxIter,
                                                   mass = mass,
                                                   verbose = verbose)
  
  print("Height embeddings complete")
  
  #Extract edge tension and strain from the network
  tension_strain_embeddings <- calc_tension_strain(g = g,
                                                   height_embeddings_df,
                                                   distance = distance, 
                                                   capacity = capacity, 
                                                   flow = flow, 
                                                   edge_name = edge_name, 
                                                   k = k)
  
  print("Strain and Tension embeddings complete")
  embeddings_data <- list(node_embeddings = height_embeddings_df, edge_embeddings = tension_strain_embeddings)
}