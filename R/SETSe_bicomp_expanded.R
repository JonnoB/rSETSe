#' Bicomponent SETS embedding expanded
#' 
#' A wrapper function that takes a prepared graph and outputs a list of the embeddings and the aggregate dynamics throughout the simulations.
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param coef_drag A numeric. This sets the multiplier of friction. Only use if you want to be annoyed and confused
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param sparse Logical. Whether sparse matrices will be used. This becomes valubale for larger networks
#' 
#' @return A list containing 3 dataframes, the dataframe of the node embeddings, edge embeddings, and network dynamics
#'@export
SETSe_bicomp_expanded <- function(g, 
                          force = "force", 
                          distance = "distance", 
                          edge_name = "edge_name",
                          k = "k",
                          tstep,
                          coef_drag = 1,
                          tol,
                          max_iter = 20000,
                          mass = 1,
                          sparse = FALSE){
  
  #seperate out the network into blocks
  List_of_BiConComps <- create_balanced_blocks(g, 
                                               force = force)
  
  #find the largest component and use that as the origin block
  OriginBlock_number <-List_of_BiConComps %>% map_dbl(~vcount(.x)) %>% which.max()
  
  #print("Giant component found")
  
  #use the largest block to set the simulation parameters k and m.
  #k needs to be sufficiently stretched to allow enough topology variation. otherwise all that happens is a 
  #surface angled in the direct of net power flow. Which is interesting but not that interesting
  OriginBlock <- SETSe_expanded(g = List_of_BiConComps[[OriginBlock_number]],
                                               force =force,
                                               distance = distance,
                                               edge_name = edge_name,
                                               tstep = tstep,
                                               tol = tol,
                                               max_iter = max_iter,
                                               coef_drag = coef_drag,
                                               mass = mass,
                                               sparse = sparse)
  
  #print("Origin block complete, beggining remaining blocks")
  
  #Calculate the height embeddings using the Orgin block as a base
  height_embeddings <- Create_stabilised_blocks_expanded(g = g,
                                                OriginBlock = OriginBlock,
                                                OriginBlock_number = OriginBlock_number,
                                                force = force,
                                                distance = distance,
                                                edge_name = edge_name,
                                                tstep = tstep,
                                                coef_drag = coef_drag,
                                                tol = tol,
                                                max_iter = max_iter,
                                                mass = mass,
                                                sparse = sparse)

  
  return(height_embeddings)
}