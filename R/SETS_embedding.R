#' calculate SETS embedding
#' 
#' A wrapper function that takes a prepared graph and outputs a list of the embeddings and the aggregate dynamics throughout the simulations.
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param flow A character string. This is the edge attribute that is the power flow on the edges.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param capacity A character string. This is the edge attribute that is the flow limit of the edges.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param coef_drag A numeric. This sets the multiplier of friction. Only use if you want to be annoyed and confused
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' @param sparse Logical. Whether sparse matrices will be used. This becomes valubale for larger networks
#' 
#' @return A list containing 3 dataframes, the dataframe of the node embeddings, edge embeddings, and network dynamics
#'@export
SETS_embedding <- function(g, 
                          force = "force", 
                          flow = "flow", 
                          distance = "distance", 
                          capacity = "capacity", 
                          edge_name = "edge_name",
                          k = "k",
                          tstep,
                          coef_drag = 1,
                          tol,
                          max_iter = 20000,
                          mass = 1,
                          #verbose = FALSE,
                          sparse = FALSE){
  
  #seperate out the network into blocks
  List_of_BiConComps <- create_balanced_blocks(g, 
                                               force = force, 
                                               flow = flow)
  
  #find the largest component and use that as the origin block
  giant_componant <-List_of_BiConComps %>% map_dbl(~vcount(.x)) %>% which.max()
  
  #print("Giant component found")
  
  #use the largest block to set the simulation parameters k and m.
  #k needs to be sufficiently stretched to allow enough topology variation. otherwise all that happens is a 
  #surface angled in the direct of net power flow. Which is interesting but not that interesting
  OriginBlock_complete <- Find_network_balance(g = List_of_BiConComps[[giant_componant]],
                                               force =force,
                                               flow = flow,
                                               distance = distance,
                                               edge_name = edge_name,
                                               tstep = tstep,
                                               tol = tol,
                                               max_iter = max_iter,
                                               coef_drag = coef_drag,
                                               mass = mass,
                                             #  verbose = verbose,
                                               sparse = sparse)
  
  #print("Origin block complete, beggining remaining blocks")
  
  #Calculate the height embeddings using the Orgin block as a base
  height_embeddings <- Create_stabilised_blocks(g = g,
                                                OriginBlock = OriginBlock_complete,
                                                OriginBlock_number = giant_componant,
                                                force = force,
                                                flow = flow,
                                                distance = distance,
                                                edge_name = edge_name,
                                                tstep = tstep,
                                                coef_drag = coef_drag,
                                                tol = tol,
                                                max_iter = max_iter,
                                                mass = mass,
                                               # verbose = verbose,
                                                sparse = sparse)
  
  # print("Height embeddings complete")
  
  #Extract edge tension and strain from the network
  tension_strain_embeddings <- calc_tension_strain(g = g,
                                                   height_embeddings$node_status,
                                                   distance = distance, 
                                                   capacity = capacity, 
                                                   flow = flow, 
                                                   edge_name = edge_name, 
                                                   k = k)
  
 # print("Strain and Tension embeddings complete")
  embeddings_data <- list(node_embeddings = height_embeddings$node_status, 
                          edge_embeddings = tension_strain_embeddings,
                          network_dynamics = height_embeddings$network_dynamics)
  
  return(embeddings_data)
}