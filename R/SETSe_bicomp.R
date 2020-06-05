#' Bi-component SETSe embedding
#' 
#' Performs the SETS embedding on a network using the bi-connected component/block-tree method. This is the most reliable method to 
#' perform SETSe embeddings and can be substantially quicker on unsual networks.
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param sparse Logical. Whether sparse matrices will be used. This becomes valubale for larger networks
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is the system absolute mean force.
#' @param hyper_iters integer. The hyper parameter that determines the number of iterations allowed to find an acceptable convergence value.
#' @param hyper_tol numeric. The convergence tolerance when trying to find the minimum value
#' @param hyper_max integer. The maximum number of iterations that the setse will go through whilst searching for the minimum.
#' @param step_size numeric. The hyper parameter that determines the log ratio search step size for auto convergence
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#'
#'@details
#' This approach can be faster for larger graphs or graphs with many nodes of degree 2, or networks with a low clustering coefficient.
#' This is becuase although the algorithm is very efficient the topology of larger graphs make them more difficult to converge.
#' Large graph tend to be made of 1 very large biconnected component and many very small biconnected components. As the mass of the 
#' system is concentrated in the major biconnected component the small ones can be knocked around by minor movements of the major. This
#' leads to long convergence times. By solving all biconnected components seperately and then resassmbling the block tree at the end,
#' the system can be converged considerably faster. In addition the smaller biconnected components iterate faster than a single large one.
#' 
#' Although the default mass is set to 1, setting mass to the absolute system force divided by the total nodes, often leads to faster convergence.
#' 
#' @return A list containing 5 dataframes.
#' \enumerate{
#'   \item The node embeddings. Includes all data on the nodes the forces exerted on them position and dynamics at simulation termination
#'   \item The network dynamics describing several key figures of the network during the convergence process, this includes the static_force
#'   \item memory_df A datframe recording ethe iteration history of the convergence of each component.
#'   \item A data frame giving the time taken for the simulation as well as the number of nodes and edges. Node and edge data is given
#'   as this may differ from the total number of nodes and edges in the network depending on the method used for convergnence.
#'   For example if SETSe_bicomp is used then some simulations may contain as little as two nodes and 1 edge
#'   \item time taken. the amount of time taken per component, includes the edge and nodes of each component
#' }
#' @export
SETSe_bicomp <- function(g, 
                          force = "force",
                          distance = "distance",
                          edge_name = "edge_name",
                          k = "k",
                          tstep = 0.02,
                          tol,
                          max_iter = 20000,
                          mass = 1,
                          sparse = FALSE,
                          sample = 100,
                          static_limit = NULL,
                          hyper_iters = 100,
                          hyper_tol  = 0.01,
                          hyper_max = 30000,
                          step_size = 0.1,
                          verbose = FALSE){
  if(verbose){print("finding biconnected components")}
  
  bigraph <- biconnected_components(g)
  
  #seperate out the network into blocks
  if(verbose){print("creating balanced blocks")}
  balanced_blocks <- create_balanced_blocks(g, 
                                            force = force,
                                            bigraph = bigraph)
  
  #find the largest component and use that as the origin block
  OriginBlock_number <-balanced_blocks %>% map_dbl(~vcount(.x)) %>% which.max()
  
  #print("Giant component found")
  
  #total in network
  total_force <- sum(abs(get.vertex.attribute(g, force)))
  
  if(!is.null(static_limit)){
    #this if statement prevents an error if the static limit is null the below returns a numeric vector of 0 length
    static_limit <- static_limit*sum(abs(get.vertex.attribute(balanced_blocks[[OriginBlock_number]], force)))/total_force
  }
  
  #calculate the parameters of the largest block
  if(verbose){print("Calculating Origin Block")}
  OriginBlock <- auto_SETSe(g = balanced_blocks[[OriginBlock_number]],
                            force = force,
                            distance = distance, 
                            edge_name = edge_name,
                            k = k,
                            tstep = tstep, 
                            tol = tol*sum(abs(get.vertex.attribute(balanced_blocks[[OriginBlock_number]], force)))/total_force, #the force has to be scaled to the component 
                            max_iter =  max_iter, 
                            mass =  mass, 
                            sparse = sparse,
                            sample = sample,
                            static_limit = static_limit,
                            hyper_iters = hyper_iters,
                            hyper_tol = hyper_tol,
                            hyper_max = hyper_max,
                            step_size = step_size,
                            verbose = verbose,
                            include_edges = FALSE )
  
  if(verbose){print("Origin block complete, beginning remaining blocks")}
  
  
  #Calculate the height embeddings using the Orgin block as a base
  embeddings_data <- Create_stabilised_blocks(g = g,
                                              OriginBlock = OriginBlock,
                                              OriginBlock_number = OriginBlock_number,
                                              force = force,
                                              k = "k",
                                              distance = distance,
                                              edge_name = edge_name,
                                              tstep = tstep,
                                              tol = tol,
                                              max_iter = max_iter,
                                              mass = mass,
                                              sparse = sparse,
                                              hyper_iters = hyper_iters,
                                              hyper_tol = hyper_tol,
                                              hyper_max = hyper_max,
                                              step_size = step_size,
                                              sample = sample,
                                              static_limit = static_limit,
                                              verbose = verbose,
                                              bigraph = bigraph,
                                              balanced_blocks = balanced_blocks)
  
  # print("Height embeddings complete")
  
  #Extract edge tension and strain from the network
  embeddings_data$edge_embeddings <- calc_tension_strain(g = g,
                                                         embeddings_data$node_embeddings,
                                                         distance = distance,
                                                         edge_name = edge_name, 
                                                         k = k)
  
  
  return(embeddings_data)
}