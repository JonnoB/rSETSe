#' SETSe embedding on each bi-connected component using SETSe_auto
#' 
#' Embeds/smooths a feature network using the SETSe algorithm automatically finding convergence parameters using a grid search. In addition it breaks
#' the network into bi-connected component solves each sub-component inidividually and re-assembles them back into a single network. 
#' This is the most reliable method to perform SETSe embeddings and can be substantially quicker on certain network topologies.
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks. 
#' Default is set to NULL and call mass_adjuster to set the mass for each biconnected component
#' @param sparse Logical. Whether sparse matrices will be used. This becomes valuable for larger networks
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is the system absolute mean force.
#' @param hyper_iters integer. The hyper parameter that determines the number of iterations allowed to find an acceptable convergence value.
#' @param hyper_tol numeric. The convergence tolerance when trying to find the minimum value
#' @param hyper_max integer. The maximum number of iterations that SETSe will go through whilst searching for the minimum.
#' @param drag_min integer. A power of ten. The lowest drag value to be used in the search
#' @param drag_max integer. A power of ten. if the drag exceeds this value the tstep is reduced
#' @param tstep_change numeric. A value between 0 and 1 that determines how much the time step will be reduced by default value is 0.5
#' @param verbose Logical. This value sets whether messages generated during the process are suppressed or not.
#' @param noisy_termination Stop the process if the static force does not monotonically decrease.
#'
#'@details
#' Embedding the network by solving each bi-connected component then re-assembling can be faster for larger graphs, graphs with many nodes of degree 2, 
#' or networks with a low clustering coefficient.
#' This is because although SETSe is very efficient the topology of larger graphs make them more difficult to converge.
#' Large graph tend to be made of 1 very large biconnected component and many very small biconnected components. As the mass of the 
#' system is concentrated in the major biconnected component smaller ones can be knocked around by minor movements of the largest component. This
#' can lead to long convergence times. By solving all biconnected components separately and then reassembling the block tree at the end,
#' the system can be converged considerably faster. 
#' 
#' Setting mass to the absolute system force divided by the total nodes, often leads to faster convergence. As such
#' When mass is left to the default of NULL, the mean absolute force value is used.
#' 
#' @return A list containing 5 dataframes.
#' \enumerate{
#'   \item The node embeddings. Includes all data on the nodes the forces exerted on them position and dynamics at simulation termination
#'   \item The network dynamics describing several key figures of the network during the convergence process, this includes the static_force
#'   \item memory_df A dataframe recording the iteration history of the convergence of each component.
#'   \item Time taken. A data frame giving the time taken for the simulation as well as the number of nodes and edges. Node and edge data is given
#'   as this may differ from the total number of nodes and edges in the network depending on the method used for convergence.
#'   For example if SETSe_bicomp is used then some simulations may contain as little as two nodes and 1 edge
#'   \item The edge embeddings. Includes all data on the edges as well as the strain and tension values.
#' }
#' @family SETSe
# @seealso \code{\link{SETSe_auto}} \code{\link{SETSe}}
#' @examples
#' set.seed(234) #set the random see for generating the network
#' g <- generate_peels_network(type = "E")
#' embeddings <- g %>%
#' #prepare the network for a binary embedding
#' prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                      force_var = "class", 
#'                      positive_value = "A") %>%
#' #embed the network using  SETSe_bicomp
#'   SETSe_bicomp(tol = 0.02)
#' @export
SETSe_bicomp <- function(g, 
                         force = "force",
                         distance = "distance",
                         edge_name = "edge_name",
                         k = "k",
                         tstep = 0.02,
                         tol,
                         max_iter = 20000,
                         mass = NULL,
                         sparse = FALSE,
                         sample = 100,
                         static_limit = NULL,
                         hyper_iters = 100,
                         hyper_tol  = 0.1,
                         hyper_max = 30000,
                         drag_min = 0.01,
                         drag_max = 100,
                         tstep_change = 0.2,
                         verbose = FALSE,
                         noisy_termination = TRUE
){
  
  if(verbose){print("finding biconnected components")}
  
  start_time_bigraph <- Sys.time()
  bigraph <- igraph::biconnected_components(g)
  if(verbose){print(paste("Biconnected components found. Time taken",
                          round(as.numeric( difftime(Sys.time(), start_time_bigraph, units = "mins")), 1),
                          "minutes."))}
  
  #if the network cannot be decomposed into biconnected components then
  #create balanced blocks throughs an error and create_Stabilised blocks throughs an error
  #This prevents that
  if(bigraph$no==1){
    if(verbose){print("Network has no bi-connected components, running auto-SETSe instead")}
    embeddings_data <- SETSe_auto(g = g,
                                  force = force,
                                  distance = distance, 
                                  edge_name = edge_name,
                                  k = k,
                                  tstep = tstep, 
                                  tol = tol, #the force has to be scaled to the component 
                                  max_iter =  max_iter, 
                                  mass =  ifelse(is.null(mass), mass_adjuster(g, force = force, resolution_limit = TRUE), mass), 
                                  sparse = sparse,
                                  sample = sample,
                                  static_limit = static_limit,
                                  hyper_iters = hyper_iters,
                                  hyper_tol = hyper_tol,
                                  hyper_max = hyper_max,
                                  drag_min = drag_min,
                                  drag_max = drag_max,
                                  tstep_change = tstep_change,
                                  verbose = verbose,
                                  include_edges = FALSE,
                                  noisy_termination = noisy_termination)
    
  } else {
    
    #separate out the network into blocks
    if(verbose){print("creating balanced blocks")}
    start_time_bb <-Sys.time()
    balanced_blocks <- create_balanced_blocks(g, 
                                              force = force,
                                              bigraph = bigraph)
    
    if(verbose){print(paste("Balanced blocks created. Time taken",
                            round(as.numeric( difftime(Sys.time(), start_time_bb, units = "mins")), 1),
                            "minutes.",
                            "Total number of bi-connected components",
                            bigraph$no))}
    
    
    
    #find the largest component and use that as the origin block
    OriginBlock_number <-balanced_blocks %>% purrr::map_dbl(~igraph::vcount(.x)) %>% which.max()
    
    #total in network
    total_force <- sum(abs(igraph::get.vertex.attribute(g, force)))
    
    if(!is.null(static_limit)){
      #If the static limit is NULL the below returns a numeric vector of 0 length. By embedding the expression in an if statement the
      #error is prevented
      static_limit <- static_limit*sum(abs(igraph::get.vertex.attribute(balanced_blocks[[OriginBlock_number]], force)))/total_force
    }
    
    #calculate the parameters of the largest block
    if(verbose){print("Calculating Origin Block")}
    
    
    
    #do special case solution for two nodes only
    if(igraph::ecount(balanced_blocks[[OriginBlock_number]])==1){
      
      Prep <- SETSe_data_prep(g = balanced_blocks[[OriginBlock_number]], 
                              force = force, 
                              distance = distance, 
                              mass = ifelse(is.null(mass), mass_adjuster(balanced_blocks[[OriginBlock_number]], 
                                                                         force = force, resolution_limit = TRUE), mass), 
                              k = k,
                              edge_name = edge_name,
                              sparse = sparse)
      
      OriginBlock <- two_node_solution(g, Prep = Prep, auto_setse_mode = TRUE)
      
    } else {
      
      start_time_origin <- Sys.time()
      OriginBlock <- SETSe_auto(g = balanced_blocks[[OriginBlock_number]],
                                force = force,
                                distance = distance, 
                                edge_name = edge_name,
                                k = k,
                                tstep = tstep, 
                                tol = tol*sum(abs(igraph::get.vertex.attribute(balanced_blocks[[OriginBlock_number]], force)))/total_force, #the force has to be scaled to the component 
                                max_iter =  max_iter, 
                                mass =  ifelse(is.null(mass), mass_adjuster(balanced_blocks[[OriginBlock_number]], 
                                                                            force = force, resolution_limit = TRUE), mass), 
                                sparse = sparse,
                                sample = sample,
                                static_limit = static_limit,
                                hyper_iters = hyper_iters,
                                hyper_tol = hyper_tol,
                                hyper_max = hyper_max,
                                drag_min = drag_min,
                                drag_max = drag_max,
                                tstep_change = tstep_change,
                                verbose = verbose,
                                include_edges = FALSE,
                                noisy_termination = noisy_termination)
      
    }
    
    
    if(verbose){print(paste("Origin block complete, time taken", 
                            round(as.numeric( difftime(Sys.time(), start_time_origin, units = "mins")), 1), 
                            "minutes. beginning remaining blocks"))}
    
    
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
                                                drag_min = drag_min,
                                                drag_max = drag_max,
                                                tstep_change = tstep_change,
                                                sample = sample,
                                                static_limit = static_limit,
                                                verbose = verbose,
                                                bigraph = bigraph,
                                                balanced_blocks = balanced_blocks,
                                                noisy_termination = noisy_termination)
    
  } 
  # print("Height embeddings complete")
  
  #Extract edge tension and strain from the network
  embeddings_data$edge_embeddings <- calc_tension_strain(g = g,
                                                         embeddings_data$node_embeddings,
                                                         distance = distance,
                                                         edge_name = edge_name, 
                                                         k = k)

  #SETSe bicomp returns the nodes in a different order to the original graph. This can cause stress if you are not aware
  #This left join ensures the correct ordering. it is relatively but only happens once and the cost is small compared to the overall operation
  #The correct edge order is returned by calc_tension_strain
  embeddings_data$node_embeddings <- tibble::tibble(node = get.vertex.attribute(g, "name")) %>%
    dplyr::left_join(embeddings_data$node_embeddings)
  
  return(embeddings_data)
}
