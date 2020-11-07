#' SETSe algorithm with automatic timestep adjustment
#' 
#' The basic SETSe function with added timestep adjustment. 
#' The time shift functionality automatically adjusts the timestep if the convergence process is noisy
#'  
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param coef_drag A numeric. 
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param two_node_solution Logical. The Newton-Raphson algo is used to find the correct angle
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is twice the system absolute mean force.
#' @param tstep_change a numeric scaler. A value between 0 and one, the fraction the new timestep will be relative to the previous one
#' this can stop the momentum of the nodes forcing a divergence, but also can slow down the process. default is TRUE.
#' 
#' @details This is the basic SETS embeddings algorithm, it outputs all elements of the embeddings as well as convergence dynamics. It is a
#' wrapper around the core SETS algorithm which requires data preparation and only produces node embeddings and network dynamics. 
#' There is little reason to use this function as \code{\link{SETSe_auto}} and \code{\link{SETSe_bicomp}} 
#' are faster and easier to use.
#' 
#' @return A list of three elements. A data frame with the height embeddings of the network, a data frame of the edge embeddings
#' as well as the convergence dynamics dataframe for the network.
#' 
#' @examples
#' biconnected_network %>%
#' prepare_SETSe_continuous(., node_names = "name", force_var = "force") %>%
#' #embed the network using SETSe
#' SETSe_shift(., k = "weight", tstep = 0.000029)
#' @seealso \code{\link{SETSe_auto}} \code{\link{SETSe_bicomp}}
#' @export

SETSe_shift <- function(g, 
                  force ="force", 
                  distance = "distance", 
                  edge_name = "edge_name",
                  k ="k",
                  tstep = 0.02, 
                  mass = 1, 
                  max_iter = 20000, 
                  coef_drag = 1, 
                  tol = 1e-6,
                  sparse = FALSE,
                  two_node_solution = TRUE,
                  sample = 1,
                  static_limit = NULL,
                  tstep_change = 0.5){
  
  #helper function that prepares the data
  Prep <- SETSe_data_prep(g = g, 
                          force = force, 
                          distance = distance, 
                          mass = mass, 
                          k = k,
                          edge_name = edge_name,
                          sparse = sparse)
  
  #do special case solution 
  if(igraph::ecount(g)==1 & two_node_solution){
    
    Out <- two_node_solution(g, Prep = Prep, auto_setse_mode = FALSE)
    
    #Solves using the iterative method.
  } else{
    
    #The core algorithm
    Out <- SETSe_core_time_shift(
      node_embeddings = Prep$node_embeddings, 
      ten_mat = Prep$ten_mat, 
      non_empty_matrix = Prep$non_empty_matrix, 
      kvect = Prep$kvect, 
      dvect = Prep$dvect, 
      mass = mass,
      tstep = tstep, 
      max_iter = max_iter, 
      coef_drag = coef_drag,
      tol = tol, 
      sparse = sparse,
      sample = sample,
      static_limit = static_limit,
      tstep_change = tstep_change) 
    
  }
  
  
  #Extract edge tension and strain from the network
  Out$edge_embeddings <- calc_tension_strain(g = g,
                                             Out$node_embeddings,
                                             distance = distance, 
                                             edge_name = edge_name, 
                                             k = k)
  
  
  return(Out)
  
}
