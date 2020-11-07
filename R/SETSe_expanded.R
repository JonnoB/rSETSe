#' SETSe embedding showing full convergence history
#' 
#' This is a special case function which keeps the history of the network dynamics. It is useful for demonstrations. 
#' or parametrising difficult networks
#' 
#' @param g An igraph object. The network
#' @param force A character string
#' @param tstep A numeric. The time in seconds that elapses between each iteration
#' @param distance A character string. The name of the graph attribute that contains the graph distance
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param mass A numeric. The mass in kg of the nodes, this is arbitrary and commonly 1 is used. 
#' @param max_iter An integer. The maximum number of iterations before terminating the simulation
#' @param coef_drag A numeric. A multiplier used to tune the damping. Generally no need to twiddle
#' @param tol A numeric. Early termination. If the dynamics of the nodes fall below this value the algorithm will be classed as 
#' "converged" and the simulation terminates.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param verbose Logical value. Whether the function should output messages or run quietly.
#' @param two_node_solution Logical. The Newton-Raphson algo is used to find the correct angle
#' 
#' @return A list of four elements. A dat frame with the height embedding of the network, a data frame of the edge embeddings, 
#' the convergence dynamics dataframe for the network as well as the search history for convergence criteria of the network
#' 
#' @examples
#' 
#' g_prep <- biconnected_network %>%
#'  prepare_SETSe_continuous(., node_names = "name", force_var = "force", k = NULL)
#'
#' #the base configuration does not work
#' divergent_result <- SETSe_expanded(g_prep, k = "weight", tstep = 0.1)
#' 
#' #with a smaller timestep the algorithm converges
#' convergent_result <- SETSe_expanded(g_prep, k = "weight", tstep = 0.01)
#' 
#' \dontrun{
#' library(ggplot2)
#' #plot the results for a given node
#' convergent_result %>%
#'  ggplot(aes(x = t, y = net_force, colour = node)) + geom_line()
#' #replot with divergent_result to see what it looks like
#' }
#' @export

SETSe_expanded <- function(g, 
                           force ="force", 
                           distance = "distance", 
                           edge_name = "edge_name",
                           k = "k",
                           tstep = 0.02, 
                           mass = 1, 
                           max_iter = 20000, 
                           coef_drag = 1, 
                           tol = 1e-6,
                           sparse = FALSE,
                           verbose = TRUE,
                           two_node_solution = TRUE#,
                        #   include_edges = FALSE
){
  #needs an edge attribute "distance"
  #needs an edge attribute "Link" for the the edge name
  #converges faster if the network has been decomposed into blocks
  #TwoNodeSolution: Logical value if true blocks that are a node pair will be solved by Newton Raphson method for speed
  
  #
  #
  # This can be merged with the regular version when appropriate
  #
  #
  
  
  
  #helper function that prepares the data
  Prep <- SETSe_data_prep(g = g, 
                          force = force, 
                          distance = distance, 
                          mass = mass, 
                          edge_name = edge_name,
                          k = k,
                          sparse = sparse)
  
  #do special case solution I should change this to a standalone function for ease of reading but it isn't important
  if(igraph::ecount(g)==1 & two_node_solution){
    
    if(Prep$node_embeddings$force[1]==0 &Prep$node_embeddings$force[2]==0){
      
      solution_angle <-0
      
    } else {
      #uses the non-linear optimiser from minpack.lm to find the solution to the two node special case, this is much faster
      solution_angle <- minpack.lm::nlsLM(Force ~ ForceV_from_angle(target_angle, k = k, d = d), 
                              start = c(target_angle = pi/4), 
                              data = list(Force = abs(Prep$node_embeddings$force[1]), k = Prep$Link$k, d = Prep$Link$distance), 
                              upper = pi/2,
                              lower = 0) %>% stats::coefficients()      
      
    }
    
    Out <- Prep$node_embeddings %>%
      dplyr::mutate(elevation = ifelse(force>0, 
                                tan(solution_angle)/2, #height above mid point
                                -tan(solution_angle)/2 ), #height below mid-point
             net_force = 0,
             acceleration = 0,
             
             net_tension = ifelse(force>0, 
                                  -abs(Prep$node_embeddings$force[1]), #height above mid point
                                  abs(Prep$node_embeddings$force[1]))
      )  %>%
      dplyr::slice(rep(1:dplyr::n(), max_iter)) %>% #repeats the rows max_iter times so that
      dplyr::group_by(node) %>%
      dplyr::mutate(Iter = 1:max_iter,
             t = (tstep*Iter)) %>%
      dplyr::ungroup %>%
      dplyr::bind_rows(Prep$node_embeddings, .)
    
  } else{
    #Solves using the iterative method.   
    Out <- SETSe_core_expanded(
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
      sparse = sparse) 
    
  }
  
  
  # if(include_edges){
  #   
  #   #Extract edge tension and strain from the network
  #   Out$edge_embeddings <- calc_tension_strain(g = g,
  #                                              Out$node_embeddings,
  #                                              distance = distance, 
  #                                              edge_name = edge_name, 
  #                                              k = k)
  #   
  # }
  # 
  
  return(Out)
  
}
