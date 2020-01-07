#' Find stabil system expanded
#' This is a special case function which keeps the history of the network dynamics. It is useful for demonstrations. or Parametrizing difficult networks
#' @param g An igraph object. The network
#' @param force A character string
#' @param flow A character string. the name of the graph attribute that contains the flow information
#' @param tstep A numeric. The time in seconds that elapses between each iteration
#' @param distance A character string. The name of the graph attribute that contains the graph distance
#' @param mass A numeric. The mass in kg of the nodes, this is arbitrary and commonly 1 is used. 
#' @param maxIter An interger. The maximum nuber of iterations before terminating the simulation
#' @param frctmultiplier A numeric. A multplier used to tune the damping. Generally no need to twiddle
#' @param tol A numeric. Early termination. If the dynamics of the nodes fall below this value the algortihm will be classed as 
#' "converged" and the simulation terminates.
#' @param verbose Logical value. Whether the function should output messages or run quietly.
#' @param two_node_solution Logical. The newton-raphson algo is used to find the correct angle
#' 
#' @export

Find_network_balance_expanded <- function(g, 
                                           force ="net_generation", 
                                           flow = "power_flow", 
                                           distance = "distance", 
                                           edge_name = "edge_name",
                                           tstep = 0.02, 
                                           mass = 1, 
                                           max_iter = 20000, 
                                           coef_drag = 1, 
                                           tol = 1e-6,
                                           sparse = FALSE,
                                           verbose = TRUE,
                                           two_node_solution = TRUE){
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
  Prep <- Prepare_data_for_find_network_balance(g = g, 
                                                force = force, 
                                                flow = flow, 
                                                distance = distance, 
                                                mass = mass, 
                                                edge_name = edge_name,
                                                sparse = sparse)
  
  #do special case solution I should change this to a standalone function for ease of reading but it isn't important
  if(nrow(Prep$Link)==1 & two_node_solution){
    
    if(Prep$node_status$force[1]==0 &Prep$node_status$force[2]==0){
      
      solution_angle <-0
      
    } else {
      #uses the non-linear optimiser from minpack.lm to find the solution to the two node special case, this is much faster
      solution_angle <- nlsLM(Force ~ ForceV_from_angle(target_angle, k = k, d = d), 
                              start = c(target_angle = pi/4), 
                              data = list(Force = abs(Prep$node_status$force[1]), k = Prep$Link$k, d = Prep$Link$distance), 
                              upper = pi/2) %>% coefficients()      
      
    }
    
    Out <- Prep$node_status %>%
      mutate(elevation = ifelse(force>0, 
                                tan(solution_angle)/2, #height above mid point
                                -tan(solution_angle)/2 ), #height below mid-point
             net_force = 0,
             acceleration = 0
      )  %>%
      slice(rep(1:n(), max_iter)) %>% #repeats the rows max_iter times so that
      group_by(node) %>%
    mutate(t = (tstep*(1:max_iter))) %>%
      ungroup

  } else{
    #Solves using the iterative method.   
    Out <- FindStabilSystem_expanded(
      node_status = Prep$node_status, 
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
      verbose = verbose) 
    
  }
  
  return(Out)
  
}
