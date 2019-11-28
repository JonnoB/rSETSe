#' Find the network balance
#' 
#' This simulates the dynamics of the network for a set number of iterations or until convergence which ever is the sooner.
#' 
#' This function is often used in conjunction with \code{Create_stabilised_blocks} and \code{create_balanced_blocks}
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param flow A character string. This is the edge attribute that is the power flow on the edges.
#' @param capacity A character string. This is the edge attribute that is the flow limit of the edges.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param maxIter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param frctmultiplier A numeric. This sets the multiplier of friction. Only use if you want to be annoyed and confused
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' @param two_node_solution Logical. The 
#' 
#' @return A data frame with the height embeddings of the network
#' @seealso \code{\link{Create_stabilised_blocks}} \code{\link{create_balanced_blocks}}
#' @export

Find_network_balance <- function(g, 
                                 force ="net_generation", 
                                 flow = "power_flow", 
                                 capacity = "capacity", 
                                 distance = "distance", 
                                 edge_name = "edge_name",
                                 tstep = 0.5, 
                                 mass = 2000, 
                                 maxIter =2000, 
                                 frctmultiplier = 1, 
                                 tol = 1e-10, 
                                 verbose = TRUE,
                                 two_node_solution = TRUE){
  #needs an edge attribute "distance"
  #needs an edge attribute "Link" for the the edge name
  #converges faster if the network has been decomposed into blocks
  #two_node_solution: Logical value if true blocks that are a node pair will be solved by Newton Raphson method for speed
  
  #helper function that prepares the data
  Prep <- Prepare_data_for_find_network_balance(g, force, flow, distance, mass, edge_name)
  
  #do special case solution
  if(nrow(Prep$Link)==1 & two_node_solution){
    
    if(Prep$NodeStatus$force[1]==0 &Prep$NodeStatus$force[2]==0){
      
      solution_angle <-0
      
    } else {
      #uses the non-linear optimiser from minpack.lm to find the solution to the two node special case, this is much faster
      solution_angle <- nlsLM(Force ~ ForceV_from_angle(target_angle, k = k, d = d), 
                              start = c(target_angle = pi/4), 
                              data = list(Force = abs(Prep$NodeStatus$force[1]), k = Prep$Link$k, d = Prep$Link$distance), 
                              upper = pi/2) %>% coefficients()      
      
    }
    
    temp <- Prep$NodeStatus %>%
      mutate(z = ifelse(force>0, 
                        tan(solution_angle)/2, #height above mid point
                        -tan(solution_angle)/2 ), #height below mid-point
             acceleration = 0,
             t = 0,
             Delta_acceleration = 0
      ) 
    
    Out <- list(results = temp, 
                NodeList = temp
    )
    #Solves using the iterative method.
  } else{
    
    Out <- FindStabilSystem(
      g = g,
      distance = distance,
      NodeStatus = Prep$NodeStatus, 
      EdgeNode = Prep$A, 
      flow = flow,
      kvect = Prep$Link$k, 
      dvect = Prep$Link$distance,
      capacity = capacity,
      edge_name = edge_name,
      tstep = tstep, 
      maxIter = maxIter, 
      frctmultiplier = frctmultiplier, 
      tol = tol, 
      verbose = verbose) 
    
  }
  
  
  return(Out)
  
}
