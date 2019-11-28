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
#' @param TwoNodeSolution A logical Value. Are two node networks solved using Newton-Raphson
#' 
#' @export
 
Find_network_balance_expanded <- function(g, force ="net_generation", flow = "power_flow", tstep = 0.5, distance = "Imp", 
                                 mass = 1, maxIter = 20000, frctmultiplier = 1, tol = 1e-10, verbose = TRUE,
                                 TwoNodeSolution = TRUE){
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
  Prep <- Prepare_data_for_find_network_balance(g, force, flow, distance, mass)
  
    
    Out <- FindStabilSystem_expanded(
      NodeStatus = Prep$NodeStatus, 
      EdgeNode = Prep$A, 
      kvect = Prep$Link$k, 
      dvect = pull(Prep$Link,distance), 
      tstep = tstep, 
      maxIter = maxIter, 
      frctmultiplier = frctmultiplier) 
    
  
  
  return(Out)
  
}
