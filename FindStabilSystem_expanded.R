#' Find stabil system expanded
#' 
#' A helper function to Find_network_balance_expanded()
#' @param NodeStatus A data frame The current dynamics and forces experienced by the node a data frame.
#' @param EdgeNode 
#' @param kvect A numeric vector of the spring stiffnesses
#' @param dvect A numeric vector of the initial distances between the nodes
#' @param tstep A numeric value. The time step, measured in seconds, that will be used to calculate the new dynamic state
#' @param maxIter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param frctmultiplier A numeric value. Used to set a multiplier on the friction value. Generally leave this alone..s.
#' 
#' @seealso [Find_network_balance_expanded()]
#' @export

FindStabilSystem_expanded <- function(NodeStatus, EdgeNode, kvect, dvect,  tstep, maxIter = 1000, frctmultiplier = 1){
  #This iterates through the graph for a set number of iterations hopefully until it converges
  #It is the expanded version of Find Stabil system, the difference with this version is that it keeps all outputs of the Calc_System_Dynamics function.
  #The result is a very long data frame of convergence information. 
  #This can be useful for styduying convergence behaviour or making animations, but can be cumbersome for large graphs over many iterations.
  
  
  NodeList <- list(NodeStatus)
  
  for(n in 1:maxIter){
    
    
    # print(NodeList[[n]])
    temp <- NodeList[[n]] %>%
      Calc_System_Dynamics(., EdgeNode, kvect, dvect, tstep, frctmultiplier)
    
    NodeList[[n + 1]] <- temp
    
  }
  

  Out <- NodeList %>% bind_rows()
  
  return(Out)
  
}