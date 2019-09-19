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