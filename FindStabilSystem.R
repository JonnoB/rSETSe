FindStabilSystem <- function(NodeStatus, EdgeNode, kvect, dvect,  tstep, maxIter = 1000, frctmultiplier = 1){
  #This iterates through the graph for a set number of iterations hopefully until it converges
  #
  
  NodeList <- list(NodeStatus)
  
  for(n in 1:maxIter){
    
    
    # print(NodeList[[n]])
    temp <- NodeList[[n]] %>%
      Calc_System_Dynamics(., EdgeNode, kvect, dvect, tstep, frctmultiplier)
    
    NodeList[[n + 1]] <- temp
    
  }
  
  #Are all yhr accellerations below tolerance?
  #ProcessFinished<- sum(temp$acceleration<tol)==nrow(temp)
  
  return(NodeList)
  
}