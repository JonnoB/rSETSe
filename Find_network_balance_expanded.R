Find_network_balance_expanded <- function(g, force ="BalencedPower", flow = "PowerFlow", tstep = 0.5, distance = "Imp", 
                                 mass = 2000, maxIter =2000, frctmultiplier = 1, tol = 1e-10, verbose = TRUE,
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
