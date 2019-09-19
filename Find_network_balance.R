Find_network_balance <- function(g, force ="BalencedPower", flow = "PowerFlow", capacity = "capacity", distance = "distance", tstep = 0.5, 
                                 mass = 2000, maxIter =2000, frctmultiplier = 1, tol = 1e-10, verbose = TRUE,
                                 TwoNodeSolution = TRUE){
  #needs an edge attribute "distance"
  #needs an edge attribute "Link" for the the edge name
  #converges faster if the network has been decomposed into blocks
  #TwoNodeSolution: Logical value if true blocks that are a node pair will be solved by Newton Raphson method for speed
  
  #helper function that prepares the data
  Prep <- Prepare_data_for_find_network_balance(g, force, flow, distance, mass)
  
  #do special case solution
  if(nrow(Prep$Link)==1 & TwoNodeSolution){
    
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
      tstep = tstep, 
      maxIter = maxIter, 
      frctmultiplier = frctmultiplier, 
      tol = tol, 
      verbose = verbose) 
    
  }
  
  
  return(Out)
  
}
