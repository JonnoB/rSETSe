FindStabilSystem <- function(g, distance, NodeStatus, EdgeNode, flow, kvect, dvect, capacity, edge_name = edge_name, 
                             tstep, maxIter = 1000, frctmultiplier = 1, 
                              tol = 1e-10, verbose = TRUE, friction_Stop= FALSE ){
  #Runs the physics model to find the convergence of the system.
  
  #friction_stop fricton is a stopping condition. defualts to FALSE. 
  NodeList <- NodeStatus
  
  #results <- as.list(rep(NA,maxIter))
  results <- matrix(data = NA, nrow = maxIter, ncol = 9) %>%
    as_tibble() %>%
    set_names(c("t", "z", "NetForce", "velocity", "acceleration", "max_accel", 
                "max_Delta_accel", "friction", "strain"))
  
  Iter <- 1
  system_stable <- FALSE
  
  while((Iter <= maxIter) & !system_stable ){
    
    # print(NodeList[[n]])
    temp <- NodeList %>%
      Calc_System_Dynamics(., EdgeNode, kvect, dvect, tstep, frctmultiplier)
    
    #calculates the line strain each round
    line_strain <- Calc_line_strain(g, 
                                    solved_height_df = temp, 
                                    distance = distance, 
                                    capacity = capacity, 
                                    flow = flow,
                                    edge_name = edge_name
                                    )
   
    #results[[Iter]] 
    results[Iter,]<- temp %>%
      group_by(t) %>%
      summarise(z = mean(abs(z)),
                NetForce = mean(abs(NetForce)),
                velocity = mean(abs(velocity)),
                acceleration = mean(abs(acceleration)),
                max_accel = max(abs(acceleration)),
                max_Delta_accel = max(abs(Delta_acceleration)),
                friction = max(abs(friction))) %>%
      ungroup %>%
      mutate( strain =mean(line_strain$strain))

    NodeList <- temp

    system_stable <- (results$max_accel[Iter] < tol & results$max_Delta_accel[Iter] < tol)#check if system is stable
    
    if(friction_Stop){

      system_stable <- system_stable & results$friction[Iter] < tol #includes friction as a stopping condition. useful in some situations
    }
    
    
    if(verbose){

      print(paste("Iteration", Iter, 
                "strain", signif(results$strain[Iter], 3),
                  "z", signif(results$z[Iter], 3),
                  " velocity", signif(results$velocity[Iter], 3), 
                  "max acceleration", signif(results$max_accel[Iter], 3), 
                  "friction",  signif(results$friction[Iter], 3) )
      ) # print result
      
    }
    
    Iter <- Iter + 1 # add next iter

  }
  
  #empty elements are a logical NA value. otherwise they are a triple classed tibble. so can be filtered out.. tortourous and probably risky

  Out <- results %>%
    filter(!is.na(z)) %>%
    list(., NodeList)

  names(Out) <- c("results", "NodeList")
  
  return(Out)

}
