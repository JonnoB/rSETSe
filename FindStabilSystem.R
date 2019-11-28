#' Find stabil system
#' 
#' This function uses the SETS embedding to find the equilibrium position of a network or a bi-connected component
#' @param g An igraph object. The network
#' @param NodeStatus A data frame The current dynamics and forces experienced by the node a data frame.
#' @param EdgeNode 
#' @param flow A character string. This is the edge attribute that is the power flow on the edges.
#' @param kvect A numeric vector of the spring stiffnesses
#' @param dvect A numeric vector of the initial distances between the nodes
#' @param capacity A character string. This is the edge attribute that is the flow limit of the edges
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param tstep A numeric value. The time step, measured in seconds, that will be used to calculate the new dynamic state
#' @param maxIter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param frctmultiplier A numeric value. Used to set a multiplier on the friction value. Generally leave this alone..s.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' @param friction_stop Logical. Includes friction as a stopping condition. useful in some situations
#' @export


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
