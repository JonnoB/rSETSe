FindStabilSystem2 <- function(NodeStatus, EdgeNode, kvect, dvect,  tstep, maxIter = 1000, frctmultiplier = 1, 
                              tol = 1e-10, verbose = TRUE, friction_Stop= FALSE){
  #Runs the physics model to find the convergence of the system.
  
  #friction_stop fricton is a stopping condition. defualts to FALSE. 
  NodeList <- NodeStatus
  
  results <- as.list(rep(NA,maxIter))
  
  Iter <- 1
  system_stable <- FALSE
  
  while((Iter <= maxIter) & !system_stable ){
    
    
    # print(NodeList[[n]])
    temp <- NodeList %>%
      Calc_System_Dynamics(., EdgeNode, kvect, dvect, tstep, frctmultiplier)
    
    #Calculates current angle of the system
    
    results[[Iter]] <- temp %>%
      group_by(t) %>%
      summarise(z = mean(abs(z)),
                NetForce = mean(abs(NetForce)),
                velocity = mean(abs(velocity)),
                acceleration = mean(abs(acceleration)),
                max_accel = max(abs(acceleration)),
                max_Delta_accel = max(abs(Delta_acceleration)),
                friction = max(abs(friction))) %>%
      ungroup
    
  
    NodeList <- temp
    
    system_stable <- results[[Iter]]$max_accel < tol & results[[Iter]]$max_Delta_accel < tol #check if system is stable
    
    if(friction_Stop){
      
      system_stable <- system_stable & results[[Iter]]$friction < tol #includes friction as a stopping condition. useful in some situations
    }
    
    
    if(verbose){
      
      print(paste("Iteration", Iter, "z", signif(results[[Iter]]$z, 3),
                  " velocity", signif(results[[Iter]]$velocity, 3), 
                  "max acceleration", signif(results[[Iter]]$max_accel, 3), 
                  "friction",  signif(results[[Iter]]$friction, 3) )
            ) # print result
      
    }
    Iter <- Iter + 1 # add next iter

  }
  
  
  #Are all yhr accellerations below tolerance?
  #ProcessFinished<- sum(temp$acceleration<tol)==nrow(temp)
  
  
  keep_elements <-results %>% map(~length(class(.x))==3) %>% unlist #removes empty list elements due to early terminations
  #empty elements are a logical NA value. otherwise they are a triple classed tibble. so can be filtered out.. tortourous and probably risky

  
  Out <- results[keep_elements] %>% bind_rows() %>%
    list(., NodeList)
  
  names(Out) <- c("results", "NodeList")
  
  return(Out)

}
