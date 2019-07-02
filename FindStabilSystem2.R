FindStabilSystem2 <- function(NodeStatus, EdgeNode, kvect, dvect,  tstep, maxIter = 1000, frctmultiplier = 1, 
                              tol = 1e-10, verbose = TRUE, window_size = Inf,  theta_var_tol = 13e-3, friction_Stop= FALSE, g ){
  #Runs the physics model to find the convergence of the system.
  
  #friction_stop fricton is a stopping condition. defualts to FALSE. 
  NodeList <- NodeStatus
  
  #results <- as.list(rep(NA,maxIter))
  results <- matrix(data = NA, nrow = maxIter, ncol = 11) %>%
    as_tibble() %>%
    set_names(c("t", "z", "NetForce", "velocity", "acceleration", "max_accel", 
                "max_Delta_accel", "friction", "theta_rads", "theta_degs",
                "theta_sd"))
  
  Iter <- 1
  system_stable <- FALSE
  
  while((Iter <= maxIter) & !system_stable ){
    
    # print(NodeList[[n]])
    temp <- NodeList %>%
      Calc_System_Dynamics(., EdgeNode, kvect, dvect, tstep, frctmultiplier)
    
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
      ungroup
    

    #Calculates current angle of the system
    current_angle <- angle_from_solved_heights(temp, g)

    results$theta_rads[Iter] <- current_angle$theta_rads
    results$theta_degs[Iter] <- current_angle$theta_degs
    # results[[Iter]] <- results[[Iter]] %>%
    #   mutate(theta_rads = current_angle$theta_rads,
    #        theta_degs = current_angle$theta_degs)
  
    NodeList <- temp
    
    #system_stable <- results[[Iter]]$max_accel < tol & results[[Iter]]$max_Delta_accel < tol #check if system is stable

    system_stable <- (results$max_accel[Iter] < tol & results$max_Delta_accel[Iter] < tol)#check if system is stable
    
    if(friction_Stop){
      
      #system_stable <- system_stable & results[[Iter]]$friction < tol #includes friction as a stopping condition. useful in some situations
      system_stable <- system_stable & results$friction[Iter] < tol #includes friction as a stopping condition. useful in some situations
    }
    
    
    if(verbose){
      
      # print(paste("Iteration", Iter, 
      #             "z", signif(results[[Iter]]$z, 3),
      #             "theta degs", signif(results[[Iter]]$theta_degs, 3),
      #             " velocity", signif(results[[Iter]]$velocity, 3), 
      #             "max acceleration", signif(results[[Iter]]$max_accel, 3), 
      #             "friction",  signif(results[[Iter]]$friction, 3) )
      #       ) # print result
      
      print(paste("Iteration", Iter, 
                  "z", signif(results$z[Iter], 3),
                  " velocity", signif(results$velocity[Iter], 3), 
                  "max acceleration", signif(results$max_accel[Iter], 3), 
                  "friction",  signif(results$friction[Iter], 3) )
      ) # print result
      
    }
    
    Iter <- Iter + 1 # add next iter

  }
  
  
  #Are all yhr accellerations below tolerance?
  #ProcessFinished<- sum(temp$acceleration<tol)==nrow(temp)
  
  
  #keep_elements <-results %>% map(~length(class(.x))==3) %>% unlist #removes empty list elements due to early terminations
  #empty elements are a logical NA value. otherwise they are a triple classed tibble. so can be filtered out.. tortourous and probably risky

  Out <- results %>%
    filter(!is.na(z)) %>%
    list(., NodeList)
  
  # Out <- results[keep_elements] %>% bind_rows() %>%
  #   list(., NodeList)
  
  names(Out) <- c("results", "NodeList")
  
  return(Out)

}
