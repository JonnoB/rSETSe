#' Function that finds the convergence parameters of a simulation
#' 
#' 

auto_SETS2 <- function(g, 
                      force ="net_generation", 
                      flow = "power_flow", 
                      distance = "distance", 
                      edge_name = "edge_name",
                      tstep = 0.02, 
                      mass = 1, 
                      max_iter = 20000, 
                      tol = 2e-3,
                      sparse = FALSE,
                      two_node_solution = TRUE,
                      restarts = 100,
                      step_size = 0.1,
                      sample = 1){
  
  
  memory_df<-tibble(iteration = 1:(2+restarts),
                    error = NA,
                    log_ratio = NA,
                    common_drag_iter = NA,
                    direction = 1,
                    target_area = NA,
                    res_stat = NA,
                    upper = NA,
                    lower = NA,
                    best_log_ratio =NA)
  #set initial round data
  memory_df$log_ratio[1] <- 0
  memory_df$common_drag_iter[1] <- 0
  memory_df$error[1] <-Inf#log10( 2-common_tol)
  memory_df$best_log_ratio[1] <-0
  memory_df$direction[1] = 1 #increasing
  memory_df$target_area[1] <- FALSE
  memory_df$res_stat[1] <- 2
  #These two are the limits of the log ratio. they are set to +/- infinity
  memory_df$upper[1] <- Inf
  memory_df$lower[1] <- -Inf
  memory_df$upper[2] <- Inf
  memory_df$lower[2] <- -Inf
  
  
  drag_iter<- 1
  
  memory_df$log_ratio[drag_iter+1] <- 1#-memory_df$log_ratio[drag_iter] - (memory_df$error[drag_iter])* 1 *direction_change
  memory_df$common_drag_iter[drag_iter+1] <- 10^( -memory_df$log_ratio[drag_iter+1]) * common_time
  
  while((drag_iter <= restarts) & (memory_df$res_stat[drag_iter]>common_tol)){
    drag_iter <- drag_iter+1      
    
    embeddings_data <- Find_network_balance(g, 
                                            force =force,
                                            flow = flow,
                                            distance = distance,
                                            edge_name = edge_name,
                                            tstep =  tstep,
                                            tol = tol,
                                            max_iter = max_iter,
                                            coef_drag =  memory_df$common_drag_iter[drag_iter],
                                            mass = mass,
                                            sample = sample,
                                            sparse = sparse,
                                            two_node_solution = two_node_solution
    )
    
    node_embeds <- embeddings_data$node_status
    
    memory_df$res_stat[drag_iter] <- sum(abs(node_embeds$static_force))
    
    #set the residual static force to 2 if it exceeds this value
    memory_df$res_stat[drag_iter] <- ifelse(memory_df$res_stat[drag_iter]>2, 2 ,memory_df$res_stat[drag_iter])
    
    #Is the algo in the convex area?
    memory_df$target_area[drag_iter] <- memory_df$res_stat[drag_iter] < 2
    
    #Log error is negative when small
    memory_df$error[drag_iter] <- log10(memory_df$res_stat[drag_iter] - common_tol) 
    
    #setting the limits
    if(min(memory_df$res_stat, na.rm = TRUE)<2){
      
      
      #ARE THESE EQUALITES THE RIGHT WAY ROUND?
      #THere is an issue when there are two unique values and one of them is NA
      min_error <- which.min(memory_df$error)
      upper_ratios <- memory_df$log_ratio[memory_df$log_ratio[min_error] >= memory_df$log_ratio] %>% unique()
      lower_ratios  <- memory_df$log_ratio[memory_df$log_ratio[min_error] <= memory_df$log_ratio] %>% unique()
      
      memory_df$upper[drag_iter+1] <-   ifelse(!is.finite( upper_ratios[rank(-upper_ratios, ties.method =  "first")==2]), 
                                               -Inf, 
                                               upper_ratios[rank(-upper_ratios, ties.method =  "first")==2])  
      # memory_df$upper[drag_iter+1] <-   ifelse(sum(is.finite(upper_ratios))==0, 
      #                                          -Inf, 
      #                                          upper_ratios[rank(-upper_ratios, ties.method =  "first")==2])
      
      memory_df$lower[drag_iter+1] <- ifelse(!is.finite(lower_ratios[rank(lower_ratios, ties.method =  "first" )==2]), 
                                             Inf, 
                                             lower_ratios[rank(lower_ratios, ties.method =  "first")==2])
      
      memory_df$best_log_ratio[drag_iter] <- memory_df$log_ratio[min_error]
    } else {
      memory_df$upper[drag_iter+1] <- memory_df$upper[drag_iter] 
      memory_df$lower[drag_iter+1] <-  memory_df$lower[drag_iter]
      memory_df$best_log_ratio[drag_iter] <- memory_df$best_log_ratio[drag_iter-1]
      
    }
    
    #if you are in the target zone updates happpen adaptively if not using a fixed step size
    if( is.finite(memory_df$upper[drag_iter]) ){
      # print("upper limit found")
      memory_df$log_ratio[drag_iter+1] <- (memory_df$upper[drag_iter] + memory_df$lower[drag_iter] )/2
      
      #if the logratios are the same break by choosing the lower limmit and the best value
      memory_df$log_ratio[drag_iter+1] <- ifelse(isTRUE(all.equal( memory_df$log_ratio[drag_iter+1],
                                                                   memory_df$log_ratio[drag_iter])), 
                                                 (memory_df$lower[drag_iter]+ memory_df$best_log_ratio[drag_iter])/2,
                                                 memory_df$log_ratio[drag_iter+1] )
      
    } else {
      #Before the non trivial upper and lower bound are found simple search is performed to find the target zone
      #fixed rate search across the plateau
      memory_df$log_ratio[drag_iter+1] <- memory_df$log_ratio[drag_iter] - step_size*memory_df$direction[drag_iter]
      
    }
    
    memory_df$common_drag_iter[drag_iter+1] <- 10^( -memory_df$log_ratio[drag_iter+1]) * common_time
    
    
    
    
    
    message_val <- ifelse(memory_df$direction[drag_iter] < 1, "accuracy increasing",  "accuracy decreasing")
    
    # print(paste("Iteration", drag_iter, message_val,  
    #             memory_df$common_drag_iter[drag_iter + 1], 
    #             "static force", memory_df$res_stat[drag_iter]))
    
  }
  
  
  
  memory_df <- memory_df %>%
    filter(complete.cases(.))
  
  embeddings_data$memory_df <- memory_df
  
  return(embeddings_data)
}