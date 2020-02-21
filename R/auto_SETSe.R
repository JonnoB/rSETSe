#' Auto SETSe
#'
#' Uses a grid search and a binary search to find appropriate convergence conditions.
#' 
#' This function is pretty useful, it takes advantage of the linear relationship between the timestep and the coefficient of drag
#' to search along the log linear line formed by tstep/coef_drag to find the convergence conditions.
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param flow A character string. This is the edge attribute that is the power flow on the edges.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param capacity A character string. This is the edge attribute that is the flow limit of the edges.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param hyper_iters integer. The hyper parameter that determines the number of iterations allowed to find an acceptable convergence value.
#' @param step_size numeric. The hyper parameter that determines the log ratio search step size for auto convergence
#' @param hyper_tol numeric. The convergence tolerance when trying to find the minimum value
#' @param hyper_max integer. The maximum number of iterations that the setse will go through whilst searching for the minimum.
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#'  
#' @return A list of four elements. A data frame with the height embeddings of the network, a data frame of the edge embeddings, 
#' the convergence dynamics dataframe for the network as well as the search history for convergence criteria of the network
#' @export

auto_SETSe <- function(g, 
                        force ="net_generation", 
                        flow = "power_flow", 
                        distance = "distance", 
                        capacity = "capacity", 
                        edge_name = "edge_name",
                        k = "k",
                        tstep = 0.02, 
                        mass = 1, 
                        max_iter = 100000, 
                        tol = 2e-3,
                        sparse = FALSE,
                        hyper_iters = 100,
                        hyper_tol  = 0.01,
                        hyper_max = 30000,
                        step_size = 0.1,
                        sample = 1,
                        verbose = FALSE){
  
  
  memory_df<-tibble(iteration = 1:(2+hyper_iters),
                    error = NA,
                    perc_change = NA,
                    log_ratio = NA,
                    common_drag_iter = NA,
                    direction = 1,
                    target_area = NA,
                    res_stat = NA,
                    upper = NA,
                    lower = NA,
                    best_log_ratio =NA,
                    stable = NA)
  #set initial round data
  memory_df$log_ratio[1] <- 0
  memory_df$common_drag_iter[1] <- 0
  memory_df$error[1] <-Inf
  memory_df$best_log_ratio[1] <-0
  memory_df$direction[1] = 1 #increasing
  memory_df$target_area[1] <- FALSE
  memory_df$res_stat[1] <- 2
  #These two are the limits of the log ratio. they are set to +/- infinity
  memory_df$upper[1] <- Inf
  memory_df$lower[1] <- -Inf
  memory_df$upper[2] <- Inf
  memory_df$lower[2] <- -Inf
  #This calculates the percentage change it is set to 1 by default for all values of res_stat =2
  memory_df$perc_change[1] <- 1
  memory_df$perc_change[2] <- 1
  #This is whether the iteration is stable or not
  memory_df$stable[1] <- FALSE
  memory_df$stable[2] <- FALSE
  
  drag_iter<- 1
  
  memory_df$log_ratio[drag_iter+1] <- 1#-memory_df$log_ratio[drag_iter] - (memory_df$error[drag_iter])* 1 *direction_change
  memory_df$common_drag_iter[drag_iter+1] <- 10^( -memory_df$log_ratio[drag_iter+1]) * tstep
  
  #Prep the data before going into the converger
  Prep <- SETSe_data_prep(g = g, 
                          force = force, 
                          flow = flow, 
                          distance = distance, 
                          mass = mass, 
                          edge_name = edge_name,
                          sparse = sparse)
  #The number of iterations has to be smaller than the hyper_iters variable AND the residual static force has to be bigger
  #than the tolerance AND the last two rounds cannot both be stable
  while((drag_iter <= hyper_iters) & 
        ifelse(is.na(memory_df$res_stat[drag_iter]>tol), TRUE, memory_df$res_stat[drag_iter]>tol) & 
        !all(memory_df$stable[c(drag_iter, drag_iter-1)])){
    
    drag_iter <- drag_iter+1      
    # print( memory_df$common_drag_iter[drag_iter])
    embeddings_data <- SETSe_core(
      node_embeddings = Prep$node_embeddings, 
      ten_mat = Prep$ten_mat, 
      non_empty_matrix = Prep$non_empty_matrix, 
      kvect = Prep$kvect, 
      dvect = Prep$dvect, 
      mass = mass,
      tstep = tstep, 
      max_iter = hyper_max, 
      coef_drag = memory_df$common_drag_iter[drag_iter],
      tol = tol,
      sparse = sparse,
      sample = sample) 
    
    node_embeds <- embeddings_data$node_embeddings
    
    memory_df$res_stat[drag_iter] <- sum(abs(node_embeds$static_force))
    
    #set the residual static force to 2 if it exceeds this value
    memory_df$res_stat[drag_iter] <- ifelse(memory_df$res_stat[drag_iter]>2, 2 ,memory_df$res_stat[drag_iter])
    
    #Is the algo in the convex area?
    memory_df$target_area[drag_iter] <- memory_df$res_stat[drag_iter] < 2
    
    #stores a temporary error value
    temp_error <- memory_df$res_stat[drag_iter] - tol
    
    #Log error is negative when small, absolute difference is used to prevent NaNs if the true error is smaller than the tolerance
    memory_df$error[drag_iter] <- log10(abs(temp_error))
    
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
      
      memory_df$perc_change[drag_iter] <- 1 - temp_error/10^memory_df$error[drag_iter-1]
      
    } else {
      memory_df$upper[drag_iter+1] <- memory_df$upper[drag_iter] 
      memory_df$lower[drag_iter+1] <-  memory_df$lower[drag_iter]
      memory_df$best_log_ratio[drag_iter] <- memory_df$best_log_ratio[drag_iter-1]
      
      #The change is twice the hyper tolerance 
      memory_df$perc_change[drag_iter] <- 2*hyper_tol
      
    }
    #complete the stability part of the dataframe. This replaces any infinite values with 2*tol so that an error won't be thrown
    
    memory_df$stable[drag_iter] <- ifelse(is.finite(memory_df$perc_change[drag_iter]), 
                                          memory_df$perc_change[drag_iter],
                                          hyper_tol*2 ) < hyper_tol
    
    ###
    #if you are in the target zone updates happpen adaptively if not using a fixed step size
    ###
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
    
    #Calculate the common drag iteration for this round
    memory_df$common_drag_iter[drag_iter+1] <- 10^( -memory_df$log_ratio[drag_iter+1]) * tstep
    
    
    
    
    if(verbose){
      message_val <- ifelse(memory_df$direction[drag_iter] < 1, "accuracy increasing",  "accuracy decreasing")
      # print(c((drag_iter <= hyper_iters), (memory_df$res_stat[drag_iter]>tol), !all(memory_df$stable[c(drag_iter, drag_iter-1)])))
      print(paste("Iteration", drag_iter, message_val,
                  memory_df$common_drag_iter[drag_iter + 1],
                  "static force", memory_df$res_stat[drag_iter]))
    }
  }
  
  #Keep only value that actually have a residual force
  memory_df <- memory_df %>%
    filter(!is.na(res_stat))
  
  #If the smallest residual force is not less than the tolerance even after the minimum point hase been found
  #Then run the setse algo again using the best log ratio but for the maximum number of iterations
  if(min(memory_df$res_stat)>tol){
    print("Minimum tolerance not exceeded, running SETSe on best parameters")
    embeddings_data <- SETSe_core(
      node_embeddings = Prep$node_embeddings, 
      ten_mat = Prep$ten_mat, 
      non_empty_matrix = Prep$non_empty_matrix, 
      kvect = Prep$kvect, 
      dvect = Prep$dvect, 
      mass = mass,
      tstep = tstep, 
      max_iter = max_iter, 
      coef_drag = memory_df$common_drag_iter[nrow(memory_df)], #uses the best/last coefficient of drag
      tol = tol,
      sparse = sparse,
      sample = sample) 
    
    
  }
  
  
  #Put in edge embeddings
  embeddings_data$edge_embeddings <- calc_tension_strain(g = g,
                                                         embeddings_data$node_embeddings,
                                                         distance = distance, 
                                                         capacity = capacity, 
                                                         flow = flow, 
                                                         edge_name = edge_name, 
                                                         k = k)
  #Add the search record
  embeddings_data$memory_df <- memory_df
  
  return(embeddings_data)
}