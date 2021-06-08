#' SETSe embedding with automatic drag and timestep selection for high-dimensional feature vectors
#'
#' Uses a grid search and a binary search to find appropriate convergence conditions.
#' 
#' @param g An igraph object
#' @param force A character vector. These are the nodes attributes that contain the force the nodes exert on the network.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. This is k for the moment don't change it.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param hyper_iters integer. The hyper parameter that determines the number of iterations allowed to find an acceptable convergence value.
#' @param drag_min integer. A power of ten. The lowest drag value to be used in the search
#' @param drag_max integer. A power of ten. if the drag exceeds this value the tstep is reduced
#' @param tstep_change numeric. A value between 0 and 1 that determines how much the time step will be reduced by default value is 0.5
#' @param hyper_tol numeric. The convergence tolerance when trying to find the minimum value
#' @param hyper_max integer. The maximum number of iterations that SETSe will go through whilst searching for the minimum.
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#'  prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#'  If left blank the static limit is the system absolute mean force.
#' @param verbose Logical. This value sets whether messages generated during the process are suppressed or not.
#' @param include_edges logical. An optional variable on whether to calculate the edge tension and strain. Default is TRUE.
#'  included for ease of integration into the bicomponent functions.
#' @param noisy_termination Stop the process if the static force does not monotonically decrease.
#' 
#' @details This is one of the most commonly used SETSe functions. It automatically selects the convergence time-step and drag values
#'  to ensure efficient convergence.
#' 
#' The noisy_termination parameter is used as in some cases the convergence process can get stuck in the noisy zone of SETSe space.
#' To prevent this the process is stopped early if the static force does not monotonically decrease.  On large networks this 
#' greatly speeds up the search for good parameter values. It increases the chance of successful convergence. 
#' More detail on auto-SETSe can be found in the paper "The spring bounces back" (Bourne 2020).
#' 
#' @return A list of four elements. A data frame with the height embeddings of the network, a data frame of the edge embeddings, 
#' the convergence dynamics dataframe for the network as well as the search history for convergence criteria of the network
#' 
#' @family SETSe
#' @examples
#' 
#' set.seed(234) #set the random see for generating the network
#' 
#' g <- generate_peels_network(type = "E")
#' 
#' g_prep <- g %>%
#' #prepare the network for a binary embedding
#' prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                      force_var = "class", 
#'                      positive_value = "A")
#'                      
#' #embed the network using auto SETSe with default settings
#' embeddings <- SETSe_auto(g_prep)
#' 
#' @export

SETSe_auto_hd <- function(g, 
                       force ="force", 
                       distance = "distance", 
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
                       drag_min = 0.01,
                       drag_max = 100,
                       tstep_change = 0.2,
                       sample = 100,
                       static_limit = NULL,
                       verbose = FALSE,
                       include_edges = TRUE,
                       noisy_termination = TRUE){
  
  #The more negative the log ratio becomes the larger the drag ratio becomes
  
  #mass <- ifelse(is.null(mass), mass_adjuster(g, force = force, resolution_limit = TRUE), mass)
  
  if(verbose){print("prepping dataset")}
  #Prep the data before going into the converger
  Prep <- SETSe_data_prep_hd(g = g, 
                          force = force, 
                          distance = distance, 
                          mass = mass, 
                          edge_name = edge_name,
                          k = k,
                          sparse = sparse)

  #calculate the  default static limit if one is not provided.
  #By default this is the absolute static_force at initiation
  if(is.null(static_limit)){
    static_limit <- force %>%
      purrr::map(~igraph::vertex_attr(g, .x)) %>% unlist %>% {sum(abs(.))}
  }
  
  #This is a safety feature!
  #makes sure the static limit is not smaller than the machine precision. If it is the static limit is set to half
  #the tolerance making the simulation terminate at the first opportunity
  #print(paste("static limit greater than eps", static_limit> .Machine$double.eps^0.5))
  res_stat_limit <- ifelse(static_limit> .Machine$double.eps^0.5, static_limit, tol/2)
  
  #This print out can be deleted, it is only here for debugging reasons
  #print(paste("static_limit", res_stat_limit ))
  
  memory_df<-tibble::tibble(iteration = 1:(2+hyper_iters),
                    error = NA,
                    perc_change = NA,
                    log_ratio = NA,
                    common_drag_iter = NA,
                    tstep = NA,
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
  memory_df$tstep[1] <- tstep
  memory_df$error[1] <- drag_max
  memory_df$best_log_ratio[1] <- 0
  memory_df$direction[1] = 1 #increasing
  memory_df$target_area[1] <- FALSE
  memory_df$res_stat[1] <- ifelse(res_stat_limit<tol, 2*tol, res_stat_limit) #this ensures setse is run at least once
  #The above conduition is only activated on components with very small amounts of force
  #These two are the limits of the log ratio. they are set to +/- infinity
  memory_df$upper[1] <- -log10(drag_max/tstep)
  memory_df$lower[1] <-  Inf
  memory_df$upper[2] <- -log10(drag_max/tstep)
  memory_df$lower[2] <- Inf
  #This calculates the percentage change it is set to 1 by default for all values of res_stat =2
  memory_df$perc_change[1] <- 1
  memory_df$perc_change[2] <- 1
  #This is whether the iteration is stable or not
  memory_df$stable[1] <- FALSE
  memory_df$stable[2] <- FALSE
  
  drag_iter<- 1
  
  #The vector of drag values to search through
  #If the target zone is not found the time step is reduced and the search repeated.
  #This process continues until the target zone is found or hyper_iters < drag_iter
  drag_values <- seq(round(log10(drag_min)), round(log10(drag_max)), by = 1)
  drag_value_index <- 1

  #the tstep value that is passed to setse_core.
  #This value is adapted if no suitable drag coeffeicient is found
  tstep_adapt <- tstep
  
  #The ratio of the drag coefficient over the timestep, the result is then log10'ened
  memory_df$log_ratio[drag_iter+1] <- log10(tstep) - drag_values[drag_value_index]  
  #The first drag value is simply the lowest drag value
  memory_df$common_drag_iter[drag_iter+1] <- 10^drag_values[drag_value_index]
  
  min_error <- which.min(memory_df$error)
  
  if(verbose){print("beginning embeddings search")}
  
  #For the loop to continue the following conditions must all hold
  #1 The number of iterations the loop has gone through has to be smaller than the hyper_iters variable 
  #2 The residual static force has to be bigger than the minimum system tolerance
  #3 The las two rounds cannot both be stable
  
  #The below commented code is useful in some debugging situations
  
  # print(paste("Iters smaller than max",(drag_iter <= hyper_iters),
  #         "res_stat > tol or NA", ifelse(is.na(memory_df$res_stat[drag_iter]>tol), TRUE, memory_df$res_stat[drag_iter]>tol),
  #         "sims unstable", !all(memory_df$stable[c(drag_iter, drag_iter-1)]))
  #       )
  
  # while(drag_iter<5){  
  while((drag_iter <= hyper_iters) &
        ifelse(is.na(memory_df$res_stat[drag_iter]>tol), TRUE, memory_df$res_stat[drag_iter]>tol) &
        !all(memory_df$stable[c(drag_iter, drag_iter-1)])
  ){
    drag_iter <- drag_iter+1    

    #If the drag exceeds the max reset to minimum and divide the time by value x
    #This should never be needed but just in case it is here.
    #can be removed after testing to make everything easier
    if(memory_df$common_drag_iter[drag_iter] > drag_max){
      
      drag_value_index <- 1
      tstep_adapt <- tstep_adapt*tstep_change
      
      #move on to the next drag value and log ratio
      memory_df$log_ratio[drag_iter+1] <- drag_values[drag_value_index] - log10(tstep_adapt)
      
      memory_df$common_drag_iter[drag_iter] <- 10^drag_values[drag_value_index]
      
    } 
    
    # print( memory_df$common_drag_iter[drag_iter])
    embeddings_data <- SETSe_core_hd(
      node_embeddings = Prep$node_embeddings, 
      ten_mat = Prep$ten_mat, 
      non_empty_matrix = Prep$non_empty_matrix, 
      kvect = Prep$kvect, 
      dvect = Prep$dvect, 
      mass = mass,
      tstep = tstep_adapt, 
      max_iter = hyper_max, 
      coef_drag = memory_df$common_drag_iter[drag_iter],
      tol = tol,
      sparse = sparse,
      sample = sample,
      static_limit = static_limit,
      noisy_termination = noisy_termination) 
    
    node_embeds <- embeddings_data$node_embeddings
    
    memory_df$tstep[drag_iter] <- tstep_adapt
    
    memory_df$res_stat[drag_iter] <- sum(abs(node_embeds %>% dplyr::select(dplyr::starts_with("static_force_"))))
    #In certain circumstances SETSe core terminates straight away producing NA values.
    #The below line prevents NA values causing issues by ensureing that the params only produces a max res_stat
    memory_df$res_stat[drag_iter] <- ifelse(is.na(memory_df$res_stat[drag_iter]), res_stat_limit, memory_df$res_stat[drag_iter])
    
    #set the residual static force to res_stat_limit if it exceeds this value
    memory_df$res_stat[drag_iter] <- ifelse(memory_df$res_stat[drag_iter]>res_stat_limit, 
                                            res_stat_limit ,
                                            memory_df$res_stat[drag_iter])
    
    #when noisy_termination =TRUE. processes that terminate due to noisy behaviour are treated as if they had exceeded the
    #static limit
    #note noisy behaviour is a system that does not have a monotonic reduction in static force
    if(noisy_termination & nrow(embeddings_data$network_dynamics)>1){

      monotonic_decrease <- all(embeddings_data$network_dynamics$static_force[2:nrow(embeddings_data$network_dynamics)]<
      embeddings_data$network_dynamics$static_force[1:(nrow(embeddings_data$network_dynamics)-1)])
      
      memory_df$res_stat[drag_iter] <-  ifelse(monotonic_decrease, memory_df$res_stat[drag_iter] , res_stat_limit)
    }
    
    #Is the algo in the convex area?
    memory_df$target_area[drag_iter] <- memory_df$res_stat[drag_iter] < res_stat_limit
    
    #stores a temporary error value
    temp_error <- memory_df$res_stat[drag_iter] - tol
    
    #Log error is negative when small. Absolute difference is used to prevent NaNs if the true error is smaller than the tolerance
    memory_df$error[drag_iter] <- log10(abs(temp_error))
    
    #setting the limits
    #this is triggered when at least one hyperiteration is less than res_stat_limit
    if(sum(memory_df$target_area, na.rm = T)>=1){
      #print("A")
      #Identify row containing the minimum error
      min_error <- which.min(memory_df$error)
      #find log ratio that give the lowest error. there may be more than one
      upper_ratios <- memory_df$log_ratio[memory_df$log_ratio[min_error] >= memory_df$log_ratio] %>% unique()
      #find the log ratio that gives the maximum error aka the res_stat_limit, there are almost certainly more than one of those
      lower_ratios  <- memory_df$log_ratio[memory_df$log_ratio[min_error] <= memory_df$log_ratio] %>% unique()
      
      #This test returns a vector of length zero if the best value is the only unique value.
      #In that case another test needs to be performed.
      upper_check_1 <- !is.finite( upper_ratios[rank(-upper_ratios, ties.method =  "first")==2])
      #this checks the previous check is actually longer than 1
      upper_check_2 <- length(upper_check_1)>0
      upper_result_check <- ifelse(upper_check_2, upper_check_1, FALSE)
      
      #this is the same for the lower section
      lower_check_1 <- !is.finite( lower_ratios[rank(lower_ratios, ties.method =  "first" )==2])
      lower_check_2 <- length(lower_check_1)>0
      lower_result_check <- ifelse(lower_check_2, lower_check_1, FALSE)
      
      memory_df$upper[drag_iter+1] <-  ifelse(upper_result_check, 
                                              # ifelse(!is.finite(upper_result_check), 
                                              -log10(drag_max/tstep), 
                                              upper_ratios[rank(-upper_ratios, ties.method =  "first")==2])  
      
      memory_df$lower[drag_iter+1] <- ifelse(lower_result_check, 
                                             Inf, 
                                             lower_ratios[rank(lower_ratios, ties.method =  "first")==2])
      
      memory_df$best_log_ratio[drag_iter] <- memory_df$log_ratio[min_error]
      
      # percent change in res stat
      memory_df$perc_change[drag_iter] <- (memory_df$res_stat[drag_iter]-min(memory_df$res_stat[1:(drag_iter-1)]))/min(memory_df$res_stat[1:(drag_iter-1)])
      
    } else {
      #print("B")
      #finally if the search has never been in the target zone the upper and lower limits are the same Inf as before
      memory_df$upper[drag_iter+1] <- memory_df$upper[drag_iter] 
      memory_df$lower[drag_iter+1] <-  memory_df$lower[drag_iter]
      memory_df$best_log_ratio[drag_iter] <- memory_df$best_log_ratio[drag_iter-1]
      
      #The change is twice the hyper tolerance 
      memory_df$perc_change[drag_iter] <- 2*hyper_tol
      
    }
    #complete the stability part of the dataframe. This replaces any infinite values with 2*tol so that an error won't be thrown
    
    #If the change in res stat is smaller than the hyper tolerance the system is considered stable.
    #This is only ever used INSIDE the target zone
    memory_df$stable[drag_iter] <- ifelse(is.finite(memory_df$perc_change[drag_iter]), 
                                          abs(memory_df$perc_change[drag_iter]), #The absolute percent change otherwise large negative values count as stable
                                          hyper_tol*2 ) < hyper_tol
    
    ###
    #if you are in the target zone updates happen adaptively if not using a fixed step size
    ###
    if( sum(memory_df$target_area, na.rm = T)>0 ){
      #print("1")
      # print("upper limit found")
      #This makes the next drag ratio the mean of the upper and lower limit.
      #The first time the drag is in the target area the lower limit will be -Inf
      memory_df$log_ratio[drag_iter+1] <- (memory_df$upper[drag_iter+1] + memory_df$lower[drag_iter+1] )/2
      #print(memory_df$log_ratio[drag_iter+1] )
      #if the new log ratio has been used before, break the tie by using the mid point between
      #the best log ratio and the lower log ratio
      memory_df$log_ratio[drag_iter+1] <- ifelse(memory_df$log_ratio[drag_iter+1] %in% memory_df$log_ratio[memory_df$iteration <=drag_iter & 
                                                                                                             memory_df$tstep==tstep_adapt], 
                                                 (memory_df$lower[drag_iter] + memory_df$best_log_ratio[drag_iter])/2,
                                                 memory_df$log_ratio[drag_iter+1] )
      #print(memory_df$log_ratio[drag_iter+1] )
      #If the lower limit has not been found the log ratio will be infinite. 
      #In that case subtract the mean of the best value and the upper limit from the best value otherwise do nothing
      memory_df$log_ratio[drag_iter+1] <- ifelse(is.finite(memory_df$log_ratio[drag_iter+1]),
                                                 memory_df$log_ratio[drag_iter+1],
                                                 memory_df$best_log_ratio[drag_iter] - (memory_df$best_log_ratio[drag_iter] + memory_df$upper[drag_iter] )/2)
      
      #Calculate the common drag iteration for this round
      memory_df$common_drag_iter[drag_iter+1] <- 10^( -memory_df$log_ratio[drag_iter+1]) * tstep_adapt
      #print(memory_df$log_ratio[drag_iter+1] )
    } else {
      #print("2")
      #Before the non trivial upper and lower bound are found simple search is performed to find the target zone
      #fixed rate search across the plateau
      
        drag_value_index <-  drag_value_index + 1
        
        #If the index exceeds the total number of elements in the drag range
        #reset the to minimum drag value and reduce the timstep
 
        #if(10^drag_values[drag_value_index] > drag_max){
        if(drag_value_index > length(drag_values)) {
         
          drag_value_index <- 1
          tstep_adapt <- tstep_adapt*tstep_change
        }
        
      #move on to the next drag value and log ratio
      memory_df$log_ratio[drag_iter+1] <-  log10(tstep_adapt) - drag_values[drag_value_index] 
      
      memory_df$common_drag_iter[drag_iter+1] <- 10^(drag_values[drag_value_index])
    }
    

    if(verbose){
      # message_val <- ifelse(memory_df$res_stat[drag_iter] < memory_df$res_stat[drag_iter-1], 
      #                       "static force decreasing",  
      #                       "static force increasing or stable")
      # print(c((drag_iter <= hyper_iters), (memory_df$res_stat[drag_iter]>tol), !all(memory_df$stable[c(drag_iter, drag_iter-1)])))
      print(paste("Iteration",drag_iter, "drag value",
                  signif(memory_df$common_drag_iter[drag_iter],3),
                  "tstep", memory_df$tstep[drag_iter],
                  #message_val, 
                  "static force is",
                  signif(memory_df$res_stat[drag_iter],3), "target is", signif(tol,3)))
      
     # print(paste("log ratio",memory_df$log_ratio[drag_iter+1]  ,"next drag", memory_df$common_drag_iter[drag_iter+1]))
    }
  }
  
  #Keep only value that actually have a residual force
  memory_df <- memory_df %>%
    dplyr::filter(!is.na(res_stat))
  
  
  #If the smallest residual force is not less than the tolerance even after the minimum point hase been found
  #Then run the setse algo again using the best log ratio but for the maximum number of iterations
  #Ten percent is added to the drag score as the found value is likely to be noisy. This will speed up the convergence
  #If it slows it down by over damping it won't be too bad... I am open to changing this!
  if(min(memory_df$res_stat)>tol){
    
    #if any value has reduced the static_force use the best one.
    #otherwise use the last drag value as it doesn't matter anyway

    if(memory_df$best_log_ratio[nrow(memory_df)]==0){
      
      drag_val <- memory_df$common_drag_iter[nrow(memory_df)]
    } else{
      
      drag_val  <-10^( -memory_df$best_log_ratio[nrow(memory_df)]) * tstep_adapt
      
    }
    
    
    print(paste0("Minimum tolerance not exceeded, running SETSe on best parameters, drag value ", signif(drag_val, 3)))
    embeddings_data <- SETSe_core_time_shift_hd(
      node_embeddings = Prep$node_embeddings, 
      ten_mat = Prep$ten_mat, 
      non_empty_matrix = Prep$non_empty_matrix, 
      kvect = Prep$kvect, 
      dvect = Prep$dvect, 
      mass = mass,
      tstep = tstep_adapt, 
      max_iter = max_iter, 
      coef_drag = drag_val, #uses the best/last coefficient of drag
      tol = tol,
      sparse = sparse,
      sample = sample,
      static_limit = static_limit,
      tstep_change = 0.5,
      dynamic_reset = TRUE) #This is false as something has to be produced by the algo. I can change later 
    
    
    
    
  }
  
  #This is TRUE in almost all cases, but for bicomp is set to false.
  if(include_edges){  
    #Put in edge embeddings
    embeddings_data$edge_embeddings <- calc_tension_strain_hd(g = g,
                                                           embeddings_data$node_embeddings,
                                                           distance = distance, 
                                                           edge_name = edge_name, 
                                                           k = k)
  }
  
  #Add the search record
  embeddings_data$memory_df <- memory_df
  
  return(embeddings_data)
}
