#' Auto SETSe
#'
#' Uses a grid search and a binary search to find appropriate convergence conditions.
#' 
#' @details This function is pretty useful, it takes advantage of the linear relationship between the timestep and the coefficient of drag
#' to search along the log linear line formed by tstep/coef_drag to find the convergence conditions.
#' 
#' @param g An igraph object
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
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
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is the system absolute mean force.
#' @param inlcude_edges logical. An optional variable on wehther to calculate the edge tension and strain. Default is TRUE.
#'  included for ease of integration into the bicomponent functions.
#' @return A list of four elements. A data frame with the height embeddings of the network, a data frame of the edge embeddings, 
#' the convergence dynamics dataframe for the network as well as the search history for convergence criteria of the network
#' @seealso \code{\link{auto_SETSe}} \code{\link{SETSe_bicomp}}
#'  @examples
#' set.seed(234) #set the random see for generating the network
#' g <- generate_peels_network(type = "E")
#' embeddings <- g %>%
#' #prepare the network for a binary embedding
#' prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                      force_var = "class", 
#'                      positive_value = "A") %>%
#' #embed the network using auto setse
#'   auto_SETSe()
#' 
#' @export

auto_SETSe <- function(g, 
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
                       step_size = 1,
                       sample = 100,
                       static_limit = NULL,
                       verbose = FALSE,
                       include_edges = TRUE){
  
  
  if(verbose){print("prepping dataset")}
  #Prep the data before going into the converger
  Prep <- SETSe_data_prep(g = g, 
                          force = force, 
                          distance = distance, 
                          mass = mass, 
                          edge_name = edge_name,
                          k = k,
                          sparse = sparse)
  
  if(is.null(static_limit)){
    static_limit <- sum(abs(vertex_attr(g, force)))
  }
  
  #This is a safety feature!
  #makes sure the static limit is not smaller than the machine precision. If it is the static limit is set to half
  #the tolerance making the simulation terminate at the first opportunity
  #print(paste("static limit greater than eps", static_limit> .Machine$double.eps^0.5))
  res_stat_limit <- ifelse(static_limit> .Machine$double.eps^0.5, static_limit, tol/2)
  
  #This print out can be deleted, it is only here for debugging reasons
  #print(paste("static_limit", res_stat_limit ))
  
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
  memory_df$res_stat[1] <- ifelse(res_stat_limit<tol, 2*tol, res_stat_limit) #this ensures setse is run at least once
  #The above conduition is only activated on components with very small amounts of force
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
  
  min_error <- which.min(memory_df$error)
  
  if(verbose){print("beginning embeddings search")}
  #The number of iterations has to be smaller than the hyper_iters variable AND the residual static force has to be bigger
  #than the tolerance AND the last two rounds cannot both be stable
  
  #The below commented code is useful in some debugging situations
  
  # print(paste("Iters smaller than max",(drag_iter <= hyper_iters),
  #         "res_stat > tol or NA", ifelse(is.na(memory_df$res_stat[drag_iter]>tol), TRUE, memory_df$res_stat[drag_iter]>tol),
  #         "sims unstable", !all(memory_df$stable[c(drag_iter, drag_iter-1)]))
  #       )
  
 # while(drag_iter<9){  
  while((drag_iter <= hyper_iters) &
        ifelse(is.na(memory_df$res_stat[drag_iter]>tol), TRUE, memory_df$res_stat[drag_iter]>tol) &
        !all(memory_df$stable[c(drag_iter, drag_iter-1)])
        ){

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
      sample = sample,
      static_limit = static_limit) 
    
    node_embeds <- embeddings_data$node_embeddings
    
    memory_df$res_stat[drag_iter] <- sum(abs(node_embeds$static_force))
    
    #set the residual static force to res_stat_limit if it exceeds this value
    memory_df$res_stat[drag_iter] <- ifelse(memory_df$res_stat[drag_iter]>res_stat_limit, 
                                            res_stat_limit ,
                                            memory_df$res_stat[drag_iter])
    
    #Is the algo in the convex area?
    memory_df$target_area[drag_iter] <- memory_df$res_stat[drag_iter] < res_stat_limit
    
    #stores a temporary error value
    temp_error <- memory_df$res_stat[drag_iter] - tol
    
    #Log error is negative when small. Absolute difference is used to prevent NaNs if the true error is smaller than the tolerance
    memory_df$error[drag_iter] <- log10(abs(temp_error))
    
    #setting the limits
    #this is triggered when at least one hyperiteration is less than res_stat_limit
    if(sum(memory_df$target_area, na.rm = T)==1){
      
      #The first time the static force is in the target zone, the lower bound of drag is known, as it has to be the previous value
      #The upper bound is not yet known and can still be infinity
      
      memory_df$upper[drag_iter+1] <- Inf
      memory_df$lower[drag_iter+1] <- memory_df$log_ratio[which(memory_df$target_area)-1]
      
      memory_df$best_log_ratio[drag_iter] <- memory_df$log_ratio[min_error]
      
      # percent change in res stat
      memory_df$perc_change[drag_iter] <- (memory_df$res_stat[drag_iter]-min(memory_df$res_stat[1:(drag_iter-1)]))/min(memory_df$res_stat[1:(drag_iter-1)])
      
    } else if(sum(memory_df$target_area, na.rm = T)>1){
      #The next and subsequent iterations the upper and lower limit will be known and one of them will be adjusted as 
      #the search continues
      
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
      upper_result_check <- ifelse(upper_check_2, upper_check_1, Inf)
      
      #this is the same for the lower section
      lower_check_1 <- lower_ratios[rank(lower_ratios, ties.method =  "first" )==2]
      lower_check_2 <- length(lower_check_1)>0
      lower_result_check <- ifelse(lower_check_2, lower_check_1, Inf)
      
      memory_df$upper[drag_iter+1] <-   ifelse(!is.finite(upper_result_check), 
                                               -Inf, 
                                               upper_ratios[rank(-upper_ratios, ties.method =  "first")==2])  
      # memory_df$upper[drag_iter+1] <-   ifelse(sum(is.finite(upper_ratios))==0, 
      #                                          -Inf, 
      #                                          upper_ratios[rank(-upper_ratios, ties.method =  "first")==2])
      
      memory_df$lower[drag_iter+1] <- ifelse(!is.finite(lower_result_check), 
                                             Inf, 
                                             lower_ratios[rank(lower_ratios, ties.method =  "first")==2])
      
      memory_df$best_log_ratio[drag_iter] <- memory_df$log_ratio[min_error]
      
      # percent change in res stat
      memory_df$perc_change[drag_iter] <- (memory_df$res_stat[drag_iter]-min(memory_df$res_stat[1:(drag_iter-1)]))/min(memory_df$res_stat[1:(drag_iter-1)])
      
    } else {
      
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
    if( is.finite(memory_df$upper[drag_iter+1]) ){
      # print("upper limit found")
      memory_df$log_ratio[drag_iter+1] <- (memory_df$upper[drag_iter+1] + memory_df$lower[drag_iter+1] )/2

      #if the logratios are the same break by choosing the lower limit and the best value
      memory_df$log_ratio[drag_iter+1] <- ifelse(memory_df$log_ratio[drag_iter+1] %in%  memory_df$log_ratio[1:drag_iter], 
                                                 (memory_df$lower[drag_iter] + memory_df$best_log_ratio[drag_iter])/2,
                                                 memory_df$log_ratio[drag_iter+1] )
      
    } else {
      #Before the non trivial upper and lower bound are found simple search is performed to find the target zone
      #fixed rate search across the plateau
      memory_df$log_ratio[drag_iter+1] <- memory_df$log_ratio[drag_iter] - step_size*memory_df$direction[drag_iter]
      
    }
    
    #Calculate the common drag iteration for this round
    memory_df$common_drag_iter[drag_iter+1] <- 10^( -memory_df$log_ratio[drag_iter+1]) * tstep
    
    
    
    
    if(verbose){
      message_val <- ifelse(memory_df$res_stat[drag_iter] < memory_df$res_stat[drag_iter-1], 
                            "static force decreasing",  
                            "static force increasing or stable")
      # print(c((drag_iter <= hyper_iters), (memory_df$res_stat[drag_iter]>tol), !all(memory_df$stable[c(drag_iter, drag_iter-1)])))
      print(paste("Iteration",drag_iter, "drag value",
                signif(memory_df$common_drag_iter[drag_iter + 1],3),
                  message_val, signif(memory_df$res_stat[drag_iter],3)))
    }
  }
  
  #Keep only value that actually have a residual force
  memory_df <- memory_df %>%
    filter(!is.na(res_stat))
  
  #If the smallest residual force is not less than the tolerance even after the minimum point hase been found
  #Then run the setse algo again using the best log ratio but for the maximum number of iterations
  #Ten percent is added to the drag score as the found value is likely to be noisey. This will speed up the convergence
  #If it slows it down by over damping it won't be too bad... I am open to changing this!
  if(min(memory_df$res_stat)>tol){
    
    #if any value has reduced the static_force use the best one.
    #otherwise use the last drag value as it doesn't matter anyway
    if(memory_df$best_log_ratio[nrow(memory_df)]==0){
      
      drag_val <- memory_df$common_drag_iter[nrow(memory_df)]
    } else{
      
      drag_val  <-10^( -memory_df$best_log_ratio[nrow(memory_df)]) * tstep 
      drag_val <- drag_val * 1.1 #add tem percent to get smooth convergence
      
    }
    
    
    print(paste0("Minimum tolerance not exceeded, running SETSe on best parameters, drag value ", drag_val))
    embeddings_data <- SETSe_core(
      node_embeddings = Prep$node_embeddings, 
      ten_mat = Prep$ten_mat, 
      non_empty_matrix = Prep$non_empty_matrix, 
      kvect = Prep$kvect, 
      dvect = Prep$dvect, 
      mass = mass,
      tstep = tstep, 
      max_iter = max_iter, 
      coef_drag = drag_val, #uses the best/last coefficient of drag
      tol = tol,
      sparse = sparse,
      sample = sample,
      static_limit = static_limit) 
    
    
  }
  
  #This is TRUE in almost all cases, but for bicomp is set to false.
  if(include_edges){  
    #Put in edge embeddings
    embeddings_data$edge_embeddings <- calc_tension_strain(g = g,
                                                           embeddings_data$node_embeddings,
                                                           distance = distance, 
                                                           edge_name = edge_name, 
                                                           k = k)
  }
  
  #Add the search record
  embeddings_data$memory_df <- memory_df
  
  return(embeddings_data)
}