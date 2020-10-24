#' two node solution
#' 
#' Uses a newton-raphson to solve the two node solution. Seldom used outside of a 'SETSe' function
#' 
#' @details a helper function inside the main SETSe suite of functions. However it can be used to solve two node graphs
#' 
#' @param g the graph must be two nodes connected by a single edge
#' @param Prep The output of SETSe_data_prep. provides the inputs needed to do the two node solution
#' @param auto_setse_mode outputs an additional list element "memory_df" to make it compatible with the auto_SETSe function
#' 
#' @return 
#' A list of two elements node_embeddings and network_dynamics. if `auto_setse_mode==TRUE` 
#' then a third element is returned the memory_df dataframe
#' 
#' @export
#' 
two_node_solution <- function(g, Prep = Prep, auto_setse_mode = FALSE){
  
  start_time <- Sys.time()
  #If the the force of the two nodes is effectively identical then the solution angle is zero.
  #This also covers the case of the forces being close to 0
  #It should be noted this means there should be a minimum force in the system to account for the nodes.
  #Very large networks may have small force values if normalised to one.
  if(isTRUE(all.equal(Prep$node_embeddings$force[1], Prep$node_embeddings$force[2]))){
    
    solution_angle <-0
    
  } else {
    
    
    #uses the non-linear optimiser from minpack.lm to find the solution to the two node special case, this is much faster
    solution_angle <- minpack.lm::nlsLM(Force ~ ForceV_from_angle(target_angle, k = k, d = d), 
                            start = c(target_angle = pi/4), 
                            data = list(Force = abs(Prep$node_embeddings$force[1]), k = Prep$Link$k, d = Prep$Link$distance), 
                            upper = pi/2) %>% stats::coefficients()      
    

    
  }
  
  stop_time <- Sys.time()
  
  temp <- Prep$node_embeddings #%>%
  # mutate(elevation = ifelse(force>0, 
  #                           tan(solution_angle)/2, #height above mid point
  #                           -tan(solution_angle)/2 ), #height below mid-point
  #        net_force = 0,
  #        acceleration = 0,
  #        #         k      * the extension of the edge due to stretching      * 
  #        net_tension = ifelse(force>0, 
  #                             -abs(Prep$node_embeddings$force[1]), 
  #                             abs(Prep$node_embeddings$force[1]) )
  # ) 
  temp$elevation <-ifelse(temp$force>0, 
                          tan(solution_angle)/2, #height above mid point
                          -tan(solution_angle)/2 )
  temp$net_force <- 0
  temp$acceleration <- 0
  temp$net_tension <- ifelse(temp$force>0, 
                             -abs(Prep$node_embeddings$force[1]), 
                             abs(Prep$node_embeddings$force[1])
  )

  
  Out <- list(network_dynamics = tibble::tibble(t = 0, 
                                        Iter = 0,
                                        static_force = 0, 
                                        kinetic_force = 0), 
              node_embeddings = temp,
              time_taken = tibble::tibble(time_diff = stop_time - start_time, nodes = 2, edges = 1)
  )
  
  if(auto_setse_mode){
    #the memory_df data frame
    #everything is basically NA but it allows easier post processing
    Out$memory_df<-tibble::tibble(iteration = 1,
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
  }
  
  return(Out)
  
}