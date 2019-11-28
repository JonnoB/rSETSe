#' Calculate system dynamics
#' 
#' This is a helper function that calculates the dynamics at time t of all the nodes in the system.
#' It is seldom called on it's own but is called by other functions.
#' 
#' The output of the function is a dataframe that is the next iteration of the NodeStatus dataframe the function recieves.
#' 
#' @param NodeStatus A data frame The current dynamics and forces experienced by the node a data frame.
#' @param EdgeNode 
#' @param kvect A numeric vector of the spring stiffnesses
#' @param dvect A numeric vector of the initial distances between the nodes
#' @param tstep A numeric value. The time step, measured in seconds, that will be used to calculate the new dynamic state
#' @param frctmultiplier A numeric value. Used to set a multiplier on the friction value. Generally leave this alone.
#'
#' @export

Calc_System_Dynamics <- function(NodeStatus, EdgeNode, kvect, dvect, tstep = 1, frctmultiplier = 1){
  
  Tension_vect <- Create_Tension_matrix(EdgeNode, NodeStatus$z, kvect, dvect) %>% rowSums()
  Friction_vect <- Calc_Damping_matrix(EdgeNode, NodeStatus$velocity, kvect, NodeStatus$mass) %>% rowMeans()
  
  NodeStatus2 <- NodeStatus %>%
    mutate(z = Distance(z, velocity, acceleration, t0 = t, t1 = t + tstep),
           velocity = velocity(velocity, acceleration, t, t + tstep),
           NetTension =  Tension_vect,
           friction = Friction_vect*frctmultiplier, #velocity*10,
           NetForce = force + NetTension - friction,
           acceleration = NetForce/mass,
           Delta_acceleration = (acceleration-NodeStatus$acceleration)/tstep, #Find change in acceleration, used in early termination
           t = t + tstep)
  
  paste("acceleration", sum(abs(NodeStatus2$acceleration)), "velocity", sum(abs(NodeStatus2$velocity)))
  
  return(NodeStatus2)
  
}