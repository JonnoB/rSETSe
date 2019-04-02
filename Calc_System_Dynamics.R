Calc_System_Dynamics <- function(NodeStatus, EdgeNode, kvect, dvect, tstep = 1, frctmultiplier = 1){
  
  Tension_vect <- Create_Tension_matrix(EdgeNode, NodeStatus$z, kvect, dvect) %>% rowSums()
  Friction_vect <- Calc_Damping_matrix(EdgeNode, NodeStatus$velocity, kvect, NodeStatus$mass) %>% rowMeans()
  
  NodeStatus2 <- NodeStatus %>%
    mutate(z = distance(z, velocity, acceleration, t0 = t, t1 = t + tstep),
           velocity = velocity(velocity, acceleration, t, t + tstep),
           NetTension =  Tension_vect,
           friction = Friction_vect*frctmultiplier, #velocity*10,
           NetForce = force + NetTension - friction,
           acceleration = NetForce/mass,
           Delta_acceleration = acceleration-NodeStatus$acceleration, #Find change in acceleration, used in early termination
           t = t + tstep)
  
  paste("acceleration", sum(abs(NodeStatus2$acceleration)), "velocity", sum(abs(NodeStatus2$velocity)))
  
  return(NodeStatus2)
  
}