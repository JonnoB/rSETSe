Find_network_balance <- function(g, force ="BalencedPower", flow = "PowerFlow", capacity = "Link.Limit",  tstep = 0.5, distance = "Imp", 
                                 mass = 2000, maxIter =2000, frctmultiplier = 1, tol = 1e-10, verbose = TRUE,
                                 TwoNodeSolution = TRUE){
  #needs an edge attribute "distance"
  
  g <- set.edge.attribute(g, "distance", value = get.edge.attribute(g, distance))
  
  A <- as_data_frame(g) %>% 
    select(Link, from, to) %>% 
    gather(key = type, Node, -Link) %>%
    arrange(Node) %>%
    mutate(value = ifelse(type =="from", 1, -1)) %>%
    ungroup %>%
    select(-type) %>%
    spread(key = Node, value, fill = 0) %>%
    arrange(Link)
  
  rowdat <- A$Link
  
  A <- A %>% select(-Link) %>%
    as.matrix()
  
  rownames(A) <-rowdat
  
  rm(rowdat)
  
  
  NodeStatus <- as_data_frame(g, what = "vertices") %>%
    select(node = name, force = force ) %>%
    mutate(
      z = 0,
      mass = mass,
      NetTension = 0, velocity = 0, 
      friction = 0,
      NetForce = force + NetTension - friction,
      acceleration = NetForce/mass,
      t = 0)
  
  Link <- as_data_frame(g)  %>%
    rename(flow = flow,
           capacity = capacity) %>%
    mutate(EdgeName = Link, LL = abs(flow)/capacity, k = Area*E/distance) %>% #This sets a floor and ceiling 
    #to the k values. the more highly loaded a line is the more it should stretch. as LL varies between 0, no loading (stiffness)
    #to 1, overload point, (most elastic). The larger kdiff is the larger the difference in elasticity for highly and lightly loaded lines.
    #Very large kdiff means very little elasticty on lightly loaded lines
    select(EdgeName, distance, k) 
  
  #do special case solution
  if(nrow(Link)==1 & TwoNodeSolution){

    #uses the non-linear optimiser from minpack.lm to find the solution to the two node special case, this is much faster
    solution_angle <- nlsLM(Force ~ ForceV_from_angle(target_angle, k = k, d = d), 
                  start = c(target_angle = pi/4), 
                  data = list(Force = abs(NodeStatus$force[1]), k = Link$k, d = Link$distance), 
                  upper = pi/2)
    
    temp <- NodeStatus %>%
      mutate(z = ifelse(force>0, 
                        tan(coefficients(solution_angle))/2, #height above mid point
                        -tan(coefficients(solution_angle))/2 ),
             acceleration = 0,
             t = 0,
             Delta_acceleration = 0) #height below mid-point
    
    Out <- list(results = temp, 
                NodeList = temp
                  )
    
    
    
  } else{

    Out <- FindStabilSystem2(NodeStatus, A, Link$k, Link$distance, tstep, maxIter, frctmultiplier, tol, verbose = verbose) 
    
  }

  
  return(Out)
  
}
