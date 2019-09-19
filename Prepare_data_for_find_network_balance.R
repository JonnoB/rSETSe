Prepare_data_for_find_network_balance <-function(g, force, flow, distance, mass){
  #this is a helper function that goes inside the the find network balance function to help make the code easier to read
  
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
      t = 0) %>%
    arrange(node)
  
  Link <- as_data_frame(g)  %>%
    rename(flow = flow) %>%
    mutate(EdgeName = Link, k = Area*E/distance) %>% #This sets a floor and ceiling 
    #to the k values. the more highly loaded a line is the more it should stretch. as LL varies between 0, no loading (stiffness)
    #to 1, overload point, (most elastic). The larger kdiff is the larger the difference in elasticity for highly and lightly loaded lines.
    #Very large kdiff means very little elasticty on lightly loaded lines
    select(EdgeName, distance, k) %>%
    arrange(EdgeName)
  
  Out <- list(NodeStatus, Link, A)
  names(Out) <- c("NodeStatus", "Link", "A")
  
  return(Out)
  
}