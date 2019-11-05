#'Preapre data for the `find network balance` function
#'
#'This is a helpfer function that makes `find_network_balance` code easier to read. It prepares the data into three datasets.
#'Whether stuff is a character string or not needs to be double checked!
#'
#'The file outputs a named list containing 
#'
#'@param g An igraph object of the network that is going to be turned into a spring embedding.
#'@param force a character string. The node attribute that contains the force information for the network.
#'@param distance a character string. The edge attribute that contains the distance of the edge.
#'@param mass a character string. the node attribute that contains the mass of the nodes.
#'@param edge_name a character string. The edge attribute that contains the names of all the edges.

Prepare_data_for_find_network_balance <-function(g, force, flow, distance, mass, edge_name = edge_name){
  #this is a helper function that goes inside the the find network balance function to help make the code easier to read
  
  g <- set.edge.attribute(g, "distance", value = get.edge.attribute(g, distance))
  
  A <- as_data_frame(g) %>% 
    select(edge_name, from, to) %>% 
    gather(key = type, Node, -edge_name) %>%
    arrange(Node) %>%
    mutate(value = ifelse(type =="from", 1, -1)) %>%
    ungroup %>%
    select(-type) %>%
    spread(key = Node, value, fill = 0) %>%
    arrange(edge_name)
  
  rowdat <- A$edge_name
  
  A <- A %>% select(-edge_name) %>%
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
    mutate(EdgeName = .data[[edge_name]], #The edge name has to be flexible so .data is used
           k = Area*E/distance) %>% #This sets a floor and ceiling 
    #to the k values. the more highly loaded a line is the more it should stretch. as LL varies between 0, no loading (stiffness)
    #to 1, overload point, (most elastic). The larger kdiff is the larger the difference in elasticity for highly and lightly loaded lines.
    #Very large kdiff means very little elasticty on lightly loaded lines
    select(EdgeName, distance, k) %>%
    arrange(EdgeName)
  
  Out <- list(NodeStatus, Link, A)
  names(Out) <- c("NodeStatus", "Link", "A")
  
  return(Out)
  
}