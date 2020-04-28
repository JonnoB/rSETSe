#'Preapre data for the `SETSe_core` function
#'
#'This is a helper function that makes the `SETSe` function code easier to read.
#'
#'The file outputs a named list containing 
#'
#' @param g An igraph object of the network that is going to be turned into a spring embedding.
#' @param force a character string. The node attribute that contains the force information for the network.
#' @param distance a character string. The edge attribute that contains the distance of the edge.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param edge_name a character string. The edge attribute that contains the names of all the edges.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @export

SETSe_data_prep <-function(g, force, distance, mass, edge_name = edge_name, k = "k", sparse = FALSE){
  #this is a helper function that goes inside the the find network balance function to help make the code easier to read
  
  #just calls distance 'distance' for simplicities sake
  g <- set.edge.attribute(g, "distance", value = get.edge.attribute(g, distance))
  
  #get the adjacency matrix
  Adj <- as_adjacency_matrix(g, sparse = T) #produces a zero 1 adjacency matrix
  #just in case
  diag(Adj) <-0
  
  node_embeddings <- as_data_frame(g, what = "vertices") %>%
    select(node = name, force = force ) %>%
    mutate(
      elevation = 0,
      net_tension = 0, 
      velocity = 0, 
      friction = 0,
      static_force = force + net_tension,
      net_force = static_force - friction,
      acceleration = net_force/mass,
      t = 0,
      Iter = 0)

  #This dataframe is actually a large proportion of the memory requirements of the function
  #However as it is only used by the two node solution it now returns NA unless the graph has two nodes.
  #This should improve efficiency in some cases
  if(ecount(g)==1){
  #Link is used for the two node solution
  Link <- as_data_frame(g)  %>%
    mutate(EdgeName = .data[[edge_name]] #The edge name has to be flexible so .data is used
    ) %>% #This sets a floor and ceiling 
    #to the k values. the more highly loaded a line is the more it should stretch. as LL varies between 0, no loading (stiffness)
    #to 1, overload point, (most elastic). The larger kdiff is the larger the difference in elasticity for highly and lightly loaded lines.
    #Very large kdiff means very little elasticty on lightly loaded lines
    select(EdgeName, distance, kvect = k) %>%
    arrange(EdgeName)
  } else{
    
    Link <- NA
  }
  
  # kmat <- as_adjacency_matrix(g, attr = k, sparse = T)
  # dmat <- as_adjacency_matrix(g, attr = distance, sparse = T)

  #This creates a 4 column matrix with the row and column position, the absolute index and transpose index of the edges in the matrix
  #This means vectors can be used for most operation greatly reducing the amount of memory required and 
  #providing a modest speed increase.
  
  Adj2 <- as(Adj, "dgTMatrix")
  
  non_empty_matrix <-cbind(i = Adj2@i + 1, j = Adj2@j + 1) %>% 
    tibble(#names = rownames(.), 
      rows = .[,1], 
      cols = .[,2],
      index = rows+ (cols-1)*ncol(Adj),
      t_index = cols + (rows-1)*ncol(Adj)) %>% {.[,2:5]} %>%
    as.matrix()

  #extract the non-zero entries of the weighted adjacency matrix for K and distance
  #k can change, whilst distance is effectively a constant
  kvect <- as_adjacency_matrix(g, attr = k, sparse = T)[non_empty_matrix[,3]]

  dvect <- as_adjacency_matrix(g, attr = distance, sparse = T)[non_empty_matrix[,3]]

  #Sparse matrix mode reduces time and memory requirements for larger matrices 100 nodes use dense. matrices of 300or more 300 I'm not sure though
  if(!sparse){
    Adj  <- as.matrix(Adj)} else 
    {
      #converts from sparse matrix of form  dgCMatrix to dsCMatrix 
      Adj <- as(Adj, "symmetricMatrix") 
    }
  
  #The outlist of the variables needed for find_stabil_system
  Out <-     list( node_embeddings,   Link,   Adj,       non_empty_matrix,   kvect,   dvect )
  names(Out) <- c("node_embeddings", "Link", "ten_mat", "non_empty_matrix", "kvect", "dvect")
  
  return(Out)
  
}