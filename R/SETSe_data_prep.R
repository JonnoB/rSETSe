#'Preapre data for the `SETSe_core` function
#'
#'This is a helper function that makes the `SETSe` function code easier to read. Seldom used otherwise
#'
#'The file outputs a named list containing 
#'
#' @param g An igraph object of the network that is going to be turned into a spring embedding.
#' @param force a character string. The node attribute that contains the force information for the network.
#' @param distance a character string. The edge attribute that contains the distance of the edge.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param edge_name a character string. The edge attribute that contains the names of all the edges.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param k A character string. This is k for the moment don't change it.
#' 
#' @details The output of the function is different depending on the number of nodes in the network. This optimises the
#' operations performed so that networks that will be solved using the two node solution only have data for the two node solution
#' whilst networks that will be solved using the SETSe algorithm do not have data required for the two node solution.
#' This minimises time and memory, which is useful when large numbers of data preparations need to be performed, for example
#' when running setse_bicomp() on a large network.
#'
#' @return Returns a list that contains all the parts needed to allow the SETSe family of functions and the two_node_solution
#' function to produce embeddings.
#' 
#' @noRd

SETSe_data_prep  <-function(g, force, distance, mass, edge_name = edge_name, k = "k", sparse = FALSE){
  #this is a helper function that goes inside the the find network balance function to help make the code easier to read
  
  #just calls distance 'distance' for simplicities sake
  g <- igraph::set.edge.attribute(g, "distance", value = igraph::get.edge.attribute(g, distance))
  
  #node_embeddings <- as_data_frame(g, what = "vertices") %>%
  #select(node = name, force = force ) %>%
  # 
  
  #
  # These commmented out chunks are me trying to remove compilation errors that show up when I am using pkgdown
  #
  # 
  # total_nodes <- igraph::vcount(g)
  # 
  # node_embeddings <-  data.frame(node = igraph::get.vertex.attribute(g, "name"),
  #                                force =  igraph::get.vertex.attribute(g, force),
  #                                elevation = rep(0, total_nodes),
  #                                net_tension = rep(0, total_nodes),
  #                                velocity = rep(0, total_nodes),
  #                                friction = rep(0, total_nodes))
    
  # node_embeddings <-  data.frame(node = igraph::get.vertex.attribute(g, "name"),
  #                                force =  igraph::get.vertex.attribute(g, force),
  #                                elevation = 0,
  #                                net_tension = 0,
  #                                velocity = 0,
  #                                friction = 0)
  node_embeddings <-  data.frame(node = igraph::get.vertex.attribute(g, name = "name"))
  
  node_embeddings$force <- igraph::get.vertex.attribute(g, name = force)
  node_embeddings$elevation <-0
  node_embeddings$net_tension <-0
  node_embeddings$velocity <-0
  node_embeddings$friction <-0

  
  node_embeddings$static_force <- node_embeddings$force
  node_embeddings$net_force <- node_embeddings$static_force
  node_embeddings$acceleration <- node_embeddings$net_force/mass
  node_embeddings$t <- 0
  node_embeddings$Iter <- 0
  
  #This dataframe is actually a large proportion of the memory requirements of the function
  #However as it is only used by the two node solution it now returns NA unless the graph has two nodes.
  #This should improve efficiency in some cases
  if(igraph::ecount(g)==1){
    #Link is used for the two node solution
    Link <- igraph::as_data_frame(g) 
    Link$EdgeName <- Link[,edge_name]
    Link <- Link[,c("EdgeName", "distance", "k")] #%>%
    # arrange(EdgeName) The arrange is removed as the edges are already ordered alphabetically
    
    Out <-     list( node_embeddings,   Link,   NA,       NA,   NA,   NA )
    
  } else{
    #Under normal circumstances this is the code that is run.
    #keepinig normal data-prep seperate from the two node prep, speeds up both cases.
    #This is especially useful for the biconnected components version which can get bogged down when there are many
    #two node solutions to be solved.
    
    #get the adjacency matrix
    Adj <- igraph::as_adjacency_matrix(g, sparse = T) #produces a zero 1 adjacency matrix
    #just in case
    diag(Adj) <-0
    
    #This creates a 4 column matrix with the row and column position, the absolute index and transpose index of the edges in the matrix
    #This means vectors can be used for most operation greatly reducing the amount of memory required and 
    #providing a modest speed increase.
    
    Adj2 <- methods::as(Adj, "dgTMatrix")
    
    non_empty_matrix <-cbind(i = Adj2@i + 1, j = Adj2@j + 1) %>% 
      tibble::tibble(#names = rownames(.), 
        rows = .[,1], 
        cols = .[,2],
        index = rows+ (cols-1)*ncol(Adj),
        t_index = cols + (rows-1)*ncol(Adj)) %>% {.[,2:5]} %>%
      as.matrix()
    
    #extract the non-zero entries of the weighted adjacency matrix for K and distance
    #k can change, whilst distance is effectively a constant
    kvect <- igraph::as_adjacency_matrix(g, attr = k, sparse = T)[non_empty_matrix[,3]]
    
    dvect <- igraph::as_adjacency_matrix(g, attr = distance, sparse = T)[non_empty_matrix[,3]]
    
    #Sparse matrix mode reduces time and memory requirements for larger matrices 100 nodes use dense. matrices of 300or more 300 I'm not sure though
    if(!sparse){
      Adj  <- as.matrix(Adj)} else 
      {
        #converts from sparse matrix of form  dgCMatrix to dsCMatrix 
        Adj <- methods::as(Adj, "symmetricMatrix") 
      }
    
    #The outlist of the variables needed for find_stabil_system
    Out <-     list( node_embeddings,   NA,   Adj,       non_empty_matrix,   kvect,   dvect )
  }
  
  names(Out) <- c("node_embeddings", "Link", "ten_mat", "non_empty_matrix", "kvect", "dvect")
  
  return(Out)
  
}
