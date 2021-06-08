#'Preapre data for the `SETSe_core_hd` function for high dimensional feature networks
#'
#'This is a helper function that makes the `SETSe` function code easier to read. Seldom used otherwise
#'
#'The file outputs a named list containing 
#'
#' @param g An igraph object of the network that is going to be turned into a spring embedding.
#' @param force a vector of character strings. THe node attributes that contain the force information for the network.
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

SETSe_data_prep_hd  <-function(g, force, distance, mass, edge_name = edge_name, k = "k", sparse = FALSE){
  #this is a helper function that goes inside the the find network balance function to help make the code easier to read
  
  #just calls distance 'distance' for simplicities sake
  g <- igraph::set.edge.attribute(g, "distance", value = igraph::get.edge.attribute(g, distance))
  
  node_embeddings <-  data.frame(node = igraph::get.vertex.attribute(g, name = "name"))
  
  #create the matrix of forces
  force_tib <- force %>% purrr::map(~{
    igraph::get.vertex.attribute(g, name = .x)})
  names(force_tib) <- force
  force_tib <- tibble::as_tibble(force_tib)
  
  #create an empty tibble so that the other SETSe elements can be added for each dimension
  empty_tib <- force_tib
  empty_tib[] <- 0
  
  node_embeddings <- dplyr::bind_cols(node_embeddings, force_tib%>% dplyr::rename_with(., .fn = ~paste0("force_", .))) 
  node_embeddings <- dplyr::bind_cols(node_embeddings, empty_tib %>% dplyr::rename_with(., .fn = ~paste0("elevation_", .)))
  node_embeddings <- dplyr::bind_cols(node_embeddings, empty_tib %>% dplyr::rename_with(., .fn = ~paste0("net_tension_", .)))
  node_embeddings <- dplyr::bind_cols(node_embeddings, empty_tib %>% dplyr::rename_with(., .fn = ~paste0("velocity_", .)))
  node_embeddings <- dplyr::bind_cols(node_embeddings, empty_tib %>% dplyr::rename_with(., .fn = ~paste0("friction_", .)))
  
  
  node_embeddings <- dplyr::bind_cols(node_embeddings, force_tib %>% dplyr::rename_with(., .fn = ~paste0("static_force_", .)))
  node_embeddings <- dplyr::bind_cols(node_embeddings, force_tib %>% dplyr::rename_with(., .fn = ~paste0("net_force_", .)))
  node_embeddings <- dplyr::bind_cols(node_embeddings, {force_tib/mass} %>% dplyr::rename_with(., .fn = ~paste0("acceleration_", .)))
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
