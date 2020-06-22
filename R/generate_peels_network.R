#'Create a random Peel network
#'
#' Creates an example of a network from Peel's quintet of the specified type.
#' 
#' @param type A character which is any of the capital letters A-E
#'  
#' @details  This function generates networks matching the 5 types described in Peel et al 2019(\url{https://doi.org/10.1073/pnas.1713019115}). All networks have 40 nodes,
#'  160 edges, two node classes and four node sub-classes. The connections between the are equal across all 5 types.
#'  As a result all networks generated have identical assortativity. However, as the sub-classes have different connection
#'  probability the structures produced by the networks are very different. When projected into SETSe space the network types
#'  occupy there own area, see Bourne 2020 for details
#'  
#' @return An igraph object that matches one of the 5 Peel's quintet types
#'  
#' @examples
#' g <- generate_peels_network(type = "E")
#' plot(g)
#' 
#' @export

generate_peels_network <- function(type){
  
  #number of connections between classes and sub-classes for each of the types in Peel's quintet
  #There may be a faster way of loading this data than a dataframe
  peel_connection <- structure(list(sub_class_1 = c("A_1", "A_1", "B_1", "A_1", "B_1","A_2", "A_1", "B_1", "A_2", "B_2"), 
                                    sub_class_2 = c("A_1", "B_1", "B_1", "A_2", "A_2", "A_2", "B_2", "B_2", "B_2", "B_2"), 
                                    class_1 = structure(c(1L,1L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L), .Label = c("A", "B"), class = "factor"), 
                                    class_2 = structure(c(1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L), .Label = c("A", "B"), class = "factor"), 
                                    A = c(10, 20,10, 20, 20, 10, 20, 20, 20, 10), 
                                    B = c(38, 20, 0, 2, 20, 0, 20, 2, 20, 38), 
                                    C = c(38, 0, 10, 2, 0, 0, 0, 20, 80, 10),
                                    D = c(10, 0, 10, 20, 0, 10, 0, 20, 80, 10), 
                                    E = c(38,0, 0, 2, 80, 0, 0, 2, 0, 38), 
                                    position.x = c(1, 1, 11, 1, 11, 21, 1, 11, 21, 31), 
                                    position.y = c(1, 11, 11, 21, 21, 21, 31, 31, 31, 31)), 
                               row.names = c(NA, -10L), class = c("tbl_df","tbl", "data.frame"))
  
  
  full_matrix <- matrix(NA, ncol = 40, nrow = 40 )
  submatrix_size <- 10^2 
  numeric_submatrix  <- matrix(1:submatrix_size, nrow = 10)  
  
  
  peel_connection_2 <- peel_connection %>% dplyr::rename(edges = all_of(type))
  
  
  #The loop that randomly samples the adjacency matrix of peels quintet 
  for(n in 1:nrow(peel_connection_2)){
    #in the loop
    # print(n)
    df <- peel_connection_2 %>% dplyr::slice(n)
    
    logic_matrix  <- matrix(FALSE, nrow = 10, ncol = 10)
    
    #if a subclass is being sampled only the upper triangle can be used to prevent double links and self links
    if(df$sub_class_1 == df$sub_class_2){
      sample_vect <-  numeric_submatrix[upper.tri(numeric_submatrix)]
    } else{
      
      sample_vect <- 1:submatrix_size
      
    }
    
    #Convert the sampled matrix elements to TRUE, this will be used for the edges
    logic_matrix[sample(sample_vect, df$edges , replace = FALSE)] <- TRUE
    
    x.start <- dplyr::pull(df, position.x)
    y.start <- dplyr::pull(df, position.y)
    
    full_matrix[x.start:(x.start+ncol(logic_matrix)-1), y.start:(y.start+nrow(logic_matrix)-1)] <- logic_matrix
    
    #END LOOP
  }
  
  full_matrix[lower.tri(full_matrix, diag = TRUE)] <- NA
  
  #classs
  colnames(full_matrix) <- rep(c("A", "B", "A", "B"), c(10, 10, 10, 10))
  #subclass
  rownames(full_matrix) <-  rep(c("A_1", "B_1", "A_2", "B_2"), c(10, 10, 10, 10))
  
  spec_g <-igraph::graph_from_adjacency_matrix(full_matrix, mode = "undirected", diag = FALSE, add.colnames ="class", add.rownames = "sub_class")
  
  return(spec_g)
  
  
}
