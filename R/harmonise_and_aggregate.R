#'Harmonise and aggregate
#'
#'A helper function of create_stabilised_blocks_expanded poss also create_stabilised_blocks
#'This is used to quickly aggregate the biconnected component values removing double nodes
#'It is much faster than using the dplyr summarise function.
#'It works as the position of all the nodes are already known.
#'The function outputs a vector
#'
#' @param vect a numeric vector, this is the vector that requires harmonising and aggregating 
#' @param component_adjust_mat a numeric matrix. the output of the adjust_components function
#' @param harmonise_vector logical. should the vector be harmonised or just aggregated
#' @param colsum a logical. Either the function does a col sum or a col means.
#'
harmonise_and_aggregate <- function(vect, component_adjust_mat, harmonise_vector  = TRUE, colsum = TRUE){
  #create a vector of the new harmonised values
if(harmonise_vector){  
  
  harmonised <- vect +  as.vector(component_adjust_mat$ceiling %*% vect) - 
    as.vector(component_adjust_mat$floor %*% vect)
} else{
  
  harmonised <- vect
  
}
  #create matrix of the sum of harmonised node_component values per node.
  harmonised_mat <- {Diagonal(x = harmonised) %*% component_adjust_mat$aggregation_matrix} 
  
  if(is.null(colsum)){
    #pass through useful when aggregation is done elsewhere
    Out <- harmonised
    
  } else if(colsum){
    
    Out <- harmonised_mat %>% Matrix::colSums(.)
    
  } else{
    
    Out <- harmonised_mat %>% Matrix::colMeans(.)
    
  } 
  
  return(Out)
  
}