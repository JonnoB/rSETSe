#' Calculate line strain
#' 
#' This function calculates the line strain characteristics for a graph. It ia a helper function of 
#' find_stabil_system and is rarely called by the ser. It is a much faster version of'calc_tension_strain'.
#' 
#' @param edge_mat Numeric matrix. 
#' @param solved_height_mat A numeric matrix. This matrix is the result of the Calc_system_dynamics functions
#' @param merge_order_mat A numeric matrix. This provides the order that the vectors will be provided in 
#' 

Calc_line_strain <- function(edge_mat, solved_height_mat, merge_order_mat){
  
  ##
  #This is the matrix version
  ##
  edge_mat[,2] <- solved_height_mat[merge_order_mat[,1],2]
  edge_mat[,3] <- solved_height_mat[merge_order_mat[,2],2]
  
  
  #base is faster than dplyr, when it comes to lots of iterations it adds up
  edge_mat[,4] <- abs(edge_mat[,2]-edge_mat[,3])
  edge_mat[,5] <- sqrt(edge_mat[,4]^2 + edge_mat[,1]^2)
  edge_mat[,6]<- (edge_mat[,5]-edge_mat[,1])/edge_mat[,1]
  
  return(edge_mat)
  
}