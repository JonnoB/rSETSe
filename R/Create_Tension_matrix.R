#' Create tension matrix
#' 
#' Calculates the edge tension in the network. This is a helper function and is not called directly. Ideally the entire thing will be
#' rewritten in C++
#' 
#' @param ten_mat A dataframe. represents the edge node matrix
#' @param zvect A numeric vector. The node elevation, this has one element for every edge in the matrix
#' @param zvect_t A numeric vector. The node elevation transpose, this has one element for every edge in the matrix
#' @param kvect A numeric vector. The vector of spring coefficients between connected nodes, this has one element for every edge in the matrix
#' @param dvect A numeric vector. The original distance between adjacent nodes, this has one element for every edge in the matrix
#' @param non_empty_index non_empty_matrix A numeric matrix. This contains the node indexes in the adjacency matrix it has 4 columns
#'  the row, the column, the absolute index the transpose of the absolute index.
#'  
Create_Tension_matrix <- function(ten_mat, zvect, zvect_t, dvect, kvect, non_empty_index){

    #this takes the zvect and effectively subtracts the matrix from its transpose. However, a vector of the non-zero cells only is used
  #As the matrix is sparse this greatly reduces the calculations
  dzvect <- zvect_t - zvect #The difference in height between adjacent nodes 
  
  #the hypotenuse of the spring distance triangle
  Hvect <- sqrt(dzvect^2 + dvect^2)
  
  #the tension vector. the dZvect/Hvect is the vertical component of the tension
  ten_mat[non_empty_index] <- kvect*(Hvect-dvect)*dzvect/Hvect
  #arma::mat A = A.replace(datum::nan, 0); #Armadillo code to replace nans caused by dividing by 0. The sparse matrix question remains though
  return (ten_mat)
}
