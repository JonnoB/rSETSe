#' Create tension matrix
#' 
#' Calculates the edge tension in the network. This is a helper function and is seldom called directly
#' @param Edgenode A dataframe. represents the edge node matrix
#' @param zvect A numeric vector. The node elevation.
#' @param kvect A numeric vector. The vector of spring coefficients between connected nodes.
#' @param dvect A numeric vector. The original distance between adjacent nodes.
#' 
#' @export
Create_Tension_matrix <- function(EdgeNode, zvect, kvect, dvect){
  #Creates a matrix showing the forces exerted by the contracting spring for each edge of a node
  #positive numbers pull down (like the force of mg) and negative forces pull up
  #The v_vect, kvect and mvect all have to be ordered alphanumerically!
  #EdgeNode <- The edgenode matrix 
  #zvect the height of each node
  #kvect the vector of spring stiffness for each edge
  #dvect, the vector of horizontal distance between the springs.
  
  #example
  # SmallEdgeNode <- matrix(c(1,-1,0,
  #                           0,1,-1), nrow = 2, byrow = T)
  # Smallzvect <- c(2,-1,1)
  # Smallkvect <- c(3,4)
  # Smalldvect <- c(1,1)
  # 
  # Create_Tension_matrix(SmallEdgeNode, Smallzvect, Smallkvect, Smalldvect)
  
  #Creat the adjacency matrix
  Adj <- (t(EdgeNode) %*% diag(1, nrow = nrow(EdgeNode)) %*% EdgeNode) / (t(EdgeNode) %*% diag(1, nrow = nrow(EdgeNode)) %*% EdgeNode) 
  diag(Adj) <-0
  Adj[!is.finite(Adj)] <-0
  Adj
  
  Zmat <- Adj*zvect #The adjacency matrix weight by node height
  dZmat <- t(Zmat)- Zmat #The difference in height between adjacent nodes
  #Create the absolute K (spring stiffness) a and distance matrices
  kmat <- t(EdgeNode) %*% diag(kvect, nrow = length(kvect)) %*% EdgeNode %>% abs 
  Dmat <- t(EdgeNode) %*% diag(dvect, nrow = length(kvect)) %*% EdgeNode %>% abs
  
  Hmat <- sqrt(dZmat^2 + Dmat^2)
  
  Ften_mat <- kmat*(Hmat-Dmat)*dZmat/Hmat
  Ften_mat[!is.finite(Ften_mat)] <-0
  
  return (Ften_mat)
}
