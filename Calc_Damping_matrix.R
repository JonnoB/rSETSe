Calc_Damping_matrix <- function(EdgeNode, v_vect, kvect, mvect){
  
  #The v_vect, kvect and mvect all have to be ordered alphanumerically!
  #EdgeNode <- The edgenode matrix 
  #v_vect the velocity of each node
  #kvect the vector of spring stiffness for each edge
  #dvect, the vector of horizontal distance between the springs.
  
  #Creat the adjacency matrix
  Adj <- (t(EdgeNode) %*% diag(1, nrow = nrow(EdgeNode)) %*% EdgeNode) / (t(EdgeNode) %*% diag(1, nrow = nrow(EdgeNode)) %*% EdgeNode) 
  diag(Adj) <-0
  Adj[!is.finite(Adj)] <-0
  
  Vmat <- Adj*v_vect #The adjacency matrix weight by node height
  #Create the absolute K (spring stiffness) a and distance matrices
  kmat <- t(EdgeNode) %*% diag(kvect, nrow = length(kvect)) %*% EdgeNode %>% abs 
  
  Ften_mat<- 2*sqrt(kmat*mvect)*Vmat
  Ften_mat[!is.finite(Ften_mat)] <-0
  
  return((Ften_mat))
}