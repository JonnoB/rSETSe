#' Calculate system dynamics sparse
#' 
#' This is a helper function that calculates the dynamics at time t of all the nodes in the system.
#' It is seldom called on it's own but is called by other functions. This function comes in a version for sparse and dense matrices.
#' 
#' The output of the function is a matrix that is the next iteration of the NodeStatus dataframe the function recieves.
#' 
#' @param NodeStatus A numeric matrix frame The current dynamics and forces experienced by the nodes. The node ID is not included
#' This means the order of the matrix has to be fixed.
#' @param ten_mat A numeric matrix. This matrix must be a zero matrix or an adjacency matrix of dimensions nxn where n is the number of nodes
#' @param damp_mat Identical to the ten_mat. This is supposed to help reduce memory, But I don't think it does
#' @param kvect A numeric vector of the spring stiffnesses
#' @param dvect A numeric vector of the initial distances between the nodes
#' @param mat_size a numeric scaler. The dimensions of the square matrix
#' @param tstep A numeric value. The time step, measured in seconds, that will be used to calculate the new dynamic state
#' @param non_empty_matrix A numeric matrix. This contains the node indexes in the adjacency matrix it has 4 columns
#'  the row, the column, the absolute index the transpose of the absolute index.
#' @param frctmultiplier A numeric value. Used to set a multiplier on the friction value. Generally leave this alone.
#'
#' @export

Calc_System_Dynamics_sparse<- function(NodeStatus, ten_mat, damp_mat, kvect, dvect, mat_size, tstep = 1, non_empty_matrix, frctmultiplier = 1){
  
  #I am slightly unsure if the damping and tension are being updated in the right place. This isn't the end of the world
  #If they are slightly behind but It may make the solution slightly more stable which would give faster convergence
  
  # #Create the damping matrix
  damp_mat[non_empty_matrix[,3]]<- 2*sqrt(kvect*NodeStatus[non_empty_matrix[,1],3])*NodeStatus[non_empty_matrix[,1],5]

  NodeStatus2 <- NodeStatus
  NodeStatus2[,4] <- Matrix::rowSums(Create_Tension_matrix3(ten_mat, 
                                                     zvect = NodeStatus[non_empty_matrix[,1],2], 
                                                     zvect_t = NodeStatus[non_empty_matrix[,2],2], 
                                                     dvect, kvect,
                                                     non_empty_index = non_empty_matrix[,3])) #tension
  # NodeStatus2[,4] <- .rowSums(Create_Tension_matrix3(ten_mat, 
  #                                                    zvect = NodeStatus[non_empty_matrix[,1],2], 
  #                                                    zvect_t = NodeStatus[non_empty_matrix[,2],2], 
  #                                                    dvect, kvect,
  #                                                    non_empty_index = non_empty_matrix[,3]), 
  #                             m = mat_size, n = mat_size) #tension
  NodeStatus2[,2] <- Distance(NodeStatus[,2], NodeStatus[,5], NodeStatus[,8], t0 = NodeStatus[,10], t1 = NodeStatus[,10] + tstep) #distance
  NodeStatus2[,5] <- velocity(NodeStatus[,5], NodeStatus[,8], NodeStatus[,10], NodeStatus[,10] + tstep) #velocity

  #NodeStatus2[,6] <- .rowMeans(damp_mat, m = mat_size, n = mat_size)*frctmultiplier #friction
  NodeStatus2[,6] <- Matrix::rowMeans(damp_mat)*frctmultiplier #friction
  NodeStatus2[,7] <- NodeStatus2[,1] + NodeStatus2[,4] - NodeStatus2[,6] #Netforce
  NodeStatus2[,8] <- NodeStatus2[,7]/NodeStatus2[,3] #acceleration
  NodeStatus2[,9] <- (NodeStatus2[,8]-NodeStatus[,8])/tstep #delta acceleration...this can be removed
  NodeStatus2[,10] <- NodeStatus[,10]+tstep #time
  #This old school method is much faster than using mutate
  # NodeStatus2 <- NodeStatus
  # NodeStatus2[,2] <- Distance(NodeStatus[,2], NodeStatus[,5], NodeStatus[,8], t0 = NodeStatus[,10], t1 = NodeStatus[,10] + tstep)
  # NodeStatus2[,5] <- velocity(NodeStatus[,5], NodeStatus[,8], NodeStatus[,10], NodeStatus[,10] + tstep)
  # NodeStatus2[,4] <- .rowSums(Create_Tension_matrix3(ten_mat, 
  #                                                    zvect = NodeStatus[non_empty_matrix[,1],2], 
  #                                                    zvect_t = NodeStatus[non_empty_matrix[,2],2], 
  #                                                    dvect, kvect,
  #                                                    non_empty_index = non_empty_matrix[,3]), 
  #                             m = mat_size, n = mat_size)
  # NodeStatus2[,6] <- .rowMeans(damp_mat, m = mat_size, n = mat_size)*frctmultiplier
  # NodeStatus2[,7] <- NodeStatus2[,1] + NodeStatus2[,4] - NodeStatus2[,6]
  # NodeStatus2[,8] <- NodeStatus2[,7]/NodeStatus2[,3]
  # NodeStatus2[,9] <- (NodeStatus2[,8]-NodeStatus[,8])/tstep
  # NodeStatus2[,10] <- NodeStatus[,10]+tstep
  
    # Tension_vect <- rowSums(Create_Tension_matrix2(Adjmat, kmat, dmat, NodeStatus$z))
  # Friction_vect <- rowMeans(Calc_Damping_matrix2(Adjmat, NodeStatus$velocity, kmat, NodeStatus$mass))
  # 
  # #This old school method is much faster than using mutate
  # NodeStatus2 <- NodeStatus
  # NodeStatus2$z <- Distance(NodeStatus$z, NodeStatus$velocity, NodeStatus$acceleration, t0 = NodeStatus$t, t1 = NodeStatus$t + tstep)
  # NodeStatus2$velocity <- velocity(NodeStatus$velocity, NodeStatus$acceleration, NodeStatus$t, NodeStatus$t + tstep)
  # NodeStatus2$NetTension <- Tension_vect
  # NodeStatus2$friction <- Friction_vect*frctmultiplier
  # NodeStatus2$NetForce <- NodeStatus2$force + NodeStatus2$NetTension - NodeStatus2$friction
  # NodeStatus2$acceleration <- NodeStatus2$NetForce/NodeStatus2$mass
  # NodeStatus2$Delta_acceleration <- (NodeStatus2$acceleration-NodeStatus$acceleration)/tstep
  # NodeStatus2$t <- NodeStatus2$t+tstep

  
  #paste("acceleration", sum(abs(NodeStatus2$acceleration)), "velocity", sum(abs(NodeStatus2$velocity)))
  
  return(NodeStatus2)
  
}