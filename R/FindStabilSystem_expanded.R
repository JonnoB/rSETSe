#' Find stabil system expanded
#' 
#' This function uses the SETS embedding to find the equilibrium position of a network or a bi-connected component.
#' It produces the node status for every iteration of the process. useful for analysis and finding good starting parameters
#' However uses a lot more memory.
#' 
#' @param node_status A data frame The current dynamics and forces experienced by the node a data frame.
#' @param ten_mat A data frame The current dynamics and forces experienced by the node a data frame.
#' @param non_empty_matrix A numeric matrix. contains the index of the non-empty cells in the adjacency matrix. see details.
#' @param kvect A numeric vector of the spring stiffnesses
#' @param dvect A numeric vector of the initial distances between the nodes
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param tstep A numeric value. The time step, measured in seconds, that will be used to calculate the new dynamic state
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param coef_drag A numeric value. Used to set a multiplier on the friction value. Generally leave this alone..s.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @export
#' 
#' @details The non_empty matrixhe row column position absolute index and transpose index of the edges in the matrix
#' This means vectors can be used for most operation greatly reducing the amount of memory required and 
#' providing a modest speed increase. The non_empty_matrix is propduced by the 'Prepare_data_for_find_network_balance' function.
#'
# Strips out all pre processing to make it as efficient and simple as possible

#caoacity, edge_name and flow are no longer used. If the preprocessing is all done in prep then distance can also be removed

FindStabilSystem_expanded <- function(node_status, ten_mat, non_empty_matrix, kvect, dvect, mass,
                             tstep, max_iter = 1000, coef_drag = 1, 
                             tol = 1e-10, sparse = FALSE){
  #Runs the physics model to find the convergence of the system.
  
  #vectors are used throughout instead of a single matrix as it turns out they are faster due to less indexing and use much less RAM.
  
  #friction_stop fricton is a stopping condition. defualts to FALSE. 
  NodeList <- node_status[,-1]
  force <- NodeList[,1]
  elevation <-NodeList[,2]
  net_tension <-NodeList[,3]
  velocity <- NodeList[,4]
  friction <- NodeList[,5]
  static_force <-NodeList[,6]
  net_force <- NodeList[,7]
  acceleration <- NodeList[,8]
  
  #gets the dimensions of the matrix for bare bones column sum
  
  non_empty_vect <- non_empty_matrix[,1]
  non_empty_t_vect <- non_empty_matrix[,2]
  non_empty_index_vect <- non_empty_matrix[,3]
  #This dataframe is one of the final outputs of the function, it is premade for memory allocation
  total_nodes <- nrow(NodeList)
  results <- matrix(data = NA, nrow = total_nodes*(max_iter), ncol = ncol(NodeList))
  colnames(results) <- colnames(NodeList)
  
  one_vect <- rep(1, nrow(NodeList))
  
  Iter <- 1
  system_stable <- FALSE
  
  while((Iter <= max_iter) & !system_stable ){
    
    #calculate the system dynamics. Either sparse or dense mode
    #sparse or dense mode chosen by user on basis of network size and of course sparsity
    
    #The code is not put in sub-functions as this creates memory management problems and half the time
    #the program runs can be spent auto calling gc(). This reduces the copying of data...I think
    #It overwirtes the preious values but doesn't create anything else
    #NodeList2 <- NodeList
    
    #####
    #create the tension matrix
    #####
    #dz is the change in eleveation
    dzvect <- elevation[non_empty_t_vect] - elevation[non_empty_vect] #The difference in height between adjacent nodes 
    
    #the hypotenuse of the spring distance triangle
    Hvect <- sqrt(dzvect^2 + dvect^2)
    
    #the tension vector. the dZvect/Hvect is the vertical component of the tension
    ten_mat[non_empty_index_vect] <- kvect*(Hvect-dvect)*dzvect/Hvect
    
    #The remaining dynamics are calculated here
    
    elevation <- velocity*tstep +0.5*acceleration*tstep*tstep + elevation #Distance/elevation s1 = ut+0.5at^2+s0
    velocity <- velocity + acceleration*tstep #velocity v1 = v0 +at
    static_force <- force + net_tension #static force 
    
    if(sparse){
      #This uses the matrix row aggregation functions which can be used on sparse matrices. This is faster and much more memory
      #efficient for large matrices
      net_tension <- Matrix::rowSums(ten_mat) #tension
    }else{
      #This uses the standard dense matrices, this is faster for smaller matrices.
      net_tension <- ten_mat %*% one_vect  #.rowSums(ten_mat, m = m[1], n = m[1]) #tension
    }
    friction <- coef_drag*velocity #friction of an object in a viscous fluid under laminar flow
    net_force <- static_force - friction #net force
    acceleration <- net_force/mass #acceleration

    #add in the values into the holding matrix, this will slow down the algo
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),1] <- force
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),2] <- elevation
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),3] <- net_tension
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),4] <- velocity
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),5] <- friction
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),6] <- static_force 
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),7] <- net_force
    results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),8] <- acceleration
    

    
    
    #check if system is stableising the acceleration and max acceleration
    
    if(!is.finite(sum(abs(results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),6])))| 
       !is.finite(sum(abs(0.5*mass*results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),4]/tstep)))){ #if there are infinte values terminate early
      system_stable <- TRUE
      message("System has diverged. Process terminating")
    } else{
      system_stable <- (sum(abs(results[((Iter-1)*total_nodes+1):((Iter)*total_nodes),6])) < tol)
    }

    Iter <- Iter + 1 # add next iter
    
  }
  
  #Early termination causes NA values. These are removed by the below code
  #
  
  Iter_vect <- rep(1:(Iter-1), each = total_nodes)# The vector cannot be directly generated in mutate I don't know why.
  #Iter vect needs the iteration number to be repeated x times then the next iteration repeated x times. This is
  #different to the node name repetition
  results <- as_tibble(results) %>%
    mutate(node = node_status[rep(1:total_nodes, times = max_iter),"node"], #the nodename sequence needs to be repeated x times
           Iter = Iter_vect,
           t = tstep*Iter) %>%
    filter(complete.cases(.))
  

  Out <-bind_rows(node_status,
                  results )   #list(as_tibble(network_dynamics), )
  #names(Out) <- c("network_dynamics", "node_status")
  
  return(Out)
  
}