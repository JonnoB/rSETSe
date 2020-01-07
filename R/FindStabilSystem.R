#' Find stabil system
#' 
#' This function uses the SETS embedding to find the equilibrium position of a network or a bi-connected component
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
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' @export
#' 
#' @details The non_empty matrixhe row column position absolute index and transpose index of the edges in the matrix
#' This means vectors can be used for most operation greatly reducing the amount of memory required and 
#' providing a modest speed increase. The non_empty_matrix is propduced by the 'Prepare_data_for_find_network_balance' function.
#'
# Strips out all pre processing to make it as efficient and simple as possible

#caoacity, edge_name and flow are no longer used. If the preprocessing is all done in prep then distance can also be removed

FindStabilSystem <- function(node_status, ten_mat, non_empty_matrix, kvect, dvect, mass,
                              tstep, max_iter = 1000, coef_drag = 1, 
                              tol = 1e-10, sparse = FALSE, verbose = FALSE){
  #Runs the physics model to find the convergence of the system.
  
  #friction_stop fricton is a stopping condition. defualts to FALSE. 
  NodeList <- NodeList2 <-as.matrix(node_status[,-1])
  
  #gets the dimensions of the matrix for bare bones column sum
  m <- dim(NodeList)
  
  #This dataframe is one of the final outputs of the function, it is premade for memory allocation
  network_dynamics <- matrix(data = NA, nrow = max_iter, ncol = 6) %>%
    as_tibble() %>%
    set_names(c("Iter","t", "static_force", "kinetic_force", "potential_energy", "kinetic_energy")) %>%
    as.matrix()
  
  Iter <- 1
  system_stable <- FALSE
  
  while((Iter <= max_iter) & !system_stable ){
    
    #calculate the system dynamics. Either sparse or dense mode
    #sparse or dense mode chosen by user on basis of network size and of course sparsity
    
    #The code is not put in sub-functions as this creates memory management problems and half the time
    #the program runs can be spent auto calling gc(). This reduces the copying of data...I think
    #It overwirtes the preious values but doesn't create anything else
    NodeList2 <- NodeList
    #####
    #create the tension matrix
    #####
    #dz is the change in eleveation
    dzvect <- NodeList[non_empty_matrix[,2],2] - NodeList[non_empty_matrix[,1],2] #The difference in height between adjacent nodes 
    
    #the hypotenuse of the spring distance triangle
    Hvect <- sqrt(dzvect^2 + dvect^2)
    
    #the tension vector. the dZvect/Hvect is the vertical component of the tension
    ten_mat[non_empty_matrix[,3]] <- kvect*(Hvect-dvect)*dzvect/Hvect
    
    ####
    ## Create the Damping matrix
    ###
    #damp_mat[non_empty_matrix[,3]]<- 2*sqrt(kvect*NodeList[non_empty_matrix[,1],3])*NodeList[non_empty_matrix[,1],5]
    
    if(sparse){
      #This uses the matrix row aggregation functions which can be used on sparse matrices. This is faster and much more memory
      #efficient for large matrices
      NodeList2[,3] <- Matrix::rowSums(ten_mat) #tension
    }else{
      #This uses the standard dense matrices, this is faster for smaller matrices.
      NodeList2[,3] <- .rowSums(ten_mat, m = m[1], n = m[1]) #tension
    }
    #The remaining dynamics are calculated here
    
    #If these equations of motion work well then the distance and velocity equations can be removed
    NodeList2[,2] <- NodeList[,4]*tstep +0.5*NodeList[,8]*tstep*tstep + NodeList[,2] #Distance/elevation s1 = ut+0.5at^2+s0
    NodeList2[,4] <- NodeList[,4] + NodeList[,8]*tstep #velocity v1 = v0 +at
    NodeList2[,6] <- NodeList[,1] + NodeList[,3] #static force 
    NodeList2[,5] <- coef_drag*NodeList2[,4] #friction of an object in a viscous fluid under laminar flow
    NodeList2[,7] <- NodeList2[,6] - NodeList2[,5] #net force
    NodeList2[,8] <- NodeList2[,7]/mass #acceleration
    NodeList2[,9] <- NodeList[,9] + tstep #current time #This may not be neccessary but doesn't really hurt
    
    
    network_dynamics[Iter,]<-  c(Iter, #Iteration
                                 Iter*tstep, #time in seconds
                                 sum(abs(NodeList2[,6])),  #static force. The force exerted on the node
                                 sum(abs(0.5*mass*NodeList2[,4]/tstep)), #kinetic_force 
                                 sum( 0.5*kvect*(Hvect-dvect)^2),     #spring potential_energy
                                 sum(0.5*mass*NodeList2[,4]^2)    #kinetic_energy
    ) 
    
    NodeList <- NodeList2
    
    #check if system is stableusing the acceleration and max acceleration
    
    if(!is.finite(network_dynamics[Iter,3])| !is.finite(network_dynamics[Iter,4])){ #if there are infinte values terminate early
      system_stable <- TRUE
    } else{
      system_stable <- (network_dynamics[Iter,3] < tol)
    }
    
    #The system is too fast to print out iterations this should probably be removed
    # if(verbose){
    #   
    #   print(paste("Iteration", Iter,
    #               "net force", signif(network_dynamics[Iter,3], 3),
    #               "kinetic force", signif(network_dynamics[Iter,3], 3))
    #   ) # print result
    #   
    # }
    
    Iter <- Iter + 1 # add next iter
    
  }
  
  #Early termination causes NA values. These are removed by the below code
  #
  network_dynamics <- as_tibble(network_dynamics) %>%
    filter(complete.cases(.))
  
  Out <- list(as_tibble(network_dynamics), bind_cols(node_status[,"node",drop=FALSE] , as_tibble(NodeList)))
  names(Out) <- c("network_dynamics", "node_status")
  
  return(Out)
  
}