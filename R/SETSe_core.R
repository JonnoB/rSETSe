#' SETSe Core
#' 
#' This is the SETse core algorithm. It runs the SETSe model to find the equilibrium position of the network
#' 
#' This function is usally run inside a more easy to use function such as The SETSe function, SETse_bicomp or auto_SETSe. These
#' wrapper functions make the application of the SETse algorithm more straight foreword. However, this function is included
#' for completeness and to allow ground up experiments
#' 
#' @param node_embeddings A data frame The current dynamics and forces experienced by the node a data frame.
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
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @export
#' 
#' @details The non_empty matrix contains the row, column position and absolute index and transpose index of the edges in the matrix
#' This means vectors can be used for most operations greatly reducing the amount of memory required and 
#' providing a modest speed increase. The non_empty_matrix is produced by the 'Prepare_data_for_find_network_balance' function.
#'
# Strips out all pre processing to make it as efficient and simple as possible

SETSe_core <- function(node_embeddings, ten_mat, non_empty_matrix, kvect, dvect, mass,
                                  tstep, max_iter = 1000, coef_drag = 1, 
                                  tol = 2e-3, sparse = FALSE,
                             sample = 1){
  #Runs the physics model to find the convergence of the system.
  
  #vectors are used throughout instead of a single matrix as it turns out they are faster due to less indexing and use much less RAM.
  
  #friction_stop fricton is a stopping condition. defualts to FALSE. 
  NodeList <- node_embeddings[,-1]
  force <- NodeList[,1]
  elevation <-NodeList[,2]
  net_tension <-NodeList[,3]
  velocity <- NodeList[,4]
  friction <- NodeList[,5]
  static_force <-NodeList[,6]
  net_force <- NodeList[,7]
  acceleration <- NodeList[,8]
  
  #The static limit is 10 times the static force
  #Sometimes numbers can explode then converge, but whatever I don't care about them
  static_limit <- sum(abs(force))*10
  
  #gets the dimensions of the matrix for bare bones column sum
  
  non_empty_vect <- non_empty_matrix[,1]
  non_empty_t_vect <- non_empty_matrix[,2]
  non_empty_index_vect <- non_empty_matrix[,3]
  #This dataframe is one of the final outputs of the function, it is premade for memory allocation
  network_dynamics <- matrix(data = NA, nrow = max_iter/sample, ncol = 6) %>%
    as_tibble() %>%
    set_names(c("Iter","t", "static_force", "kinetic_force", "potential_energy", "kinetic_energy")) %>%
    as.matrix()
  
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
    # NodeList[,9] <- NodeList[,9] + tstep #current time #This may not be neccessary but doesn't really hurt
    
    if((Iter %% sample)==0){
      network_dynamics[Iter/sample,]<-  c(Iter, #Iteration
                                          Iter*tstep, #time in seconds
                                          sum(abs(static_force)),  #static force. The force exerted on the node
                                          sum(abs(0.5*mass*velocity/tstep)), #kinetic_force 
                                          sum( 0.5*kvect*(Hvect-dvect)^2),     #spring potential_energy
                                          sum(0.5*mass*velocity^2)    #kinetic_energy
      ) 
      
    
      #check if system is stable using static force
      #If static force is not finite or exceeds the static limit then the process terminates early
      if(!is.finite(network_dynamics[Iter/sample,3])| network_dynamics[Iter/sample,3]>static_limit){ #if there are infinte values terminate early
        system_stable <- TRUE
       # print(network_dynamics[Iter/sample,3])
      } else{
        system_stable <- (network_dynamics[Iter/sample,3] < tol)
      }
      
    }
    
    Iter <- Iter + 1 # add next iter
    
  }
  
  #Early termination causes NA values. These are removed by the below code
  #
  network_dynamics <- as_tibble(network_dynamics) %>%
    filter(complete.cases(.))
  #combine all the vectors together again into a tibble
  Out <- list(network_dynamics = as_tibble(network_dynamics), 
              node_embeddings = bind_cols(node_embeddings[,"node",drop=FALSE] , 
                                      tibble(  force = force,
                                               elevation = as.vector(elevation),
                                               net_tension = as.vector(net_tension),
                                               velocity = as.vector(velocity),
                                               friction = as.vector(friction),
                                               static_force = as.vector(static_force),
                                               net_force = as.vector(net_force),
                                               acceleration = as.vector(acceleration),
                                               t = tstep*(Iter-1),
                                               Iter = Iter-1))) #1 needs to be subtracted from the total as the final thing
  #in the loop is to add 1 to the iteration
  
  return(Out)
  
}