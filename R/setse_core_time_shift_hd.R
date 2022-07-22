#' setse Core with time shift for high dimensional feature networks
#'
#' Internal function. This is a variant of the SETse core algorithm. It runs the setse model to find the equilibrium position of the network,
#' it changes the time step if the algo is in the noisy zone. The function is called by auto_setse
#'
#' @param node_embeddings A data frame The current dynamics and forces experienced by the node a data frame.
#' @param ten_mat A data frame The current dynamics and forces experienced by the node a data frame.
#' @param non_empty_matrix A numeric matrix. contains the index of the non-empty cells in the adjacency matrix. see details.
#' @param kvect A numeric vector of the spring stiffnesses
#' @param dvect A numeric vector of the initial distances between the nodes
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param tstep A numeric value. The time step, measured in seconds, that will be used to calculate the new dynamic state
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param coef_drag A numeric value. Used to set a multiplier on the friction value. This is usualy determined by setse_auto
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample.
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is the system absolute mean force.
#' @param tstep_change a numeric scaler. A value between 0 and one, the fraction the new timestep will be relative to the previous one
#' @param dynamic_reset Logical. Whether the dynamic portion of the emebeddings is reset to zero when the timestep is changed
#' this can stop the momentum of the nodes forcing a divergence, but also can slow down the process. default is TRUE.
#' @param timeshift Logical. Whether the algorithm will adjust the timestep when in the noisy convergence zone.
#' This variable is generally handled by other SETSe functions
#' @param noisy_termination Logical. Whether the algorithm will terminate if the convergence enters the noisy zone. This
#' value takes precedence over timeshift and can be chosen by the user.
#' @param verbose_reporting Logical. Whether convergence details are reporting each sample.
#' This is useful for large networks which can take a long time to converge, but for smaller ones can be turned off.
#'
#' @details
#' This function is usally run inside a more easy to use function such as The setse function, SETse_bicomp or setse_auto. These
#' wrapper functions make the application of the SETse algorithm more straight foreword. However, this function is included
#' for completeness and to allow ground up experiments
#'
#' The non_empty matrix contains the row, column position and absolute index and transpose index of the edges in the matrix
#' This means vectors can be used for most operations greatly reducing the amount of memory required and
#' providing a modest speed increase. The non_empty_matrix is produced by the 'Prepare_data_for_find_network_balance' function.
#'
#' @return A list of three dataframes
#' \enumerate{
#'   \item The network dynamics describing several key figures of the network during the convergence process, this includes the static_force
#'   \item The node embeddings. Includes all data on the nodes the forces exerted on them position and dynamics at simulation termination
#'   \item A data frame giving the time taken for the simulation as well as the number of nodes and edges. Node and edge data is given
#'   as this may differ from the total number of nodes and edges in the network depending on the method used for convergnence.
#'   For example if setse_bicomp is used then some simulations may contain as little as two nodes and 1 edge
#' }
#'
#' @noRd
# Strips out all pre processing to make it as efficient and simple as possible

setse_core_time_shift_hd <- function(node_embeddings,
                       ten_mat,
                       non_empty_matrix,
                       kvect,
                       dvect,
                       mass,
                       tstep,
                       max_iter = 1000,
                       coef_drag = 1,
                       tol = 2e-3,
                       sparse = FALSE,
                       sample = 1,
                       static_limit = NULL,
                       tstep_change = 0.5,
                       dynamic_reset = TRUE,
                       timeshift,
                       noisy_termination,
                       verbose = TRUE){
  #Runs the physics model to find the convergence of the system.

  #vectors are used throughout instead of a single matrix as it turns out they are faster due to less indexing and use much less RAM.

  #These have to be matrices if there is a mutli-variable option
  NodeList <- node_embeddings[,-1]
  force <- NodeList %>% dplyr::select(dplyr::starts_with("force_")) %>% as.matrix()
  elevation <- NodeList %>% dplyr::select(dplyr::starts_with("elevation_")) %>% as.matrix()
  net_tension <-NodeList %>% dplyr::select(dplyr::starts_with("net_tension_")) %>% as.matrix()
  velocity <- NodeList %>% dplyr::select(dplyr::starts_with("velocity_")) %>% as.matrix()
  friction <- NodeList %>% dplyr::select(dplyr::starts_with("friction_")) %>% as.matrix()
  static_force <-NodeList %>% dplyr::select(dplyr::starts_with("static_force_")) %>% as.matrix()
  net_force <- NodeList %>% dplyr::select(dplyr::starts_with("net_force_")) %>% as.matrix()
  acceleration <- NodeList %>% dplyr::select(dplyr::starts_with("acceleration_")) %>% as.matrix()

  if(sparse){
    ten_mat <- methods::as(ten_mat, "dgTMatrix") # this is done as Dgt allows direct insertion of tension without indexing. It
    #is much faster than the standard format which does require indexing. This is despite dgt being slower to sum the columns
  }

  #The default value for the static limit if null is the sum of the absolute force.
  #This value is chosen because with good parameterization the static force never exceeds the starting amount.
  if(is.null(static_limit)){
    static_limit <- sum(abs(force))
  }

  #gets the dimensions of the matrix for bare bones column sum

  non_empty_vect <- non_empty_matrix[,1]
  non_empty_t_vect <- non_empty_matrix[,2]
  non_empty_index_vect <- non_empty_matrix[,3]
  non_empty_t_index_vect <- non_empty_matrix[,4]

  #This dataframe is one of the final outputs of the function, it is premade for memory allocation.
  #Although it would be faster to use vectors, the matrix is used only a fraction of the iteration, so has
  #very little impact on speed.
  network_dynamics <- matrix(data = NA, nrow = max_iter/sample, ncol = 6) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    purrr::set_names(c("Iter","t", "static_force", "kinetic_force", "potential_energy", "kinetic_energy")) %>%
    as.matrix()

  network_dynamics_initial_value <- network_dynamics[1:2, ,drop = FALSE]
  network_dynamics_initial_value[1,] <- c(0, 0, sum(abs(force)), 0,0,0)
  network_dynamics_initial_value <- network_dynamics_initial_value[1 , ,drop = FALSE]

  one_vect <- rep(1, nrow(NodeList))

  Iter <- 1
  current_time <- 0
  is_noisy <- FALSE
  system_stable <- FALSE

  #get the time the algo starts
  start_time <- Sys.time()

  #This should never need to be triggered, but is here as a fail safe
  if(timeshift & noisy_termination){
    print("noisy termination and timeshift conflict. Noisy termination takes precedence setting timeshift to FALSE")
    timeshift = FALSE
  }

  if(sparse){
  result_list <- core_while_loop_sparse_cpp(max_iter, sample, network_dynamics, elevation, non_empty_t_vect-1, non_empty_vect-1,
                                            dvect, velocity, acceleration, static_force,force, net_tension, tstep, ten_mat, kvect,
                                            mass, dynamic_reset, tstep_change, net_force, coef_drag, static_limit, tol, timeshift,
                                            noisy_termination, verbose_reporting = verbose)

  } else {

  result_list <- core_while_loop_dense_cpp(max_iter, sample, network_dynamics,
                                          elevation, non_empty_t_vect-1, non_empty_vect-1, dvect,
                                          velocity, acceleration, static_force, force, net_tension,  tstep,
                                          as.matrix(ten_mat), kvect,  mass,
                                          dynamic_reset,  tstep_change, net_force,  coef_drag,  static_limit,
                                          tol, non_empty_index_vect-1, timeshift,
                                          noisy_termination, verbose_reporting = verbose)

  }


  stop_time <- Sys.time()

  time_taken_df <- tibble::tibble(time_diff = stop_time-start_time,
                          nodes = nrow(node_embeddings),
                          edges = length(kvect))

  # The matrices produced by my  rcpp armadillo code do not have column headers, these are added here
  colnames(result_list$network_dynamics) <- colnames(network_dynamics_initial_value)
  #Early termination causes NA values. These are removed by the below code
  #

  tstep = result_list$tstep
  Iter = result_list$Iter

 node_embeddings_names <- colnames(node_embeddings)
 node_embeddings <- dplyr::bind_cols(node_embeddings[,"node",drop=FALSE] ,
                         #combine all the values for each dimension into a single tibble
                         list(force,
                              result_list$elevation,
                              result_list$net_tension,
                              result_list$velocity,
                              result_list$friction,
                              result_list$static_force,
                              result_list$net_force,
                              result_list$acceleration) %>%
                           do.call(what = cbind, args = .) %>%
                           dplyr::as_tibble() %>%
                           dplyr::mutate(t = tstep*(Iter-1),
                                         Iter = Iter-1)

  )

 colnames(node_embeddings) <- node_embeddings_names


  Out <- list(network_dynamics = dplyr::bind_rows(as.data.frame(network_dynamics_initial_value), as.data.frame(result_list$network_dynamics)),
              node_embeddings = node_embeddings ,  #1 needs to be subtracted from the total as the final thing
              #in the loop is to add 1 to the iteration
              time_taken = time_taken_df #This is a diff time object!
              )

  return(Out)

}
