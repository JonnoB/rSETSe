#' Create stabilised blocks
#' 
#' An internal function. This function is called by SETSe_bicomp and performs auto_SETSe on all the bi-connected components
#' of the network. This function is rarely called directly.
#' 
#' @param g An igraph object
#' @param OriginBlock An Igraph object. This is created by Create_balanced_blocks function
#' @param OriginBlock_number An integer. this is the origin block chosen from the
#' create_stable_blocks function. Usually this will be the largest block.
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param k A character string. name of the spring constant variable
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param sample Integer. The dynamics will be stored only if the iteration number is a multiple of the sample. 
#'  This can greatly reduce the size of the results file for large numbers of iterations. Must be a multiple of the max_iter
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is twice the system absolute mean force.
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' @param hyper_iters integer. The hyper parameter that determines the number of iterations allowed to find an acceptable convergence value.
#' @param drag_min integer. A power of ten. The lowest drag value to be used in the search
#' @param drag_max integer. A power of ten. if the drag exceeds this value the tstep is reduced
#' @param tstep_change numeric. A value between 0 and 1 that determines how much the time step will be reduced by default value is 0.5
#' @param bigraph A list. the list of biconnected components produced by the biconnected_components function.
#'  This function take a non trivial amount of time on large graphs so this pass through minimises the function being called.
#' @param hyper_tol numeric. The convergence tolerance when trying to find the minimum value
#' @param hyper_max integer. The maximum number of iterations that the setse will go through whilst searching for the minimum.
#' @param balanced_blocks A list 
#' @param noisey_termination Stop the process if the static force does not monotonically decrease.
#' 
#' @details This function isn't really supposed to be used apart from as a sub-routine of the SETSe biconnected component method.
#' 
#' @seealso \code{\link{SETSe_bicomp}} \code{\link{SETSe}}
#' @return A dataframe with the height embeddings of the network
#' 
#' 
#' @noRd
Create_stabilised_blocks <- function(g, 
                                     OriginBlock, 
                                     OriginBlock_number, 
                                     force ="force",
                                     edge_name = "edge_name",
                                     k = "k",
                                     tstep=0.1, 
                                     tol = 1e-10, 
                                     distance, 
                                     max_iter, 
                                     mass,
                                     sparse,
                                     sample = 100,
                                     static_limit = NULL,
                                     hyper_iters = 100,
                                     hyper_tol  = 0.01,
                                     hyper_max = 30000,
                                     drag_min = drag_min,
                                     drag_max = drag_max,
                                     tstep_change = tstep_change,
                                     verbose = FALSE,
                                     bigraph = bigraph,
                                     balanced_blocks = balanced_blocks,
                                     noisey_termination = TRUE){

  #remove the Origin block so it doesn't have to be calculated again
  BlockNumbers <-(1:length(balanced_blocks))[-OriginBlock_number]
  
  #total in network
  total_force <- sum(abs(igraph::get.vertex.attribute(g, force)))
  
  #Get timing for rest of process
  start_time_all_other_blocks <- Sys.time()
  
  StabilModels <- BlockNumbers %>% 
    purrr::map(~{

      #sub tol can be extremely small. if it is then I say that it is zero and effectively the nodes are left unconverged
      #again, it is worth considering the magnitude of force in the network
      sub_tol <- tol*sum(abs(igraph::get.vertex.attribute(balanced_blocks[[.x]], force)))/total_force
      sub_tol_larger_than_limit <- sub_tol> .Machine$double.eps^0.5
      sub_tol <- ifelse(sub_tol_larger_than_limit, sub_tol, tol)
      
      if(!is.null(static_limit)){
        
      sub_static_limit <- static_limit*sum(abs(igraph::get.vertex.attribute(balanced_blocks[[.x]], force)))/total_force

      } else {
        
        sub_static_limit <- sum(abs(igraph::vertex_attr(g, force))) #NULL
        
      }
      
      #makes sure the static limit is not smaller than the machine precision. if it is bad things happen
      sub_static_larger_than_limit <- sub_static_limit > .Machine$double.eps^0.5
      sub_static_limit <- ifelse(sub_static_larger_than_limit, sub_static_limit, static_limit)
      
      #do special case solution for two nodes only
      if(igraph::ecount(balanced_blocks[[.x]])==1){
        
        Prep <- SETSe_data_prep(g = balanced_blocks[[.x]], 
                                force = force, 
                                distance = distance, 
                                mass = ifelse(is.null(mass), mass_adjuster(balanced_blocks[[.x]], 
                                                                           force = "force", resolution_limit = TRUE), mass), 
                                k = k,
                                edge_name = edge_name,
                                sparse = sparse)
        
        Out <- two_node_solution(g, Prep = Prep, auto_setse_mode = TRUE)
        
        #S if the subtolerance and the sub static values are larger than the sqrt of the machine eps
        #Then solve using the auto-SETSe method.
      }  else if(sub_tol_larger_than_limit  & sub_tol_larger_than_limit){
        
        print(paste("Block", .x, "of", max(BlockNumbers),  "has more than two nodes. Running auto-setse"))
        
        start_time_block <- Sys.time()
        
        Out <- auto_SETSe(balanced_blocks[[.x]],
                          force = force,
                          distance = distance, 
                          edge_name = edge_name,
                          k = k,
                          tstep = tstep, 
                          tol = sub_tol, #the force has to be scaled to the component                           
                          max_iter =  max_iter, 
                          mass =  ifelse(is.null(mass), mass_adjuster(balanced_blocks[[.x]], 
                                                                      force = "force", resolution_limit = TRUE), mass), 
                          sparse = sparse,
                          sample = sample,
                          static_limit = sub_static_limit,
                          hyper_iters = hyper_iters,
                          hyper_tol = hyper_tol,
                          hyper_max = hyper_max,
                          drag_max = drag_max,
                          tstep_change = tstep_change,
                          verbose = verbose,
                          include_edges = FALSE,
                          noisey_termination = noisey_termination
        )
        
        #Print the details of the completed block
        #This is not done for biconns of size 2 they are so fast and numerous
        print(paste("Block", .x, "of", max(BlockNumbers), 
                    "complete. bicomp has", 
                    igraph::vcount(balanced_blocks[[.x]]), "nodes.",
                    " Time taken",
                    round(as.numeric( difftime(Sys.time(), start_time_block, units = "mins")), 1),
                    " minutes"))
        #Otherwise the values are smaller than the machine tolerance and should be left unconverged.
        #For very small values, uncoverged is approximately equal to the converged value
      } else {
        
        print(paste("Block", .x, "of", max(BlockNumbers),  
                    "has more than two nodes but approximately 0 force, no convergence necessary, continuing to next block"))
        
        Prep <- SETSe_data_prep(g = balanced_blocks[[.x]], 
                                force = force, 
                                distance = distance, 
                                mass = ifelse(is.null(mass), mass_adjuster(balanced_blocks[[.x]], 
                                                                           force = "force", resolution_limit = TRUE), mass), 
                                k = k,
                                edge_name = edge_name,
                                sparse = sparse)
        
        
        Out <- list(network_dynamics = tibble::tibble(t = 0, 
                                              Iter = 0,
                                              static_force = 0, 
                                              kinetic_force = 0), 
                    node_embeddings = Prep$node_embeddings,
                    #NA value is included as the timing is not applicable.
                    #It also indicates that the forces were so small as to not rquire convergence
                    time_taken = tibble::tibble(time_diff = NA, nodes = igraph::vcount(balanced_blocks[[.x]]), 
                                        edges =  igraph::ecount(balanced_blocks[[.x]])),
                    memory_df = tibble::tibble(iteration = 1,
                                       error = NA,
                                       perc_change = NA,
                                       log_ratio = NA,
                                       common_drag_iter = NA,
                                       tstep = NA,
                                       direction = 1,
                                       target_area = NA,
                                       res_stat = NA,
                                       upper = NA,
                                       lower = NA,
                                       best_log_ratio =NA,
                                       stable = NA))
        
        
      }

      return(Out)
      
    })
  
  
  if(verbose){print(paste("All Other blocks completed, time taken", 
                          round(as.numeric( difftime(Sys.time(), start_time_all_other_blocks, units = "mins")), 1), 
                          "minutes."))}

  if(verbose){print("Re-assembling network")}
  
  start_reassembly <- Sys.time()
  #extract the articulation nodes
  ArticulationVect <- igraph::get.vertex.attribute(g, "name", bigraph$articulation_points)
  
  #place all nodes relative to the origin
  relative_blocks <- 1:length(StabilModels) %>% 
    purrr::map_df(~{
      #print(.x) #It is a bit annoying and pointless now
      temp <- StabilModels[[.x]]$node_embeddings
      temp$Reference_ID <- .x
      
      return(temp)
      
    }) %>%
    dplyr::bind_rows(OriginBlock$node_embeddings %>% 
                dplyr::mutate(Reference_ID = 0)) %>%
    dplyr::mutate(Articulation_node = (node %in% ArticulationVect ))
  
  #get the network_dynamics dataframe for the total calculation
  network_dynamics <- 1:length(StabilModels) %>% 
    purrr::map_df(~{
      temp <- StabilModels[[.x]]$network_dynamics
      temp$component <-  BlockNumbers[.x]
      
      return(temp)
      
    }) %>%
    dplyr::bind_rows(OriginBlock$network_dynamics %>% dplyr::mutate(component = OriginBlock_number))  #%>%
    #It can also be useful to get the individual component values out.
    # group_by(Iter) %>%
    # summarise_all(sum) %>%
    # mutate(t = Iter*tstep) %>%
    # select(-component)
  
  
  #gets back the memory of the convergence process.
  #This is useful for knowing what logratio works for various network families
  memory_df <- 1:length(StabilModels) %>% 
    purrr::map_df(~{
      temp <- StabilModels[[.x]]$memory_df 
      temp$component <-  BlockNumbers[.x] #puts in the original block ordering so we can know which block was which
      return(temp)
    }) %>%
    dplyr::bind_rows(OriginBlock$memory_df %>% dplyr::mutate(component = OriginBlock_number)) 
  
  
  time_taken_df <- 1:length(StabilModels) %>% 
    purrr::map_df(~{
      temp <- StabilModels[[.x]]$time_taken
      temp$component <-  BlockNumbers[.x] #puts in the original block ordering so we can know which block was which
      return(temp)
    }) %>%
    dplyr::bind_rows( OriginBlock$time_taken %>%
                 dplyr::mutate(component = OriginBlock_number)) 
  
  #The biconnected components are converted to absolute values from relative ones
  node_embeddings <- fix_elevation_to_origin(relative_blocks, ArticulationVect) 
  
  print("Removing multiple articulation nodes")
  #this bind rows takes place as, there are vastly more non-articulation nodes than 
  #articulation nodes, this can mean and absolutely massive number of groups which is very slow to summarise
  #by aggregating only the necessary nodes it will be much faster
  node_embeddings <- dplyr::bind_rows(node_embeddings[!(node_embeddings$node %in% ArticulationVect),],
                               remove_articulation_duplicates(node_embeddings, ArticulationVect))  %>%
    #friction is different depending on the biconnected component so is meaningless in the overall analysis
    #Poissibly could use the drag for the major biconn
    dplyr::mutate(
      friction = NA,#coef_drag * velocity,
      static_force = force + net_tension,
      net_force = NA,#static_force - friction,
      acceleration = net_force/ifelse(is.null(mass), mass_adjuster(balanced_blocks[[OriginBlock_number]], 
                                                                   force = "force", resolution_limit = TRUE), mass),
      t = 1,
      t = tstep * Iter) %>%
    dplyr::arrange(node)
  
  if(verbose){print(paste("Re-assembly complete, time taken", 
                          round(as.numeric( difftime(Sys.time(), start_reassembly, units = "mins")), 1), 
                          "minutes."))}

  
  #combine the node_embeddings and the network_dynamics into a single list
  Out <- list(node_embeddings = node_embeddings, 
              network_dynamics = network_dynamics,
              memory_df = memory_df,
              time_taken = time_taken_df)
  
  return(Out)
  
}