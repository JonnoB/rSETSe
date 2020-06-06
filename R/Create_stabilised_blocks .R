#' Create stabilised blocks
#' 
#' decomposes the network into bi-connected components using articulation points. This speeds up the convergence process
#' and reduces the chances of the SETS algorithm diverging on certain classes of network. This function is rarely called
#' directly and is usually called by SETSe_bicomp
#' 
#' @param g An igraph object
#' @param Origin block
#' @param OriginBlock_number An integer. this is the origin block chosen from the
#' create_stable_blocks function. Usually this will be the largest block.
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param sparse Logical. Whether or not the function should be run using sparse matrices. must match the actual matrix, this could prob be automated
#' @param static_limit Numeric. The maximum value the static force can reach before the algorithm terminates early. This
#' prevents calculation in a diverging system. The value should be set to some multiple greater than one of the force in the system.
#' If left blank the static limit is twice the system absolute mean force.
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' @param hyper_iters integer. The hyper parameter that determines the number of iterations allowed to find an acceptable convergence value.
#' @param step_size numeric. The hyper parameter that determines the log ratio search step size for auto convergence
#' @param hyper_tol numeric. The convergence tolerance when trying to find the minimum value
#' @param hyper_max integer. The maximum number of iterations that the setse will go through whilst searching for the minimum.
#' 
#' @seealso \code{\link{SETSe_bicomp}} \code{\link{SETSe}}
#' @return A dataframe with the height embeddings of the network
#' 
#' 
#' @export
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
                                     step_size = 0.1,
                                     verbose = FALSE,
                                     bigraph = bigraph,
                                     balanced_blocks = balanced_blocks){

  #remove the Origin block so it doesn't have to be calculated again
  BlockNumbers <-(1:length(balanced_blocks))[-OriginBlock_number]
  
  #total in network
  total_force <- sum(abs(get.vertex.attribute(g, force)))
  
  StabilModels <- BlockNumbers %>% 
    map(~{
      print(paste("Block", .x, "of", max(BlockNumbers), "bicomp has", vcount(balanced_blocks[[.x]]), "nodes."))
      
      #sub tol can be extremely small. if it is then I say that it is zero and effectively the nodes are left unconverged
      #again, it is worth considering the magnitude of force in the network
      sub_tol <- tol*sum(abs(get.vertex.attribute(balanced_blocks[[.x]], force)))/total_force
      sub_tol <- ifelse(sub_tol> .Machine$double.eps^0.5, sub_tol, tol)
      
      if(!is.null(static_limit)){
        
      sub_static_limit <- static_limit*sum(abs(get.vertex.attribute(balanced_blocks[[.x]], force)))/total_force
      #makes sure the static limit is not smaller than the machine precision. if it is bad things happen
      sub_static_limit <- ifelse(sub_static_limit> .Machine$double.eps^0.5, sub_static_limit, static_limit)
      } else {
        
        sub_static_limit <- NULL
        
      }
      
      #do special case solution for two nodes only
      if(ecount(balanced_blocks[[.x]])==1){
        
        Prep <- SETSe_data_prep(g = balanced_blocks[[.x]], 
                                force = force, 
                                distance = distance, 
                                mass = mass, 
                                k = k,
                                edge_name = edge_name,
                                sparse = sparse)
        
        Out <- two_node_solution(g, Prep = Prep, auto_setse_mode = TRUE)
        
        #Solves using the iterative method.
      } else {
        
        
        Out <- auto_SETSe(balanced_blocks[[.x]],
                          force = force,
                          distance = distance, 
                          edge_name = edge_name,
                          k = k,
                          tstep = tstep, 
                          tol = sub_tol, #the force has to be scaled to the component                           
                          max_iter =  max_iter, 
                          mass =  mass, 
                          sparse = sparse,
                          sample = sample,
                          static_limit = sub_static_limit,
                          hyper_iters = hyper_iters,
                          hyper_tol = hyper_tol,
                          hyper_max = hyper_max,
                          step_size = step_size,
                          verbose = verbose,
                          include_edges = FALSE
                          )
        
      }

      #print if the print requirement is on otherwise silent
      #if(!verbose){print(paste("Block" ,.x, "of", max(BlockNumbers) ,"termination", nrow(Out$network_dynamics) )) }
      
      return(Out)
      
    })

  #extract the articulation nodes
  ArticulationVect <- get.vertex.attribute(g, "name", bigraph$articulation_points)
  
  #place all nodes relative to the origin
  relative_blocks <- 1:length(StabilModels) %>% 
    map_df(~{
      #print(.x) #It is a bit annoying and pointless now
      StabilModels[[.x]]$node_embeddings %>%
        mutate(Reference_ID = .x)
      
    }) %>%
    bind_rows(OriginBlock$node_embeddings %>% 
                mutate(Reference_ID = 0)) %>%
    mutate(Articulation_node = (node %in% ArticulationVect ))
  
  #get the network_dynamics dataframe for the total calculation
  network_dynamics <- 1:length(StabilModels) %>% 
    map_df(~{
      StabilModels[[.x]]$network_dynamics %>%
        mutate(component = BlockNumbers[.x])
      
    }) %>%
    bind_rows(OriginBlock$network_dynamics %>% mutate(component = OriginBlock_number))  %>%
  #It can also be useful to get the individual component values out.
    group_by(Iter) %>%
    summarise_all(sum) %>%
    mutate(t = Iter/tstep) %>%
    select(-component)


  #gets back the memory of the convergence process.
  #This is useful for knowing what logratio works for various network families
  memory_df <- 1:length(StabilModels) %>% 
    map_df(~{
      StabilModels[[.x]]$memory_df %>%
        mutate(component = BlockNumbers[.x]) #puts in the original block ordering so we can know which block was which
      
    }) %>%
    bind_rows(OriginBlock$memory_df %>% mutate(component = OriginBlock_number)) 


time_taken_df <- 1:length(StabilModels) %>% 
    map_df(~{
      StabilModels[[.x]]$time_taken %>%
        mutate(component = BlockNumbers[.x]) #puts in the original block ordering so we can know which block was which
      
    }) %>%
    bind_rows( OriginBlock$time_taken %>%
              mutate(component = OriginBlock_number)) 
  
#  test <- fix_z_to_origin(relative_blocks, ArticulationVect) #this is just to see what is added and subtracted
  #The height of each node relative to the origin and normalised
  # node_embeddings <- relative_blocks %>% mutate(elevation_diff = 1,
  #   elevation2 = pull(test, elevation),
  #                                            elevation_diff = elevation - elevation2)
  
  # component_adjust_mat <- adjust_components(g, max_iter =max(relative_blocks$Iter),
  #                                           force = force, flow = flow)
  
  print("Removing multiple articulation nodes")
  node_embeddings <- fix_z_to_origin(relative_blocks, ArticulationVect) %>%
    group_by(node) %>%
    summarise(Iter = first(Iter),
      force = sum(force),
              elevation = first(elevation),
              net_tension = sum(net_tension),
              velocity = sum(velocity)) %>% #the articulation nodes appear multiple times this removes them
    ungroup %>%
    #friction is different depending on the biconnected component so is meaningless in the overall analysis
    #Poissibly could use the drag for the major biconn
    mutate(
      friction = NA,#coef_drag * velocity,
      static_force = force + net_tension,
      net_force = NA,#static_force - friction,
      acceleration = net_force/mass,
      t = 1,
      t = tstep * Iter)
  
  # node_embeddings <- fix_z_to_origin(relative_blocks, ArticulationVect) %>%
  #   group_by(node) %>%
  #   summarise_all(mean) %>%
  #   mutate(Articulation_node = Articulation_node==1)
  
  #combine the node_embeddings and the network_dynamics into a single list
  Out <- list(node_embeddings = node_embeddings, 
              network_dynamics = network_dynamics,
              memory_df = memory_df,
              time_taken = time_taken_df)
  
  return(Out)
  
}