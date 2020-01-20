#' Create stabilosed blocks
#' 
#' decomposes the network into bi-connected components using articulation points. This speeds up the convergence process
#' and reduces the chances of the SETS algorithm diverging
#' 
#' @param g An igraph object
#' @param Origin block
#' @param OriginBlock_number An integer. this is the origin block chosen from the
#' create_stable_blocks function. Usually this will be the largest block.
#' @param force A character string. This is the node attribute that contains the force the nodes exert on the network.
#' @param flow A character string. This is the edge attribute that is the power flow on the edges.
#' @param edge_name A character string. This is the edge attribute that contains the edge_name of the edges.
#' @param tstep A numeric. The time interval used to iterate through the network dynamics.
#' @param coef_drag A numeric. This sets the multiplier of friction. Only use if you want to be annoyed and confused
#' @param tol A numeric. The tolerance factor for early stopping.
#' @param distance A character string. The edge attribute that contains the original/horizontal distance between nodes.
#' @param max_iter An integer. The maximum number of iterations before stopping. Larger networks usually need more iterations.
#' @param mass A numeric. This is the mass constant of the nodes in normalised networks this is set to 1.
#' @param verbose Logical. This value sets whether messages generated during the process are supressed or not.
#' 
#' @seealso \code{\link{Create_stabilised_blocks}} \code{\link{Find_network_balance}}
#' @return A dataframe with the height embeddings of the network
#' 
#' 
#' @export
Create_stabilised_blocks <- function(g, 
                                     OriginBlock, 
                                     OriginBlock_number, 
                                     force ="net_generation", 
                                     flow = "power_flow", 
                                     edge_name = "edge_name",  
                                     tstep=0.1, 
                                     coef_drag = coef_drag,
                                     tol = 1e-10, 
                                     distance, 
                                     max_iter, 
                                     mass,
                                     sparse,
                                     sample = 1){
  
  #Seperate out the graph into balanced blocks
  #This step will have already been done, but it is fast and simplifies the requirements for the function
  List_of_BiConComps <- create_balanced_blocks(g, force = force, flow = flow)
  
  #remove the Origin block so it doesn't have to be calculated again
  BlockNumbers <-(1:length(List_of_BiConComps))[-OriginBlock_number]
  
  #total in network
  total_force <- sum(abs(get.vertex.attribute(g, force)))
  
  StabilModels <- BlockNumbers %>% 
    map(~{

      Out <- Find_network_balance(List_of_BiConComps[[.x]],
                                  force =force, 
                                  flow = flow, 
                                  tstep = tstep, 
                                  coef_drag = coef_drag,
                                  tol = tol*sum(abs(get.vertex.attribute(List_of_BiConComps[[.x]], force)))/total_force, #the force has to be scaled to the component 
                                  distance = distance, 
                                  edge_name = edge_name,
                                  max_iter =  max_iter, 
                                  mass =  mass, 
                                  sparse = sparse,
                                  sample = sample
                                  )

      #print if the print requirement is on otherwise silent
      #if(!verbose){print(paste("Block" ,.x, "of", max(BlockNumbers) ,"termination", nrow(Out$network_dynamics) )) }
      
      return(Out)
      
    })
  
  #get the block tree of the graph
  Block_tree <- biconnected_components(g)
  
  #extract the articulation nodes
  ArticulationVect <- get.vertex.attribute(g, "name", Block_tree$articulation_points)
  
  #place all nodes relative to the origin
  relative_blocks <- 1:length(StabilModels) %>% 
    map_df(~{
      #print(.x) #It is a bit annoying and pointless now
      StabilModels[[.x]]$node_status %>%
        mutate(Reference_ID = .x)
      
    }) %>%
    bind_rows(OriginBlock$node_status %>% 
                mutate(Reference_ID = 0)) %>%
    mutate(Articulation_node = (node %in% ArticulationVect ))
  
  #get the network_dynamics dataframe for the total calculation
  network_dynamics <- 1:length(StabilModels) %>% 
    map_df(~{
      StabilModels[[.x]]$network_dynamics %>%
        mutate(component = .x)
      
    }) %>%
    bind_rows(OriginBlock$network_dynamics %>% mutate(component = OriginBlock_number))  %>%
  #It can also be useful to get the individual component values out.
    group_by(Iter) %>%
    summarise_all(sum) %>%
    mutate(t = Iter/tstep) %>%
    select(-component)
  
#  test <- fix_z_to_origin(relative_blocks, ArticulationVect) #this is just to see what is added and subtracted
  #The height of each node relative to the origin and normalised
  # node_status <- relative_blocks %>% mutate(elevation_diff = 1,
  #   elevation2 = pull(test, elevation),
  #                                            elevation_diff = elevation - elevation2)
  
  # component_adjust_mat <- adjust_components(g, max_iter =max(relative_blocks$Iter),
  #                                           force = force, flow = flow)
  
  node_status <- fix_z_to_origin(relative_blocks, ArticulationVect) %>%
    group_by(node) %>%
    summarise(Iter = first(Iter),
      force = sum(force),
              elevation = first(elevation),
              net_tension = sum(net_tension),
              velocity = sum(velocity)) %>% #the articulation nodes appear multiple times this removes them
    ungroup %>%
    mutate(
      friction = coef_drag * velocity,
      static_force = force + net_tension,
      net_force = static_force - friction,
      acceleration = net_force/mass,
      t = 1,
      t = tstep * Iter)
  
  # node_status <- fix_z_to_origin(relative_blocks, ArticulationVect) %>%
  #   group_by(node) %>%
  #   summarise_all(mean) %>%
  #   mutate(Articulation_node = Articulation_node==1)
  
  #combine the node_status and the network_dynamics into a single list
  Out <- list(node_status = node_status, network_dynamics = network_dynamics)
  
  return(Out)
  
}