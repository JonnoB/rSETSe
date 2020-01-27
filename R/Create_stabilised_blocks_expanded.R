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
Create_stabilised_blocks_expanded <- function(g, 
                                     OriginBlock, 
                                     OriginBlock_number, 
                                     force ="net_generation", 
                                     flow = "power_flow", 
                                     edge_name = "edge_name",  
                                     tstep=0.1, 
                                     coef_drag = coef_drag,
                                     tol = 0, 
                                     distance, 
                                     max_iter, 
                                     mass, 
                                     verbose = FALSE,
                                     sparse){
  
  #Seperate out the graph into balanced blocks
  #This step will have already been done, but it is fast and simplifies the requirements for the function
  List_of_BiConComps <- create_balanced_blocks(g, force = force, flow = flow)
  
  #remove the Origin block so it doesn't have to be calculated again
  BlockNumbers <-(1:length(List_of_BiConComps))[-OriginBlock_number]
  
  StabilModels <- BlockNumbers %>% 
    map(~{

      Out <- SETSe_expanded(List_of_BiConComps[[.x]],
                                  force =force, 
                                  flow = flow, 
                                  tstep = tstep, 
                                  coef_drag = coef_drag,
                                  tol = tol, 
                                  distance = distance, 
                                  edge_name = edge_name,
                                  max_iter =  max_iter, 
                                  mass =  mass, 
                                  sparse = sparse,
                                  verbose = verbose)

      return(Out)
      
    })
  
  #get the block tree of the graph
  Block_tree <- biconnected_components(g)
  
  #extract the articulation nodes
  ArticulationVect <- names(Block_tree$articulation_points)
  
  full_list <- list()
  full_list[BlockNumbers] <- StabilModels
  full_list[[OriginBlock_number]] <- OriginBlock
  message("make relative blocks")

  
 #mark each node with the correct component ID
  relative_blocks <- 1:length(full_list) %>% 
    map_df(~{
      
     # print(nrow(full_list[[.x]])/(max_iter+1))
      full_list[[.x]] %>%
        mutate(component = .x)
      
    }) %>%
    ungroup %>%
    mutate(Articulation_node = (node %in% ArticulationVect )) %>%
    arrange(Iter, node) #arrange to be in the same order as the block diagram

  
  message("creating adjustment matrices")
  component_adjust_mat <- adjust_components(g, max_iter = max_iter,
                                            force = force, flow = flow)
  

  #print(table(relative_blocks$Iter))
  #print((dim(component_adjust_mat$floor)))
  #print(table(relative_blocks$node))
  message("calculate adjusted values and aggregate")
  # node_emebeddings <- relative_blocks %>%
  #   mutate(elevation = elevation +  as.vector(component_adjust_mat$ceiling %*% relative_blocks$elevation) - 
  #            as.vector(component_adjust_mat$floor %*% relative_blocks$elevation))
  
  #This aggregates the key components of the matrix. by the appropriate metric. all other metrics can just be created from these
  # node_emebeddings <- tibble(
  #   node = rep(component_adjust_mat$node_order, times = max(relative_blocks$Iter+1)),
  #   force = harmonise_and_aggregate(relative_blocks$force, component_adjust_mat, harmonise_vector = FALSE, colsum = TRUE),
  #   elevation = harmonise_and_aggregate(relative_blocks$elevation, component_adjust_mat, harmonise_vector = TRUE, colsum = FALSE),
  #   net_tension = harmonise_and_aggregate(relative_blocks$net_tension, component_adjust_mat, harmonise_vector = FALSE, colsum = TRUE),
  #   velocity = harmonise_and_aggregate(relative_blocks$velocity, component_adjust_mat, harmonise_vector = TRUE, colsum = TRUE)) %>%
  #   mutate(
  #     friction = coef_drag*velocity,
  #     static_force = force+net_tension,
  #     net_force = static_force - friction,
  #     acceleration = net_force/mass,
  #     t = 1,
  #     Iter = rep(0:(max(relative_blocks$Iter)), each = vcount(g)),
  #     t = tstep*Iter)

  
  #Change the node height in all blocks to be correct relative to the originblock
  #place all nodes relative to the origin
  #The height of each node relative to the origin and normalised  
#  message("aggregating across iterations and nodes")
  node_emebeddings <- relative_blocks %>%
   # filter(node=="A") %>%
  # filter(Iter ==(max_iter)) %>% #Delete this once the aggregation is sorted!
    mutate(elevation = harmonise_and_aggregate(relative_blocks$elevation, component_adjust_mat, harmonise_vector = TRUE, colsum = NULL),
           velocity = harmonise_and_aggregate(relative_blocks$velocity, component_adjust_mat, harmonise_vector = TRUE, colsum = NULL)) %>%
    group_by(Iter, node) %>%
    summarise(force = sum(force),
              elevation = first(elevation),
              net_tension = sum(net_tension),
              velocity = sum(velocity)) %>% #the articulation nodes appear multiple times this removes them
  #  mutate(Articulation_node = Articulation_node==1)    %>%
    ungroup %>%
      mutate(
        friction = coef_drag*velocity,
        static_force = force+net_tension,
        net_force = static_force - friction,
        acceleration = net_force/mass,
        t = 1,
        Iter = rep(0:(max(relative_blocks$Iter)), each = vcount(g)),#vcount(g)),
        t = tstep*Iter)
  
  
  return(node_emebeddings)
  
}