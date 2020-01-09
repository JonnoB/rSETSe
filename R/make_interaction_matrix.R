#'Make interaction network
#'
#'This is a helper function for the recombine function used by stabilised block, It creates a sparse adjacency matrix
#'that has the interaction between nodes and the floor/ceiling nodes of the articulation nodes on the shortest path
#'between the origin node in the origin component and all the other nodes in the network.
#'
#'Using the relative blocks is probably not optimal this can be changed when the process works as expected

make_interaction_matrix <- function(target_df, List_of_BiConComps){
  #the edges in the meta graph, these are only the articulation nodes
  
  node_component_df <- 1:length(List_of_BiConComps) %>%
    map_df(~{
      
      as_data_frame(List_of_BiConComps[[.x]], what = "vertices") %>%
        mutate(component = .x)
      
    }) %>%
    rename(node = name)
  
  active_edges <- node_component_df %>%
    left_join(., target_df, by = "component") %>%
    mutate(active_reference = 1) %>%
    filter(complete.cases(.)) %>%
    mutate(node_comp_rows = paste(node, component, sep ="-"),
           node_comp_columns = paste(ref_node, ref_component, sep ="-")) %>%
    select(node_comp_rows, node_comp_columns, active_reference) %>%
    select(1:2)#remove missing these occur in the terminating components
  
  #the nodes in the meta graph, this is all the node/component relationships
  non_active_edges <- node_component_df %>%
    mutate(ref_node = node,
           ref_component = component,
           active_reference = 0) %>%
    mutate(node_comp_rows = paste(node, component, sep ="-"),
           node_comp_columns = paste(ref_node, ref_component, sep ="-")) %>%
    arrange(node) %>%
    select(node_comp_rows, node_comp_columns, active_reference) %>%
    select(1:2)
  
  #The sparse matrix of the  floor/ceiling interaction
  sparse_matrix_result <- graph_from_data_frame(active_edges, #The edges in the meta-graph
                                                directed = T,
                                                vertices = non_active_edges #all the node-component pairs, these acts as the nodes of the metagraph
                                                ) %>%
    as_adjacency_matrix() 
  
  return(sparse_matrix_result)
}
