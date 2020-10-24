#'Make interaction network
#'
#'This is a helper function and is called by 'adjust_components' which is itself called by 'Create_stabilised_blocks_expanded' 
#'used by stabilised block, It creates a sparse adjacency matrix that has the interaction between nodes and the 
#'floor/ceiling nodes of the articulation nodes on the shortest path between the origin node in the origin 
#'component and all the other nodes in the network.
#'
#' @param target_df A dataframe... I need to check what of
#' @param node_component_df A dataframe representing the embedded nodes of a bi-conneceted component
#'
#'@details The networks created are bipartite
#'
#'Using the relative blocks is probably not optimal this can be changed when the process works as expected
#'

make_interaction_matrix <- function(target_df, node_component_df){
  #the edges in the meta graph, these are only the articulation nodes
  
  
  active_edges <- node_component_df %>%
    dplyr::left_join(., target_df, by = "component") %>%
    dplyr::mutate(active_reference = 1) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(node_comp_rows = paste(node, component, sep ="-"),
           node_comp_columns = paste(ref_node, ref_component, sep ="-")) %>%
    dplyr::select(node_comp_rows, node_comp_columns, active_reference) %>%
    dplyr::select(1:2)#remove missing these occur in the terminating components
  
  #the nodes in the meta graph, this is all the node/component relationships
  non_active_edges <- node_component_df %>%
    dplyr::mutate(ref_node = node,
           ref_component = component,
           active_reference = 0) %>%
    dplyr::mutate(node_comp_rows = paste(node, component, sep ="-"),
           node_comp_columns = paste(ref_node, ref_component, sep ="-")) %>%
    dplyr::arrange(node) %>%
    dplyr::select(node_comp_rows, node_comp_columns, active_reference) %>%
    dplyr::select(1:2)
  
  #The sparse matrix of the  floor/ceiling interaction
  sparse_matrix_result <- igraph::graph_from_data_frame(active_edges, #The edges in the meta-graph
                                                directed = T,
                                                vertices = non_active_edges #all the node-component pairs, these acts as the nodes of the metagraph
                                                ) %>%
    igraph::as_adjacency_matrix() 
  
  return(sparse_matrix_result)
}
