

adjust_components <- function(g, max_iter, force, flow){
  
  
  #The origin block is the largest connected component
  
  List_of_BiConComps <- create_balanced_blocks(g, 
                                               force = force, 
                                               flow = flow)
  
  #find the largest component and use that as the origin block
  OriginBlock_number <-List_of_BiConComps %>% map_dbl(~vcount(.x)) %>% which.max()
  
  biconnected_g_info <- biconnected.components(g)
  
  #the articulations nodes and their graph vertex ids
  articulation_nodes_df <-tibble(name = names(biconnected_g_info$articulation_points), id = biconnected_g_info$articulation_points)
  
  ##get the first node that is NOT an articulation point
  #This is done my numbering all the nodes from 1 to the total nodes in the component. then finding the names of the node in the origin block
  #that are also biconnected components. These nodes are removed and the the node with the smallest ID is the origin point from which
  #distances and shortes paths will be measured.
  origin_block_base_node_name <-get.vertex.attribute(g,"name")[-which(get.vertex.attribute(List_of_BiConComps[[OriginBlock_number]],"name") %in% articulation_nodes_df$name)] %>% 
    min
  
  #The distance from the base node is calculated here in a named vector,
  distance_from_origin <- distances(g)[articulation_nodes_df$id, origin_block_base_node_name]
  
  #articulation nodes in each component are found here by looping through each 
  #component and seeing which of the articulation nodes are present
  node_bicomp_relation <- 1:length(biconnected_g_info$components) %>% 
    map_df(~{
      
      articulation_nodes_df %>%
        filter(id %in% biconnected_g_info$components[[.x]]) %>%
        mutate(component = .x)
      
      
    }) %>%
    left_join(., tibble(name = names(distance_from_origin), distance = distance_from_origin), by = "name") %>%#add in distance from origin data
    group_by(component) %>%
    mutate(is_floor = (min(distance)== distance) & OriginBlock_number != component ) %>% #mark the node that is closest to the origin in each component
    ungroup
  
  #Using the shortest path from the origin node to the articulation point find all the floor and ceiling nodes
  #for a given destination point each component has a pair of articulations points
  #one is the floor the other is the ceiling.
  #the floor node is the node that is closest to the origin the ceiling is the node that is furthest away.
  
  #calculates the shortest path for every floor node in each component
  #As the components are bi-components are bi-connected and the distance is from a single point
  #components can only have a single and unique floor node. However, that floor node
  #can be the floor node of multiple components
  bicomp_shortest_paths <- shortest_paths(g, 
                                          from = origin_block_base_node_name, 
                                          to = unique(node_bicomp_relation$name[node_bicomp_relation$is_floor]),
                                          output = "both")
  
  #Which compenents are on the shortest path of the component being checked component being checked?
  #This is important as a node can be a floor for multiple components. We don't want to re-level for all the components
  #the node is part of. We only want components for which there are two nodes on the shortest path as well as the terminating 
  #node becuase that is in the final component
  
  #This requires another loop. which checks to see if the component is part of the shortest path route for the nodes
  
  floor_df_2 <- ceiling_df_2 <- node_bicomp_relation  %>%
    mutate(path_component = NA,
           path_component2 = NA,
           target_node_name = NA) %>% slice(0)
  
  for(i in 1:length(bicomp_shortest_paths$vpath) ){
    target_shortest_path <-i
    
    #component_specific_floor_ceiling_nodes <- bicomp_shortest_paths$vpath[[target_shortest_path]]
    node_path <- bicomp_shortest_paths$vpath[[target_shortest_path]]
    edge_path <- bicomp_shortest_paths$epath[[target_shortest_path]] #this may not be necessary and then output can be "vpath" only
    
    component_specific_floor_ceiling_nodes_2 <- names(bicomp_shortest_paths$vpath[[target_shortest_path]])[-1]
    
    template_nodes_df <-node_bicomp_relation %>%
      ungroup() %>%
      group_by(component) %>%
      mutate(path_component = (name %in% component_specific_floor_ceiling_nodes_2)*1,
             path_component2 = sum(path_component)) %>%
      ungroup %>%
      mutate(target_node_name = names(node_path[length(node_path)]))
    
    floor_df <- template_nodes_df  %>%
      filter(path_component==1, #the node has to be on the shortest path
             (path_component2 > 1 | #And the component the node is in needs to have a floor and ceiling
                name == target_node_name) &  # or be the target node
               is_floor)#and all the nodes need to be floors
    
    #This finds the ceilings using the floor dataframe to subset the data
    ceiling_df <- template_nodes_df %>%
      filter((component %in% floor_df$component) | #The component has to be one of the components that contains a floor
               component ==OriginBlock_number, #or it has to be the origin block component
             id %in% floor_df$id, #And the node id has to be a node that is also a floor
             !is_floor) #But the node cannot be in a component where it is actually a floor
    
    floor_df_2 <- bind_rows(floor_df_2, floor_df)
    ceiling_df_2 <- bind_rows(ceiling_df_2, ceiling_df)
  }
  
  #these two peieces of code create the data frames that provide the references for all
  
  #The filtering logic is a bit tricky here might require various networks to test on
  ceiling_df <- node_bicomp_relation %>%
    left_join(ceiling_df_2 %>% select(ref_node = name, ref_component = component, name = target_node_name) , by = "name") %>%
    filter(component !=  OriginBlock_number,
           component != ref_component) %>% #no ceiling can be in the same component as the target component
    distinct(component, ref_node, ref_component) %>%
    make_interaction_matrix(., List_of_BiConComps) #calls the interaction function to create the final sparse matrix
  
  floor_df <- node_bicomp_relation %>%
    left_join(floor_df_2 %>% select(ref_node = name, ref_component = component, name = target_node_name) , by = "name") %>%
    filter(component !=  OriginBlock_number,
           is_floor) %>% #A floor has to logically be a floor
    distinct(component, ref_node, ref_component) %>%
    make_interaction_matrix(., List_of_BiConComps) #calls the interaction function to create the final sparse matrix
  
  
  #expand the sparse matrix to block diagnonal form. This means each time iteration is only multiplied by itself. 
  block_diag_floor <- bdiag(lapply(1:(max_iter+1), function(n){floor_df}))
  block_diag_ceiling <- bdiag(lapply(1:(max_iter+1), function(n){ceiling_df}))
  
  return(list(floor = block_diag_floor, ceiling = block_diag_ceiling))
  
}