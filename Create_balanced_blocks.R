Create_balanced_blocks <- function(g, force = "BalencedPower"){
  #This function creates a list of biconnected components or blocks.
  #These blocks are balanced such that the connecting vertices contain all the power of the missing part of the network
  #balancing prevents the network reaching a steady state non-zero velocity.
  
  bigraph <- biconnected_components(g)
  
  ArticulationPoints <- get.vertex.attribute(g, "name", bigraph$articulation_points)
  
  
  List_of_BiConComps <-1:length(bigraph$components) %>%
    map(~{
      Comp_num <- .x
      
      #nodes in the current block
      Nodes_in_j <- get.vertex.attribute(g, "name", bigraph$components[[Comp_num]]) 
      
      #The nodes in the current block that are articulation points. 
      #There will be at least one articulation point, unless the whol network is one block
      #In a two node block that does not terminate at an end both nodes are articulation points.
      component_art_points <- Nodes_in_j[Nodes_in_j %in% ArticulationPoints]
      
      biconnected_component <- as_data_frame(g, what = "vertices") %>%
        filter(name %in% Nodes_in_j )
      
      #This is used in the case that the block is only two nodes connected by a single edge
      edge_df <- as_data_frame(g) %>%
        filter(from %in% Nodes_in_j,to %in% Nodes_in_j ) 

      
      balanced_component_df <- 1:length(component_art_points) %>%
        map_df(~{
          
          CurrArt <- component_art_points[.x]
          
          membership_df <- delete.vertices(g, CurrArt) %>% components() %>% .$membership %>%
            tibble(Nodes = names(.), membership = .)
          
          #Identify the components that contain nodes from the same block as the articulation point. These need to be discarded
          discard_components <- membership_df %>% 
            filter(Nodes %in% Nodes_in_j) %>% pull(membership) %>% unique #the component membership of the remaining nodes
          
          #remove all nodes on the component with nodes from the target block.
          membership_df %>%
            filter( (membership %in% discard_components)) %>% pull(Nodes) %>% {get.vertex.attribute(g, "name") %in% .} %>%
            (1:vcount(g))[.] %>%
            delete_vertices(g, v = .) %>% 
            get.vertex.attribute(., force) %>% sum  %>% #get the net power of the rest of the network
            tibble(name = CurrArt, 
                   AuxPower = .) #The Auxilary power is all the power not including the node itself
        }) %>%
        left_join(biconnected_component, ., by = "name") %>%
        mutate(temp = ifelse(is.na(AuxPower), !!sym(force), AuxPower)) #
      
      #makes the nodes of the two node 1 edge block have the same power as that which flows between them.
      #this forces a balance
      #The nodes are by default in alphabetical order so the polarity of the power flow will always be correct
      if(nrow(edge_df)==1){
        balanced_component_df <- balanced_component_df %>%
          mutate(temp = edge_df %>%
                   pull(., PowerFlow) %>% c(., -.) %>% round(., 10)) #rounding stops tiny values
      }
      
      Component_j <- {!(get.vertex.attribute(g, "name") %in% balanced_component_df$name)} %>%
        (1:vcount(g))[.] %>%
        delete.vertices(g,.) 
      
      balanced_component <- tibble(name = get.vertex.attribute(Component_j, "name")) %>%
        left_join(balanced_component_df, by = "name") %>%
        pull(temp) %>%
        set.vertex.attribute(Component_j, name = force, value = .)
      
      return(balanced_component)
    })
  
}
