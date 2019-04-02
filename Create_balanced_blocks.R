Create_balanced_blocks <- function(g){
  #This function creates a list of biconnected components or blocks.
  #These blocks are balanced such that the connecting vertices contain all the power of the missing part of the network
  #balancing prevents the network reaching a steady state non-zero velocity.
  
  bigraph <- biconnected_components(g)
  
  ArticulationPoints <- get.vertex.attribute(g, "name", bigraph$articulation_points)
  
  
  List_of_BiConComps <-1:length(bigraph$components) %>%
    map(~{
      Comp_num <- .x
      
      Nodes_in_j <- get.vertex.attribute(g, "name", bigraph$components[[Comp_num]]) 
      
      component_art_points <- Nodes_in_j[Nodes_in_j %in% ArticulationPoints]
      
      biconnected_component <- as_data_frame(g, what = "vertices") %>%
        filter(name %in% Nodes_in_j )
      
      balanced_component_df <- 1:length(component_art_points) %>%
        map_df(~{
          
          CurrArt <- component_art_points[.x]
          
          membership_df <- delete.vertices(g, CurrArt) %>% components() %>% .$membership %>%
            tibble(Nodes = names(.), membership = .)
          
          discard_components <- membership_df %>% #Identify the nodes in the target component if articulation poin i is removed
            filter(Nodes %in% Nodes_in_j) %>% pull(membership) %>% unique 
          
          #remove all nodes on the component with nodes from the target block.
          membership_df %>%
            filter( (membership %in% discard_components)) %>% pull(Nodes) %>% {get.vertex.attribute(g, "name") %in% .} %>%
            (1:vcount(g))[.] %>%
            delete_vertices(g, v = .) %>% 
            get.vertex.attribute(., "BalencedPower") %>% sum  %>% #get the net power of the rest of the network
            tibble(name = CurrArt, AuxPower = .)
        }) %>%
        left_join(biconnected_component, ., by = "name") %>%
        mutate(#AuxPower = ifelse(is.na(AuxPower), 0, AuxPower),
               BalencedPower = ifelse(is.na(AuxPower), BalencedPower, AuxPower))
      
      
      Component_j <- {!(get.vertex.attribute(g, "name") %in% balanced_component_df$name)} %>%
        (1:vcount(g))[.] %>%
        delete.vertices(g,.) 
      
      balanced_component <- tibble(name = get.vertex.attribute(Component_j, "name")) %>%
        left_join(balanced_component_df, by = "name") %>%
        pull(BalencedPower) %>%
        set.vertex.attribute(Component_j, name = "BalencedPower", value = .)
      
      return(balanced_component)
    })
  
}
