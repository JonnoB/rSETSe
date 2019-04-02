fix_z_to_origin <- function(relative_blocks, ArticulationVect){
  #This function makes the z values of the different blocks al relative to the same fixed point
  #a single component is chosen, usually the largest and this is given the ID = 0.
  #All the relative distance from origin is then fixed to that block and propergated through
  #the rest of the block tree using the articulation points.
  #relative blocks a data frame containing the outputs from the 
  #ArticulationVect. The articulation nodes of the network.
  
  
  #This is not correct
  # the articulation points are not giving the right answer
  
  target_blocks <- 0
  absolute_blocks <- relative_blocks %>% filter(Reference_ID ==-1)
  Articulation_df <- relative_blocks %>% filter(Reference_ID ==-1)
  n <- 0
  
  next_group <- relative_blocks %>% 
    filter(Reference_ID %in% target_blocks) %>%
    mutate(z = z )
  
  while(n <= nrow(Articulation_df)){
    
    print(n)
    #add new absolute references to the dataframe
    absolute_blocks <- next_group %>%
      bind_rows(absolute_blocks,.)
    
    #add new articulation nodes to the dataframe
    Articulation_df <- next_group %>%
      filter(Articulation_node) %>%
      bind_rows(Articulation_df,.)
    
    #remove all new absolute references from the dataframe
    relative_blocks <- relative_blocks %>%
      filter(!(Reference_ID %in% target_blocks))
    
    n <- n+1
    
    Art_n <-Articulation_df$node[n]
    
    #sumtract art_n relative from all zed scores
    #ass art_n abs to all values
    
    #This has to be done by block not all blocks together
    target_blocks <- relative_blocks %>%
      filter(node == Art_n) %>% pull(Reference_ID)
    
    next_group <- target_blocks %>%
      map_df(~{
      
        temp <- relative_blocks %>% 
          filter(Reference_ID %in% .x) 
        
        local_origin <- temp %>%
          filter(node == Art_n) %>% pull(z) #find the local origin which is the articulation node n
          
        Block_abs_ref <- temp %>% 
          mutate(z = z -local_origin + Articulation_df$z[n]) #add the global absolute reference
        
        return(Block_abs_ref)
          
      }
    )
    
    
    
  }
  
  return(absolute_blocks)
  
}