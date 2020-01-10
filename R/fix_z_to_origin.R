#' fix elevation to origin
#' 
#' This function makes the elevation of the different blocks all relative to the same fixed point, 
#' this allows easy use of bi-connected components This is a helper function only.
#' z is used instead of elevation for legacy issues. This can be changed when I can be bothered
#' 
#' 
#' The function needs to be renamed to make it more SETSe friendly something like fix_elevation_to_origin. But I need to see what it will impact first
#' 
#' @param relative_blocks a data frame containing the outputs from the 
#' @param ArticulationVect The articulation nodes of the network.

fix_z_to_origin <- function(relative_blocks, ArticulationVect){

  #This function needs some further explaining!
  
  target_blocks <- 0
  absolute_blocks <- relative_blocks %>% filter(Reference_ID ==-1)
  Articulation_df <- relative_blocks %>% filter(Reference_ID ==-1)
  n <- 0
  
  next_group <- relative_blocks %>% 
    filter(Reference_ID %in% target_blocks) %>%
    mutate(elevation = elevation ) #can this be removed? It doesn't seem to do anything
  
  while(n <= nrow(Articulation_df)){
    
    #add new absolute references to the dataframe
    absolute_blocks <- next_group %>%
      bind_rows(absolute_blocks,.) #this should be pre-allocated
    
    #add new articulation nodes to the dataframe
    Articulation_df <- next_group %>%
      filter(Articulation_node) %>%
      bind_rows(Articulation_df,.)
    
    #remove all new absolute references from the dataframe
    relative_blocks <- relative_blocks %>%
      filter(!(Reference_ID %in% target_blocks))
    
    n <- n+1
    #get the nth articulation node in the queue
    Art_n <-Articulation_df$node[n]
    
    #subtract art_n relative from all elevation scores
    #ass art_n abs to all values
    
    #This has to be done by block not all blocks together
    target_blocks <- relative_blocks %>%
      filter(node == Art_n) %>% pull(Reference_ID)
    
    next_group <- target_blocks %>%
      map_df(~{
      
        temp <- relative_blocks %>% 
          filter(Reference_ID %in% .x) 
        
        local_origin <- temp %>%
          filter(node == Art_n) %>% pull(elevation) #find the local origin which is the articulation node n
          
        Block_abs_ref <- temp %>% 
          mutate(elevation = elevation -local_origin + Articulation_df$elevation[n]) #add the global absolute reference
        
        return(Block_abs_ref)
          
      }
    )
    
    
    
  }
  
  return(absolute_blocks)
  
}