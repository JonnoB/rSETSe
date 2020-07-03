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
  absolute_blocks <- Articulation_df <- relative_blocks[relative_blocks$Reference_ID ==-1,]
  
  #This vector stors articulation nodes that have already been used,
  #This stops the same node being recycled which causes errors
  #The vector begins as an empty vector
  used_articulation_nodes<- as.numeric(list())
  
  next_group <- relative_blocks[relative_blocks$Reference_ID %in% target_blocks,]
  
  
  # vectors that are used
  # reference_id <- relative_blocks$Reference_ID
  # node_vect <- relative_blocks$node
  # articulation_df <- relative_blocks$Articulation_node
  
  #The previous while used the rows in Articulation_df but this changed and so could be confusing. As there can only be the same 
  #amount of loops as there are articulation nodes. n must be equal to or smaller than the number of articulation nodes - 1
  #nrow(Articulation_df)
  #while(n <= (length(ArticulationVect)-2)){
  
  while(length(used_articulation_nodes)< length(ArticulationVect)){
    #  while(Art_n != 190){
    
    #add new absolute references to the dataframe
    absolute_blocks <- bind_rows(absolute_blocks, next_group)
    
    #add new articulation nodes to the dataframe
    {
      Articulation_df <- bind_rows(Articulation_df, next_group[next_group$Articulation_node,] )
      Articulation_df <-Articulation_df[!(Articulation_df$node %in% used_articulation_nodes),]  #remove all previously used articulation nodes
    }
    
    #The removal of previously used articulation nodes is important as all articulation nodes are in the network at least twice
    
    #remove all new absolute references from the dataframe
    relative_blocks <- relative_blocks[!(relative_blocks$Reference_ID %in% target_blocks),] 
    
    #get the next articulation node in the queue
    Art_n <-Articulation_df$node[1] #
    #print(Art_n)
    #subtract art_n relative from all elevation scores
    #ass art_n abs to all values
    
    #This has to be done by block not all blocks together
    target_blocks <- relative_blocks$Reference_ID[relative_blocks$node == Art_n] #%>%
    #filter(node == Art_n) %>% pull(Reference_ID)
    
    next_group <- map_df(.x = target_blocks, ~{
      
      temp <- relative_blocks[relative_blocks$Reference_ID == .x,] 
      
      local_origin <- temp$elevation[temp$node ==Art_n] #find the local origin which is the articulation node n
      
      temp$elevation <- temp$elevation -local_origin + Articulation_df$elevation[1] #add the global absolute reference
      
      return(temp)
      
    }
    
    )
    
    used_articulation_nodes <- c(used_articulation_nodes, Articulation_df$node[1])
    
    
  }
  #The final block is added on here.
  #if I rearrange the process I can get it all in the loop
  #But I just can't be bothered right now :'(
  absolute_blocks <- bind_rows(absolute_blocks,next_group ) 
  
  return(absolute_blocks)
  
}
