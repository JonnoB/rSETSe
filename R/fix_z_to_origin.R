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
  
  #This is an integer vector of blocks that have been converted to absolute values
  #This is one of the growing vectors.
  absolute_blocks <- c(0)
  
  #The first block to be made absolute is the origin block indexed at 0
  target_blocks <- 0
  
  #The queued articulation nodes.
  #These are articulation nodes that have been converted to absolute elevation in at least one of thier blocks
  #They can be used to make subsequent blocks absolute
  queued_articulation_nodes <- as.numeric(list())
  
  #This vector stors articulation nodes that have already been used,
  #This stops the same node being recycled which causes errors
  #The vector begins as an empty vector
  used_articulation_nodes <- as.numeric(list())
  
  #The block reference number.The origin block is always 0
  reference_id <- relative_blocks$Reference_ID
  #The node IDs
  node_vect <- relative_blocks$node
  #the aarticulation vectors
  articulation_vect <- relative_blocks$Articulation_node
  #The node elevation vector. This is converted block by block from relative elevation to 
  #absolute elevation
  elevation_vect <- relative_blocks$elevation
  
  #The while loop continues as long as used articulation nodes vector is smaller than the total number
  #of articulation nodes in the network
  while(length(used_articulation_nodes) < length(ArticulationVect)){
    #while(sum(used_articulation_nodes %in% ArticulationVect) == length(ArticulationVect)){
    #  while(Art_n != 190){
    
    #add new articulation nodes to the queue
    #It has four logical conditions
    #1 the nodes must be in the target block
    #2 The nodes must be an articulation node
    #3 the node cannot already be in the queue
    #4 cannot be a node that has already been used
    
    #The removal of previously used articulation nodes is important as all articulation nodes are in the network at least twice.
    # not removing later occurances can lead to levelling errors and crazy results
    #The 4th logical constraint is becuase the target blocks from the previous round also contain the previous rounds
    #Art_n node (active node), this means it will be added in again.
    queued_articulation_nodes <- c(queued_articulation_nodes, node_vect[reference_id %in% target_blocks & 
                                                                          articulation_vect &
                                                                          !(node_vect %in% queued_articulation_nodes) &
                                                                          !(node_vect %in% used_articulation_nodes)])
    #print(queued_articulation_nodes)
    
    #get the next articulation node in the queue
    #This is the active articulation node
    Art_n <-queued_articulation_nodes[1] #
    
    
    #add the current articulation node to the vector of used nodes
    used_articulation_nodes <- c(used_articulation_nodes, Art_n)
    #print(Art_n)
    #subtract art_n relative from all elevation scores
    #ass art_n abs to all values
    
    #The blocks that this articulation node is in excluding the current active block
    #1 The block is not already absolute
    #2 The node is the active articulation node
    #Only unique values are used, this is more a security blanket than anything. I don't know if it is necessary
    target_blocks <- unique(reference_id[!(reference_id %in% absolute_blocks) & node_vect == Art_n])
    
    
    #This is the value that the nodes in the target blocks will be adjusted by
    #It is the elevation that matches the following conditions
    #1 The reference id of the block has to already been adjusted to absolute terms
    #2 The node has to be the articulation node being adjusted
    elevation_adjust <-  unique(elevation_vect[(reference_id %in% absolute_blocks) & node_vect == Art_n])
    
    for(n in target_blocks){
      #print(n)
      #The local origin is the active articulation node's elevation in the target block
      #This is the node which is both
      #1 in the target block
      #2 the active articulation node
      local_origin <- elevation_vect[(reference_id %in% n) & node_vect == Art_n]
      
      #The elevation of all nodes in the target block have the local origin remove to place everything relative to 0
      #then the absolute height of the articulation node is added making the entire bi-connected component elevation absolute.
      elevation_vect[(reference_id %in% n)] <- elevation_vect[(reference_id %in% n)] - local_origin + elevation_adjust
      
    }
    
    #add to the absolute blocks vector
    absolute_blocks <- c(absolute_blocks, target_blocks)
    
    #The active articulation node is removed from the queue
    queued_articulation_nodes <- queued_articulation_nodes[-1]

  }
  #The elevation vector is now absolute and can be inserted back into the original datframe
  relative_blocks$elevation <- elevation_vect
  
  return(relative_blocks)
  
}
