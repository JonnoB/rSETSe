#' remove articulation duplicates
#' 
#' This is a helper function that is used to efficiently aggregate the articulation nodes after embedding
#' 
#' This function uses vectors to be as fast as possible and use as little memory as possible
#' 
#' @param node_embeddings A dataframe, produced after running fix_z_to_origin
#' @param ArticulationVect Vector of articulation nodes
#' 
#' @return A dataframe with with articulation nodes aggregated so that the dataframe has the same number of rows
#' as nodes
#' 


remove_articulation_duplicates <- function(node_embeddings, ArticulationVect) {
  
  embeds <- node_embeddings[(node_embeddings$node %in% ArticulationVect),]
  
  node_vect <- iter_vect <- force_vect <- elevation_vect <- net_tension_vect <- velocity_vect <- rep(NA, length(ArticulationVect))
  
  for(n in 1:length(ArticulationVect)){
    node_vect[n] <- ArticulationVect[n]
    iter_vect[n] <- embeds$Iter[embeds$node == ArticulationVect[n]][1]
    force_vect[n] <- sum(embeds$force[embeds$node == ArticulationVect[n]])
    elevation_vect[n] <- embeds$elevation[embeds$node == ArticulationVect[n]][1]
    net_tension_vect[n] <- sum(embeds$net_tension[embeds$node == ArticulationVect[n]])
    velocity_vect[n] <- sum(embeds$velocity[embeds$node == ArticulationVect[n]])
    
  }
  
  data.frame(
    node = node_vect,
    Iter = iter_vect,
    force = force_vect,
    elevation = elevation_vect,
    net_tension = net_tension_vect,
    velocity = velocity_vect,
    Articulation_node = TRUE)
  
}