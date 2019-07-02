angle_from_solved_heights <- function(height_solution, g = IEEE_118 ){
  #This function creates a dataframe with the network theta from a solved height dataframe.
  max_height_diff <- height_solution %>% 
    arrange(z) %>%
    slice(1, n()) 
  
  max_height_diff2 <- max(max_height_diff$z)-min(max_height_diff$z)
  
  #indexing only works on the IEEE graphs as the node name is ther vertex number
  Adjacent <- distances(g, v = max_height_diff$node[1], to = max_height_diff$node[2])
  
  Graph_theta <- atan(max_height_diff2/Adjacent)[1,1]
  
  Graph_theta*360/(2*pi)
  
  Out <- tibble(
    theta_rads= Graph_theta,
    theta_degs = Graph_theta*360/(2*pi)
  )
  
  return(Out)
}
