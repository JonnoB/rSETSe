#' Hierarchially cluster based on metric.
#' 
#' Takes a graph and returns a hierarchical cluster object.
#' 
#' @param g An igraph object. The graph to analyse
#' @param weight A character string. The graph edge attribute that will be clustered
#' @param method A character string. The clustering method
#' @export 

hclust_graph <- function(g, weight = "weight", method = "complete"){
  #Hierarchically clusters a graph using a weighted adjacency matrix
  #g, an igraph object
  #weight, the desired weight vector
  #method, the hierachical clustering method to use see hclust for details
  
  weight_vect <- get.edge.attribute(g, weight)
  
  if(is.null(weight_vect)){
    message("Weight vector NULL continuing using edge value 'weight' if available")
  }
  
  distancedf <- distances(g, weights = weight_vect) %>% as_tibble %>% mutate(from = names(.)) %>%
    gather(key = "to", value = "distance",-from)
  
  #graph clustered based edge metric as a distance. Be sure the metric is appropriate!
  distgraph <- distancedf %>% spread(to, distance) %>% select(-from) %>% as.matrix
  rownames(distgraph) <- colnames(distgraph)
  dist_hclust <- distgraph %>% as.dist(., diag = FALSE) %>% hclust(., method) 
  
  return(dist_hclust)
  
}
