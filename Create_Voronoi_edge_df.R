#' Create a network from a voronoi tessalation
#' 
#' This function allows voronoi tessalations to be converted into networks. This is no fun to do hence this function.
#' 
#' This function is useful when you want to turn a collection of spatial points into a network.
#' First the points are converted into a voronoi tesselation then this function converts them into a network.
#' 
#' @param Voronoi_sf An Sf object that describes the voronoi geometry and has a column that give the node that the cell is based around
#' @param Node_name Node_name a bare unquoted name of the Node column in the sf object
#' 
#' @export


Create_Voronoi_edge_df <- function(Voronoi_sf, Node_name){
  #Voronoi_sf an Sf object that describes the voronoi geometry and has a column that give the node that the cell is based around
  #Node_name a bare unquoted name of the Node column in the sf object
  
  #The function takes an sf object and outputs a dataframe of three columns, which are the from and to Nodes and the edge name created by
  #concatenating from and two seperated by "-">
  
  #all that is needed here is the node and the geometry of the node
  test <- Voronoi_sf %>%
    st_intersects( sparse = FALSE)
  
  name_vect <- pull(Voronoi_sf, {{Node_name}})
  
  #takes the dense intersection and produces the edge list of the network
  test2 <- as_tibble(test) %>%
    set_names(name_vect) %>%
    mutate(Node = name_vect) %>%
    gather(key = Node2, value = Link,  -{{Node_name}}) %>%
    filter(Link,  {{Node_name}} !=Node2) %>%
    graph_from_data_frame(directed = FALSE) %>%
    as_data_frame() %>%
    mutate(Link = paste(from, to, sep = "-"))
  
}