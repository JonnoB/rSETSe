#' Prepare network edges
#' 
#' This function helps prepare the network edges for embedding
#' 
#' @param g an igraph object
#' @param k The spring constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
#'  the k value will not be added to the network. This is useful k is made through some other process.
#' @param distance The spring constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
#'  the distance value will not be added to the network. This is useful distance is made through some other process.
#' @param create_edge_name Logical. Whether to create and edge name attribute or not.
#'  
#'  @details 
#'  The function prepares the edge characteristics of the network so that they can be embedded using the SETSe_ family of functions.
#'  
#'  @return 
#'  The function creates several variables
#' \itemize{
#'   \item k: The spring constant representing the stiffness of the spring. 
#'   \item distance: The minimum distance between nodes. This is the distance between the parallel planes/hyper-planes.
#'   \item edge_name: the name of the edges. it takes the form "from_to" where "from" is the origin node and "to" is the destination node using the 
#'  \link[igraph]{as_data_frame} function from igraph
#' }
#' @family prepare_setse
#' @seealso \link{SETSe}, \link{SETSe_auto}, \link{SETSe_bicomp}, \link{SETSe_auto_hd}
#' @examples
#' set.seed(234) #set the random see for generating the network
#' g <- generate_peels_network(type = "E")
#' embeddings <- g %>%
#' prepare_edges(k = 500, distance = 1) %>%
#' #prepare the network for a binary embedding
#' prepare_SETSe_binary(., node_names = "name",
#'                      force_var = "class") %>%
#' #embed the network using auto setse
#'   SETSe_auto(., force = "group_A")
#' @export
#' 
prepare_edges <- function(g, k = NULL, distance = 1, create_edge_name = TRUE){

  g_list <- igraph::as_data_frame(g, what = "both")
  
  edges_df <- g_list$edges
  
  #set k if necessary
  if(!is.null(k)){
    edges_df <- edges_df %>% 
      dplyr::mutate( k = k)
  }
  
  #set distance if necessary
  if(!is.null(distance)){
    edges_df <- edges_df %>% 
      dplyr::mutate( distance = distance)
  }
  
if(create_edge_name){
    edges_df <- edges_df %>% 
      dplyr::mutate( edge_name = paste(from, to, sep ="-"))
} 
  
  #re-construct the network with the prepared edge data
  g2 <- graph_from_data_frame(edges_df, directed = FALSE, vertices = g_list$vertices)
  
  return(g2)
  
  }