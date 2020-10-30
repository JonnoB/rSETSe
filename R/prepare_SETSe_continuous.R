#' Prepare continuous network
#' 
#' This function prepares a continuous network for SETSe projection.
#' 
#' The network takes in an igraph object and produces an undirected igraph object that can be used with SETSe/SETSe_auto for embedding.
#'  
#' @param g an igraph object
#' @param node_names a character string. A vertex attribute which contains the node names.
#' @param k The spring constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
#'  the k value will not be added to the network. This is useful k is made through some other processs.
#' @param force_var A node attribute. This is used as the force variable, it must be a numeric or integer value, it cannot have NA's
#' @param sum_to_one Logical. whether the total positive force sums to 1, if FALSE the total is the sum of the positive cases
#' @param distance a positive numeric value. The default is 1
#' 
#' @details 
#'  The function subtracts the mean from all the values so that the system is balanced. If sum_to_one is true then everything is divided by
#'  the absolute sum over two 
#'  
#'  The function adds the node attribute 'force' and the edge attribute 'k' unless k=NULL. The purpose of the function is to easily be able to 
#'  project continuous networks using SETSe. 
#'  
#'  The function creates several variables
#' \itemize{
#'   \item force: a vertex attribute representing the force produced by each node. The sum of this variable will be 0
#'   \item k: The spring constant representing the stiffness of the spring. 
#'   \item edge_name: the name of the edges. it takes the form "from_to" where "from" is the origin node and "to" is the destination node using the 
#'  \code{\link[igraph]{as_data_frame}} function from igraph
#' }
#' 
#' @return A network with the correct edge and node attributes for the embeddings process.
#'
#' @seealso \code{\link{SETSe}}, \code{\link{SETSe_auto}}, \code{\link{prepare_SETSe_binary}}
#' @examples 
#' \dontrun{
#' library(dplyr)
#' embeddings <- biconnected_network %>%
#' #prepare the network for a binary embedding
#' #k is already present in the data so is left null in the preparation function
#' prepare_SETSe_continuous(., node_names = "name", k = NULL, 
#'                         force_var = "force") %>%
#' #embed the network using auto setse
#' #in the biconnected_network dataset the edge weights are used directly as k values
#' SETSe_auto(k = "weight")
#'  }
#' @export

prepare_SETSe_continuous <- function(g, node_names, k = NULL, force_var, sum_to_one = TRUE, distance = 1){
  
  force_var_sym <- rlang::sym(force_var)
  
  g_list <-   igraph::as_data_frame(g, what = "both")
  
  edges_df <- g_list$edges %>%
    dplyr::mutate(distance = distance,
                  edge_name = paste(from, to, sep ="-"))
  
  if(!is.null(k)){
    edges_df <- edges_df %>% 
      dplyr::mutate( k = k)
  }
  
  vertices_df <- g_list$vertices %>% tibble::tibble(.)
  
  vertices_df <- vertices_df %>%
    dplyr::mutate(force = {{force_var_sym}},
                  force = force - mean(force))
  
  if(sum_to_one){
    
    vertices_df <- vertices_df %>%
      dplyr::mutate(force = force/(sum(abs(force))/2))
    
  }
  
  g_out  <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = vertices_df %>% 
                                            dplyr::select(node_names, dplyr::everything())) 
  
  return(g_out)
  
}
