#' Prepare continuous network high dimensional equivalent
#' 
#' This function prepares a continuous network for SETSe projection.
#' 
#' The network takes in an igraph object and produces an undirected igraph object that can be used with SETSe/SETSe_auto for embedding.
#'  
#' @param g an igraph object
#' @param node_names a character string. A vertex attribute which contains the node names.
#' @param k The spring constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
#'  the k value will not be added to the network. This is useful k is made through some other process.
#' @param force_var A character vector. This is the vector of node attributes to be used as the force variables.  
#'  All the attributes must be a numeric or integer value, and cannot have NA's. On a single variable embedding this is usually "force"
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
#'  \link[igraph]{as_data_frame} function from igraph
#' }
#' 
#' @return A network with the correct edge and node attributes for the embeddings process.
#' @family prepare_setse
#' @seealso \link{SETSe_auto_hd}
#' @examples 
#' embeddings <- biconnected_network %>%
#' #prepare the network for a binary embedding
#' #k is already present in the data so is left null in the preparation function
#' prepare_SETSe_continuous(., node_names = "name", k = NULL, 
#'                         force_var = "force") %>%
#' #embed the network using auto setse
#' #in the biconnected_network dataset the edge weights are used directly as k values
#' SETSe_auto(k = "weight")
#' @export

prepare_SETSe_continuous_hd <- function(g, node_names, k = NULL, force_var, sum_to_one = TRUE, distance = 1){
  #This function can be improdved by removing the edges from the loop as they only need to be done once
  #also 
  g2 <- g
  #cycles through each of the variables in the the force_var parameter
  for(x in 1:length(force_var)){
  
  force_var_sym <- rlang::sym(force_var[x])
  
  g_list <-   igraph::as_data_frame(g2, what = "both")
  
  edges_df <- g_list$edges %>%
    dplyr::mutate(distance = distance,
                  edge_name = paste(from, to, sep ="-"))
  #set k if necessary
  if(!is.null(k)){
    edges_df <- edges_df %>% 
      dplyr::mutate( k = k)
  }
  
  #extract the vertices dataframe as a tibble
  vertices_df <- g_list$vertices %>% tibble::tibble(.)
  
  #create a temporary variable 'force' ans assign the force variable to it.
  #subtract the mean force value so that the forces of the system a are balanced
  vertices_df <- vertices_df %>%
    dplyr::mutate(
      {{force_var_sym}} := {{force_var_sym}} - mean({{force_var_sym}}))
  
  #If the sum to one option is set to TRUE, divide the force by the absolute sum of the forces over 2.
  if(sum_to_one){
    
    vertices_df <- vertices_df %>%
      dplyr::mutate({{force_var_sym}} := {{force_var_sym}}/(sum(abs({{force_var_sym}}))/2))
    
  }
  
  #re-build the graph with the the force prepared
  g2  <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = vertices_df %>% 
                                            dplyr::select(dplyr::all_of(node_names), dplyr::everything())) 
  
  }
  
  return(g2)
  
}
