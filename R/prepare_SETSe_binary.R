#' Binary network prepare
#' 
#' This function prepares a binary network for SETSe projection.
#' 
#' The network takes in an igraph object and produces an undirected igraph object that can be used with SETSe/auto_SETSe/SETse_bicomp for embedding.
#'  
#' @param g an igraph object
#' @param node_names a character string. A vertex attribute which contains the node names.
#' @param k The sping constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
#'  the k value will not be added to the network. This is useful k is made through some other processs.
#' @param force_var A node attribute. This is used as the force variable, it must be a character of factor
#' @param positive_value The value in force var that will be counted as the positive value
#' @param sum_to_one Logical. whether the total positive force sums to 1, if FALSE the total is the sum of the positive cases
#' @param distance a positive numeric value. The default is 1
#' 
#' @details The function adds the node attribute 'force' and the edge attribute 'k' unless k=NULL. The purpose of the function is to easily be able to 
#'  project binary networks using SETSe. 
#'  
#'  The function creates several variables
#' \itemize{
#'   \item force: a vertex attribute representing the force produced by each node. The sum of this variable will be 0
#'   \item k: The spring constant representing the stiffness of the spring. 
#'   \item edge_name: the name of the edges. it takes the form "from_to" where "from" is the origin node and "to" is the destination node using the 
#'  \code{\link[igraph]{as_data_frame}} function from igraph
#' }
#' @examples
#' set.seed(234) #set the random see for generating the network
#' g <- generate_peels_network(type = "E")
#' embeddings <- g %>%
#' #prepare the network for a binary embedding
#' prepare_SETSe_binary(., node_names = "name", k = 1000, 
#'                      force_var = "class", 
#'                      positive_value = "A") %>%
#' #embed the network using auto setse
#'   auto_SETSe()
#' @seealso \code{\link{SETSe}}, \code{\link{auto_SETSe}}, \code{\link{SETSe_bicomp}}, \code{\link{prepare_SETSe_continuous}}
#' @export

prepare_SETSe_binary <- function(g, node_names, k = NULL, force_var, positive_value, sum_to_one = TRUE, distance = 1){

g_list <-   as_data_frame(g, what = "both")
  
  edges_df <- g_list$edges %>%
    mutate(distance = distance,
           edge_name = paste(from, to, sep ="-"))
  
  if(!is.null(k)){
    edges_df <- edges_df %>% 
    mutate( k = k)
  }
  
  vertices_df <- g_list$vertices %>% tibble
  
  outcome_var <- vertices_df %>% pull(force_var) %>% {.==positive_value}

  total_pos <- sum(outcome_var)
  total_neg <- sum(!outcome_var)
  
  vertices_df <- vertices_df %>%
    mutate(force = ifelse(outcome_var, 1/total_pos, -1/total_neg ),
           force = ifelse(rep(sum_to_one, nrow(.)), force, force*total_pos))
  
  g_out  <- graph_from_data_frame(edges_df, directed = FALSE, 
                                  vertices = vertices_df %>%
                                    select(node_names, 
                                           everything())
  ) 
  
  return(g_out)
  
}
