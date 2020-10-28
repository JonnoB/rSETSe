#' Create dataframe of node and aggregated edge embeddings
#'
#' Aggregates edge strain and tension to node level
#' 
#' @param embeddings_data A list. The output of any of the SETSe embedding functions
#' @param function_names A string vector. the names of the aggregation methods to be used
#'
#' @details Often if can be useful to have edge data at node level, an example of this would be plotting
#' the node and tension or strain. To do this requires that the edge embeddings are aggregated somehow to node level
#' and joined to the appropriate node. This function takes as an argument the output of the SETSe embedding functions
#' and any number of aggregation functions to produce a dataframe that is convenient to use.
#' 
#' @return A dataframe with node names, node force, node elevation and strain and tension aggregated useing the named functions.
#' The strain and tension columns are returned with names in the form "strain_x" where "x" is the name of the function used 
#' to aggregate. The total number of columns is dependent on the number of aggregation functions.
#' 
#' @examples
#' 
#' embeddings_data <- two_bicomponents %>%
#'prepare_SETSe_continuous(., node_names = "name", force_var = "force") %>%
#'  SETSe_auto(., k = "weight")
#'
#'out <- create_node_edge_df(embeddings_data, function_names = c("mean", "mode", "sum"))
#' 
#' @export


create_node_edge_df <- function(embeddings_data, function_names = c("mean", "median")){
  #create a named function list of the aggregation function
  function_list <- lapply(function_names, get)
  names(function_list) <- function_names
  
  embeddings_data$edge_embeddings %>% tibble::tibble() %>%
    tidyr::separate(., col = edge_name, into = c("from","to"), sep = "-") %>%
    dplyr::select(from, to, tension, strain) %>%
    # The pivot longer means the nodes at each end of the edge are both associated with the edge.
    tidyr::pivot_longer(cols = c(from, to), names_to = "node_type", values_to = "node") %>%
    dplyr::group_by(node) %>%
    # summarise only tension and strain using the functions named on input
    dplyr::summarise(dplyr::across(.cols = c(tension, strain), 
                                   .fns = function_list  )) %>% 
    dplyr::left_join(embeddings_data$node_embeddings %>% 
                       dplyr::select(node, elevation, force) %>% 
                       tibble::tibble(), 
                     by = "node")
  
  
}