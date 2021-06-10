#' Prepare categorical features for embedding
#' 
#' This function prepares a binary network for SETSe projection.
#' 
#' The network takes in an igraph object and produces an undirected igraph object that can be used with the embedding functions.
#'  
#' @param g an igraph object
#' @param node_names a character string. A vertex attribute which contains the node names.
#' @param force_var A vector of force attributes. This describes all the categorical force attributes of the network. 
#'  All named attributes must be either character or factor attributes.
#' @param sum_to_one Logical. whether the total positive force sums to 1, if FALSE the total is the sum of the positive cases
#' 
#' @details The purpose of the function is to easily be able to project categorical features using SETSe. The function creates new variables
#'  where each variable represents one level of the categorical variables. For embedding only n-1 of the levels are needed.
#'  
#'  The function creates several variables of the format "force_". Vertex attributes representing the force produced by each node 
#'  for each categorical value, there will be n of these variables representing each level of the categorical values. The variable names 
#'  will be the the name of the variable and the name of the level seperated by and underscore. For example, with a variable group and levels A and B, the created force variables will be
#'  "group_A" and "group_B" The sum of these variables will be 0.
#' 
#' @return A network with the correct node attributes for the embeddings process.
#' @family prepare_setse
#' @seealso \link{setse}, \link{setse_auto}, \link{setse_bicomp}, \link{setse_auto_hd}
#' @examples
#' set.seed(234) #set the random see for generating the network
#' g <- generate_peels_network(type = "E")
#' embeddings <- g %>%
#' prepare_edges(k = 500, distance = 1) %>%
#' #prepare the network for a binary embedding
#' prepare_categorical_force(., node_names = "name",
#'                      force_var = "class") %>%
#' #embed the network using auto_setse
#'   setse_auto(., force = "class_A")
#'
#' @export

prepare_categorical_force <- function(g, node_names, force_var,  sum_to_one = TRUE){


  
  g_list <-   igraph::as_data_frame(g, what = "both")
  
  vertices_df <- g_list$vertices %>% tibble::tibble(.)
  
  #cycles through all the categorical force variables and converts them
  one_hot_tibble<- 1:length(force_var) %>%
    purrr::map(~{
      
      force_var_sym <- rlang::sym(force_var[.x])
      
      #create a one hot vector encoded tibble for all the categorical levels
      one_hot_tibble <- vertices_df %>%
        dplyr::select({{force_var_sym}}) %>% dplyr::mutate(value = 1,
                                                           id = 1:dplyr::n()) %>%
        tidyr::pivot_wider(names_from = {{force_var}}, values_from = value, values_fill = 0) %>%
        dplyr::select(-id) %>%
        rlang::set_names(paste(force_var, names(.), sep = "_")) %>%
        dplyr::mutate(dplyr::across(.fns=~{.x-mean(.x)}))
      
      #adjust so that each side sums to one if desired
      if(sum_to_one){
        
        one_hot_tibble <- one_hot_tibble %>%
          dplyr::mutate(dplyr::across(.fns = ~{.x/(sum(abs(.x))/2)}))
        
      }
      
    }) %>%
    dplyr::bind_cols(.)
  

  
  #bind the balanced one hot vector encoded variables back onto the node dataframe
  vertices_df <- vertices_df %>% 
    dplyr::bind_cols(one_hot_tibble)
  
  #re-build the network
  g_out  <- igraph::graph_from_data_frame(g_list$edges, directed = FALSE, 
                                  vertices = vertices_df %>%
                                    dplyr::select(node_names, 
                                           dplyr::everything())
  ) 
  
  return(g_out)
  
}
