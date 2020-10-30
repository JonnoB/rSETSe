#' Create balanced blocks
#' 
#' Separates the network into a series of bi-connected components that can be solved separately. 
#' Solving smaller subgraphs using the bi-connected component method reduces the risk of network divergence.
#' This function is seldom called independently of SETSe_bicomp
#' 
#' @details When networks are separated into the bi-connected subgraphs or blocks. The overall network balance needs to be maintained. 
#' \code{create_balanced_blocks} maintains the balance by summing the net force across the all the nodes that are being removed from
#' the subgraph. Therefore a node that is an articulation point has a force value equal to the total of all the nodes on the adjacent
#' bi-connected component.
#' 
#' @param g An igraph object. The network for which embeddings will be found
#' @param force A character vector. The name of the node attribute that is the force exerted by the nodes
#' @param bigraph A list. the list of biconnected components produced by the biconnected_components function.
#'  This function take a non trivial amount of time on large graphs so this pass through minimises the function being called.
#' @return A list containing all the bi connected component where each component is balanced to have a net force of 0.
#' 
#' @examples 
#' \dontrun{
#' #create a list of balanced network using the biconnected_network dataset
#' balanced_list <-create_balanced_blocks(biconnected_network, 
#' bigraph = biconnected_components(biconnected_network))
#' 
#' #count the edges in each of the bi-components
#' sapply(balanced_list,ecount )
#' }
#' @export

create_balanced_blocks <- function(g, force = "force", bigraph = bigraph){
  #This function creates a list of biconnected components or blocks.
  #These blocks are balanced such that the connecting vertices contain all the power of the missing part of the network
  #balancing prevents the network reaching a steady state non-zero velocity.
  
  ArticulationPoints <- igraph::get.vertex.attribute(g, "name", bigraph$articulation_points) #can also use names(biconnected.components(g)$articulation_points)
  
  biconnected_component_raw <- igraph::as_data_frame(g, what = "vertices")
  edge_df_raw <- igraph::as_data_frame(g)
  
  List_of_BiConComps <-1:length(bigraph$components) %>%
    purrr::map(~{
      Comp_num <- .x
      
      #nodes in the current block
      Nodes_in_j <- igraph::get.vertex.attribute(g, "name", bigraph$components[[Comp_num]]) 
      
      #The nodes in the current block that are articulation points. 
      #There will be at least one articulation point, unless the whol network is one block
      #In a two node block that does not terminate at an end both nodes are articulation points.
      component_art_points <- Nodes_in_j[Nodes_in_j %in% ArticulationPoints]
      
      biconnected_component <- biconnected_component_raw[biconnected_component_raw$name  %in% Nodes_in_j,]
      
      balanced_component_df <- 1:length(component_art_points) %>%
        purrr::map_df(~{
          
          CurrArt <- component_art_points[.x]
          
          #delete the articulation point to split the network into 2 peices
          membership_vect <- igraph::components(delete.vertices(g, CurrArt))$membership
          
          
          #Identify the component that contains nodes from the same block as the articulation point. 
          #This component needs to be discarded.
          discard_components <- unique(membership_vect[names(membership_vect) %in% Nodes_in_j]) 
          
          #The remaining nodes are those for which we want to know the total force
          nodes_of_interest <- c(names(membership_vect[membership_vect != discard_components]), CurrArt)
          
          data.frame(name = CurrArt,
                     AuxPower = sum(biconnected_component_raw[(biconnected_component_raw$name %in% nodes_of_interest), force])) 
          
        }) %>%
        dplyr::left_join(biconnected_component, ., by = "name") %>%
        dplyr::mutate(temp = ifelse(is.na(AuxPower), !!rlang::sym(force), AuxPower)) 
      
      
      Component_j <- {!(igraph::get.vertex.attribute(g, "name") %in% balanced_component_df$name)} %>%
        (1:igraph::vcount(g))[.] %>%
        igraph::delete.vertices(g,.) 
      
      balanced_component <- data.frame(name = igraph::get.vertex.attribute(Component_j, "name")) %>%
        dplyr::left_join(balanced_component_df, by = "name") %>%
        dplyr::pull(temp) %>%
        igraph::set.vertex.attribute(Component_j, name = force, value = .)
      
      return(balanced_component)
    })
  
}
