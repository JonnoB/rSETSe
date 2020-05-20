#' Embedd categorical values
#' 
#' Performs SETSe on every level of a categorical variable. This function is experimental
#' 
#' @param g
#' @param force_var
#' @param node_names
#' @param k
#' 
#' @export
#' 
categorical_SETSe <- function(g, force_var, node_names, k = "k"){
  
  
  #do an embedding on each level, taking out only the node elevation and the tension and strain.
  #all other data is discarded
  all_levels <- unique(vertex_attr(g, force_var)) %>%
    map(~{
      print(.x)
      embeddings_data <- g %>%
        prepare_SETSe_binary(., node_names = node_names, k = NULL, force_var = force_var, positive_value = .x) %>%
        auto_SETSe()
      
      Out <-list(node_embeddings = embeddings_data$node_embeddings %>%
                   select(node, elevation)%>% 
                   mutate(level = .x),
                 edge_embeddings = embeddings_data$edge_embeddings %>%
                   select(edge_name, mean_e, tension, strain) %>% 
                   mutate(level = .x)
      )
      
    }) %>%
    transpose() %>%
    map(~bind_rows(.x))
  
  node_embeddings <- all_levels$node_embeddings %>%
    group_by(level) %>%
    mutate(elevation_adjust = elevation + min(elevation)) %>%
    group_by(node) %>%
    summarise(
      mean_elevation = mean(elevation),
      euc_elevation = sqrt(sum(elevation_adjust^2)))
  
  dist_node <- all_levels$node_embeddings %>%
    pivot_wider(., names_from = level, values_from = elevation)
  
  edge_embeddings <- all_levels$edge_embeddings %>%
    group_by(edge_name) %>%
    summarise(
      mean_tension = mean(tension),
      euc_tension = sqrt(sum(tension^2)),
      mean_strain = mean(strain),
      euc_strain = sqrt(sum(strain^2)))
  
  node_details <- edge_embeddings %>%
    separate(., col = edge_name, into = c("from", "to")) %>%
    pivot_longer(cols = from:to, , values_to = "node") %>%
    select(-name) %>%
    group_by(node) %>%
    summarise_all(mean) %>%
    left_join(., node_embeddings, by = "node") %>%
    left_join(as_data_frame(g, what = "vertices"), by = c("node" ="name")) %>%
    arrange(as.numeric(node))
  
  out <- list( node_embeddings = node_embeddings,
               edge_embeddings = edge_embeddings,
               node_details = node_details,
               dist_node = dist_node)
  
  return(out)
  
}
