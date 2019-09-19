Calc_line_strain <- function(g, solved_height_df, distance, capacity, flow = "PowerFlow"){
  #g the network graph
  # solved_height_df a data frame with the heights of each node
  #distance the unquoted name of the distance variable in graph g
#
# I am doing this strange code choice as I have problems passing the names through multiple functions
#
#
  
  
  
  
  
  # Out <- as_data_frame(g) %>% as_tibble %>%
  #   left_join(., solved_height_df %>% select(node, z), by = c("from"= "node")) %>%
  #   left_join(., solved_height_df %>% select(node, z), by = c("to"= "node")) %>%
  #   mutate(dz = abs(z.x-z.y),
  #          mean_z = (z.x+z.y)/2,
  #          H = sqrt(dz^2 +{{distance}}^2),
  #          strain = (H-{{distance}})/{{distance}},
  #         alpha = {{capacity}}/abs({{flow}}),
  #          line_load = abs({{flow}})/{{capacity}},
  #          percentile_strain = percent_rank(strain)) %>%
  #   select(Link, alpha, line_load, dz, H, strain, percentile_strain, mean_z, {{flow}})
  
  distance <- sym(distance)
  capacity <- sym(capacity)
  flow <- sym(flow)
  
  Out <- as_data_frame(g) %>% as_tibble %>%
    left_join(., solved_height_df %>% select(node, z), by = c("from"= "node")) %>%
    left_join(., solved_height_df %>% select(node, z), by = c("to"= "node")) %>%
    mutate(dz = abs(z.x-z.y),
           mean_z = (z.x+z.y)/2,
           H = sqrt(dz^2 +(!!distance)^2),
           strain = (H-!!distance)/!!distance,
           alpha = !!capacity/abs(!!flow),
           line_load = abs(!!flow)/!!capacity,
           percentile_strain = percent_rank(strain)) %>%
    select(Link, alpha, line_load, dz, H, strain, percentile_strain, mean_z, {{flow}})
  
  return(Out)
  
}