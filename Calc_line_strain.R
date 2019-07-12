Calc_line_strain <- function(g, solved_height_df, distance){
  
  line_strain <-as_data_frame(g) %>% as_tibble %>%
    left_join(., solved_height_df %>% select(node, z), by = c("from"= "node")) %>%
    left_join(., solved_height_df %>% select(node, z), by = c("to"= "node")) %>%
    mutate(dz = abs(z.x-z.y),
           mean_z = (z.x+z.y)/2,
           H = sqrt(dz^2 +{{distance}}^2),
           strain = (H-{{distance}})/{{distance}},
           alpha = Link.Limit/abs(PowerFlow),
           line_load = abs(PowerFlow)/Link.Limit,
           percentile_strain = percent_rank(strain)) %>%
    select(Link, alpha, line_load, dz, H, strain, percentile_strain, mean_z, PowerFlow)
  
}