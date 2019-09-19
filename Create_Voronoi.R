Create_Voronoi <- function(node_values, sfc_map){
  #node values a dataframe containing a Node columns a Longitude column a Latitude column as well as any desired attribute columns
  #sfc_map and sfc object describing the area that the nodes are in.
  
  #convert the dataframe to an sf object
  node_z_sf <- node_values %>%
    gather(key = type, value = value, -Node,-Longitude, -Latitude ) %>%
    group_by(type) %>%
    mutate(value = percent_rank(value)) %>% #get the percent rank of the variables
    ungroup %>%
    st_as_sf(., coords = c("Longitude", "Latitude")) 
  #add in the correct projection data to align with the map
  st_crs(node_z_sf) <- st_crs(sfc_map)
  node_z_sf <- st_transform(node_z_sf, st_crs(sfc_map))
  
  #create the tessellation
  v <- st_voronoi(st_combine(node_z_sf)) 
  
  #clip the tesselation then convert back into an sf object to allow the simple features to be joined back in
  #the key simple feature is of course z
  clipped_tess_node <- st_intersection(st_cast(v), st_union(sfc_map)) %>% 
    st_sf(geom=.) %>%
    st_join(node_z_sf) %>%
    filter(!is.na(type))
  
}
