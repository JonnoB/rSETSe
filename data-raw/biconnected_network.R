## code to prepare `two_bicomponents` dataset goes here
library(igraph)
library(dplyr)
biconnected_network <- bind_rows(expand.grid(LETTERS[1:4],LETTERS[1:4]),
                              expand.grid(LETTERS[5:7],LETTERS[5:7]),
                              tibble(Var1 = "D", Var2 = "E")) %>%
  filter(Var1 != Var2) %>%
  mutate(edge_name = paste(Var1, Var2, sep ="_"),
         weight = case_when(
           Var1 %in% LETTERS[1:4] & Var2 %in% LETTERS[1:4] ~ 1000,
           Var1 %in% LETTERS[5:7] & Var2 %in% LETTERS[5:7] ~ 500,
           TRUE ~ 100
         )) %>%
  graph_from_data_frame(., directed = F) %>%
  set_vertex_attr(., name = "force", value = (betweenness(.) -mean(betweenness(.)))) %>%
  set.vertex.attribute(., "group", value = ifelse(get.vertex.attribute(., "name") %in% LETTERS[1:4], "A", "B"))


usethis::use_data(biconnected_network, overwrite = TRUE)
