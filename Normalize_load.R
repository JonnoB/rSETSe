Normalize_load <- function(g, 
                           Generation = Generation_MW, 
                           Demand  = Load_MW,
                           Net_Generation = Net_Generation, 
                           capacity = Link.Limit,
                           EdgeName =Link, 
                           VertexName = name){
  #This function normalised the power on the nodes and edges so that the power on the nodes so that the positive generation on the nodes sums to 1.
  #The edge capacities are then normalised so that the alpha values are the same as pre normalisation
  #This process only works for DC power-flow
  
  Net_gen <- enquo(Net_Generation) %>% as_label()


  #calculate power flow to make sure flow values are correct
  SlackRef <- SlackRefFunc(g, name = as_label(enquo(VertexName)),
                           Generation = Net_gen) #find the most appropriate node to be the slack bus
   print(SlackRef)
  g2 <- PowerFlow(g, SlackRef$name, EdgeName = as_label(enquo(EdgeName)), VertexName = as_label(enquo(VertexName)),
                  Net_generation = Net_gen)

  #find total positive power
  pos_pow <- get.vertex.attribute(g2, Net_gen) %>%
    ifelse(.>0, ., 0) %>%
    sum()

  #normalise node load and re-calculate power-flow
  g_node_df <- as_data_frame(g2, what = "vertices") %>%
    mutate({{Generation}} := {{Generation}}/pos_pow,
           {{Demand}} := {{Demand}}/pos_pow,
           {{Net_Generation}} := {{Net_Generation}}/pos_pow)

  #add new edge powerflow to edge dataframe
  g3 <- g2 %>%
    set.vertex.attribute(., Net_gen, value = g_node_df[,Net_gen]) %>%
    PowerFlow(., SlackRef$name, EdgeName = as_label(enquo(EdgeName)), VertexName = as_label(enquo(VertexName)),
              Net_generation = Net_gen)

  g_edge_df <- as_data_frame(g2, what = "edges") %>%
    #calculate edge alpha values i.e. temp_alpha
    mutate(alpha_temp = ({{capacity}}/PowerFlow) %>% abs,
           PowerFlow = get.edge.attribute(g3, "PowerFlow"),
           {{capacity}} := abs(PowerFlow*alpha_temp)) %>%   #calculate new edge limits
    select(-alpha_temp)

  #add new edge limits back into normalised network
  # g_out <- g2 %>%
  #   set.edge.attribute(., as_label(enquo(capacity)), value = pull(g_edge_df, {{capacity}}))

  g_out <- graph_from_data_frame(g_edge_df, vertices = g_node_df, directed = FALSE) %>%
    BalencedGenDem(., as_label(enquo(Demand)), as_label(enquo(Generation)),
                   OutputVar  = Net_gen)

  #return normalised network
  return(g_out)
}