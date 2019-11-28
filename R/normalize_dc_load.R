#' Normalise load
#' 
#' This function normalises the load across a DC power network.
#' 
#' @param g An igraph object. Thr graph representation of the network to be normalised
#' @param demand  A character string. The node attribute containing the flow demand
#' @param generation A character string. The node attribute containing the flow generation
#' @param net_generation A character string. The node atribute containing the difference between demand and generation for 
#' the current parametrization
#' @param edge_name A character string. The edge attribbute containing the edge names
#' @param capacity A character string. The edge attribute containing the flow capacity of each edge
#' @param node_name A character string. The node attibute containing the node names
#' @param power_flow A character string. The name of the edge atribute that contains the flow data
#' @export

normalise_dc_load <- function(g, 
                           generation = "generation", 
                           demand  = "demand",
                           net_generation = "net_generation", 
                           capacity = "edge_limit",
                           edge_name = "edge_name", 
                           node_name = "node_name",
                           power_flow = "power_flow"
                           ){
  #This function normalised the power on the nodes and edges so that the power on the nodes so that the positive generation on the nodes sums to 1.
  #The edge capacities are then normalised so that the alpha values are the same as pre normalisation
  #This process only works for DC power-flow

  #calculate power flow to make sure flow values are correct
  SlackRef <- SlackRefFunc(g, name = node_name,
                           Generation = net_generation) #find the most appropriate node to be the slack bus
   print(SlackRef)
  g2 <- PowerFlow(g, SlackRef$name, 
                  EdgeName = edge_name, 
                  VertexName = node_name,
                  Net_generation = net_generation,
                  power_flow = power_flow)

  #find total positive power
  pos_pow <- get.vertex.attribute(g2, net_generation) %>%
    ifelse(.>0, ., 0) %>%
    sum()

  #normalise node load generation and net_generation relative to total positive power-flow
  #this produces a data frame that contains normalised node data
  #Bang bang ans ":=" are ued to assign the character values represented by the variable names to be applied
  g_node_df <- as_data_frame(g2, what = "vertices") %>%
    mutate(!!generation := .data[[generation]]/pos_pow,
           !!demand := .data[[demand]]/pos_pow,
           !!net_generation := .data[[net_generation]]/pos_pow)

  #given the normalised values re-calculate power flow
  g3 <- g2 %>%
    set.vertex.attribute(., net_generation, value = g_node_df[,net_generation]) %>%
    PowerFlow(., SlackRef$name, 
              EdgeName = edge_name, 
              VertexName = node_name,
              Net_generation = net_generation,
              power_flow = power_flow)
  
  #get the normalised power flow vector
  power_flow_vect <- get.edge.attribute(g3, power_flow)
  
  #This creates a data frame of normalised edge data. That is edge capacity and power flow
  g_edge_df <- as_data_frame(g2, what = "edges") %>%
    mutate(alpha_temp = abs(.data[[capacity]]/.data[[power_flow]]), #calculate the original line alpha values
           !!power_flow := power_flow_vect, #add in the new power flow
           !!capacity := abs(.data[[power_flow]]*alpha_temp)) %>%   #re-calculate edge capacity using the original alpha values and new power flow
    select(-alpha_temp) #remove the alpha values from the data frame

#re-make the network using the normalised data.
  g_out <- graph_from_data_frame(g_edge_df, vertices = g_node_df, directed = FALSE) %>%
    BalencedGenDem(., demand, 
                   generation,
                   OutputVar  = net_generation)

  #return normalised network
  return(g_out)
}
