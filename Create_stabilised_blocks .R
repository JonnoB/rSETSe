Create_stabilised_blocks <- function(g, OriginBlock, OriginBlock_number, 
                                     force ="net_generation", 
                                     flow = "power_flow", 
                                     capacity = "edge_limit", 
                                     edge_name = "edge_name",  
                                     tstep=0.1, 
                                     tol = 1e-10, 
                                     distance, 
                                     maxIter, 
                                     mass, 
                                     verbose = TRUE){
  #This function finds the z displacement of all the nodes in the network.
  #g and igraph object of the network
  #OriginBlock a dataframe, output of the Find_network_Balance function of the graph, should be stable or almost stable. mass and K of this block
    #should be used as the basis for the rest of the inputs
  #OriginBlock_number the order number of the origin block in the after creating balanced blocks using Create_balanced_blocks
  #tstep time interval for each iteration
  #tol the tolerance of the system for early termination
  #distance the name of the edge atribute that is used to define the horizontal distance between nodes, this must be meaningful for the units or 1
  #maxIter the maxmimum number of Iterations if no eaerly termination possible
  #kbase the minimum k value must be the same as used in OriginBlock
  #kdiff the range of k values, must be the same as used in OriginBlock
  #mass the mass of the nodes, must be the same as the origin block
  
  #Seperate out the graph into balanced blocks
  #This step will have already been done, but it is fast and simplifies the requirements for the function
  List_of_BiConComps <- Create_balanced_blocks(g, force = force, flow = flow)
  
  #remove the Origin block so it doesn't have to be calculated again
  BlockNumbers <-(1:length(List_of_BiConComps))[-OriginBlock_number]
  
  StabilModels <- BlockNumbers %>% 
    map(~{

      Out <- Find_network_balance(List_of_BiConComps[[.x]],
                                  force =force, 
                                  flow = flow, 
                                  capacity = capacity,  
                                  tstep = tstep, 
                                  tol = tol, 
                                  distance = distance, 
                                  edge_name = edge_name,
                                  maxIter =  maxIter, 
                                  mass =  mass, 
                                  verbose = FALSE)
      
      #print if the print requirement is on otherwise silent
      if(!verbose){print(paste("Block" ,.x, "of", max(BlockNumbers) ,"termination", nrow(Out$results) )) }
      
      return(Out)
      
    })
  
  #get the block tree of the graph
  Block_tree <- biconnected_components(g)
  
  #extract the articulation nodes
  ArticulationVect <- get.vertex.attribute(g, "name", Block_tree$articulation_points)
  
  #place all nodes relative to the origin
  relative_blocks1 <- 1:length(StabilModels) %>% 
    map_df(~{
      print(.x)
      StabilModels[[.x]]$NodeList %>%
        mutate(Reference_ID = .x)
      
    }) %>%
    bind_rows(OriginBlock$NodeList %>% 
                mutate(Reference_ID = 0)) %>%
    mutate(Articulation_node = (node %in% ArticulationVect ))
  
  
  #The height of each node relative to the origin and normalised
  final_z <-fix_z_to_origin(relative_blocks1, ArticulationVect) %>%
    group_by(node) %>%
    summarise_all(mean) %>%
    mutate(Articulation_node = Articulation_node==1)
  return(final_z)
  
}