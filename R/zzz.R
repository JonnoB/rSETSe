#This file contains the global variables.
#These variables are usually values that are created by functions and then referenced by other functions
#This leads R to think there might be an error when the R CMD checks are being performed.
#By adding these variables to the global functions this becomes clear there is no error
utils::globalVariables(c("from", 
                         "to", 
                         "node",
                         "Iter",
                         "elevation",
                         "elevation.x",
                         "elevation.y",
                         "velocity",
                         "static_force", 
                         "friction",
                         "net_force",
                         "net_tension",
                         "temp",
                         "res_stat",
                         "de",
                         "H",
                         "mean_e",
                         "tension",
                         "strain",
                         "position.x",
                         "position.y",
                         "rows", #in setse_data_prep
                         "cols", #in setse_data_prep
                         "A", #the area variable in calc_spring_area
                         "value_2", #calc_spring_area
                         "edge_name", #create_node_edge_df
                         "."))
