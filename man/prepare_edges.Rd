% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_edges.R
\name{prepare_edges}
\alias{prepare_edges}
\title{Prepare network edges}
\usage{
prepare_edges(g, k = NULL, distance = 1, create_edge_name = TRUE)
}
\arguments{
\item{g}{an igraph object}

\item{k}{The spring constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
the k value will not be added to the network. This is useful k is made through some other process.}

\item{distance}{The spring constant. This value is either a numeric value giving the spring constant for all edges or NULL. If NULL is used 
the distance value will not be added to the network. This is useful distance is made through some other process.}

\item{create_edge_name}{Logical. Whether to create and edge name attribute or not.
 
 @details 
 The function prepares the edge characteristics of the network so that they can be embedded using the SETSe_ family of functions.
 
 @return 
 The function creates several variables
\itemize{
  \item k: The spring constant representing the stiffness of the spring. 
  \item distance: The minimum distance between nodes. This is the distance between the parallel planes/hyper-planes.
  \item edge_name: the name of the edges. it takes the form "from_to" where "from" is the origin node and "to" is the destination node using the 
 \link[igraph]{as_data_frame} function from igraph
}}
}
\description{
This function helps prepare the network edges for embedding
}
\examples{
set.seed(234) #set the random see for generating the network
g <- generate_peels_network(type = "E")
embeddings <- g \%>\%
prepare_edges(k = 500, distance = 1) \%>\%
#prepare the network for a binary embedding
prepare_categorical_force(., node_names = "name",
                     force_var = "class") \%>\%
#embed the network using auto setse
  setse_auto(., force = "class_A")
}
\seealso{
\link{setse}, \link{setse_auto}, \link{setse_bicomp}, \link{setse_auto_hd}

Other prepare_setse: 
\code{\link{prepare_categorical_force}()},
\code{\link{prepare_continuous_force}()}
}
\concept{prepare_setse}
