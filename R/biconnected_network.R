#' A simple network made of three bi-connected components
#'
#' The data set can be used to explore simple different embeddings methods on a very simple graph
#'
#' @format An igraph network with 7 nodes and 19 edges which forms three biconnected components:
#' \describe{
#'   \item{edge_name}{The name of the edge connecting the two vertices}
#'   \item{weight}{The edge weight connecting the two vertices. This value is 1000 for edges connecting nodes
#'   A to D, it is 500 for edges connecting nodes E to G, it is 100 connecting nodes D and E}
#'   \item{force}{The force produced by each node. It was calculated by subtracting the mean node centrality for the network
#'   from the node centrality}
#'   \item{group}{The group each node is in. This can be used to generate force if required}
#' }
#' @examples
#' \dontrun{plot(biconnected_network)}
"biconnected_network"
