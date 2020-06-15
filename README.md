
<img src="https://github.com/JonnoB/rSETSe/blob/master/SETSe_logo.png" width="100">


# rSETSe
## An R package for embedding graphs using the SETSe algorith,

This is the R package for the Strain Elevation Tension Spring embeddgins (SETSe) algorithm. SETSe is a deterministic graph embeddings algorithm. It converts the node attributes of a graph into forces and the edge attributes into springs. The algorithm finds an equilibrium position when the forces of the nodes are balanced by the forces on the springs. A full descpription of the algorithm is given in my forthcoming paper "The spring bounces back: Introduction to Strain Elevation Tension Spring embedding for network representation"

# Installation instructions

 1. Open R/Rstudio and ensure that devtools has been installed
 1. Run the following code library(devtools); install_github("JonnoB/rSETSE")
 1. Load the Power Grid Networking package normally using library(rSETSe)
 1. All functions have help files e.g ?auto_SETSe

The package can also be downloaded or cloned then installed locally using the install function from devtools.

# Basic use

```
library(rSETSe)

#prepares a graph for embedding using SETSe
g_setse <- g %>%
  remove_small_components(.) %>%
  prepare_SETSe_continuous(., node_names = "name", k = 1000,
                           sum_to_one = FALSE)
  
#Embedds using the bi-connected auto-parametrization algorithm.
#This method is strongly reccomended, it tends to be much faster and almost always converges
SETSe_bicomp(g,_setse,
             tol = sum(abs(vertex_attr(g_setse, "force")))/1000,
             hyper_tol = 0.1,
             hyper_iters = 3000,
             verbose = T)
```

The code used in this package is part of my PhD and a fundemental component of three projects

* The spring bounces back ([code](https://github.com/JonnoB/SETSe_assortativity_and_clusters))
* Power grid robusteness ([code](https://github.com/JonnoB/setse_and_network_robustness/edit/master/README.md))

# N.B

The package is still underdevleopment and does not have any examples of code. For current examples of code see the links above.
