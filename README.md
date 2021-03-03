# rSETSe <img src='man/figures/SETSe_logo.png' align="right" height="300" />


<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/rsetse)](https://cran.r-project.org/package=rsetse)
[![R build status](https://github.com/JonnoB/rSETSe/workflows/R-CMD-check/badge.svg)](https://github.com/JonnoB/rSETSe/actions)
[![Travis build status](https://travis-ci.com/JonnoB/rSETSe.svg?branch=master)](https://travis-ci.com/JonnoB/rSETSe)
<!-- badges: end -->


## An R package for embedding graphs using the SETSe algorithm

This is the R package for the Strain Elevation Tension Spring embeddings (SETSe) algorithm. SETSe is a deterministic graph embeddings algorithm. It converts the node attributes of a graph into forces and the edge attributes into springs. The algorithm finds an equilibrium position when the forces of the nodes are balanced by the forces on the springs. A full description of the algorithm is given in "The spring bounces back: Introduction to Strain Elevation Tension Spring embedding for network representation" ([Bourne 2020](https://doi.org/10.1007/s41109-020-00329-4)). There is a website for the package  providing documentation and vignettes at https://jonnob.github.io/rSETSe/index.html . This is a very niche package so please feel free to reach out to me on twitter or through email with questions.

# Installation instructions

The package is available on CRAN and can be installed by running `install.packages("rsetse")`.Alternatively it can be installed from github using the below method.

 1. Open R/Rstudio and ensure that devtools has been installed
 1. Run the following code library(devtools); install_github("JonnoB/rSETSe")
 1. Load the package normally using library(rsetse)
 1. All functions have help files e.g ?SETSe_auto

The package can also be downloaded or cloned then installed locally using the install function from devtools.

# Basic use

```
library(rSETSe)

#prepares a graph for embedding using SETSe
set.seed(234) #set the random see for generating the network
g <- generate_peels_network(type = "E") %>%
  #prepare the network for a binary embedding
  prepare_SETSe_binary(., node_names = "name", k = 1000, 
                       force_var = "class", 
                       positive_value = "A") 
                       
#Embedds using the bi-connected auto-parametrization algorithm.
#This method is strongly reccomended, it tends to be much faster and almost always converges
embeddings <- SETSe_bicomp(g,
                           tol = sum(abs(vertex_attr(g, "force")))/1000,
                           hyper_tol = 0.1,
                           hyper_iters = 3000,
                           verbose = T)

```

# Cite

To cite rsetse in publications use: Bourne, J. The spring bounces back: introducing the strain elevation tension spring embedding algorithm for network representation. Appl Netw Sci 5, 88 (2020). https://doi.org/10.1007/s41109-020-00329-4
