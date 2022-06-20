# EpiLine - Estimating epi-curves and distirbutions from case line list data


## Installation
This is a R package and during the package build the Stan code is compiled. To build this package, clone the repository and then call `R CMD INSTALL --no-multiarch --with-keep.source $FULL_PATH_REPO_DIR`, where `$FULL_PATH_REPO_DIR` is the full path to the directory where the respository was cloned to. The package require `rstan`, `Rcpp`, `rstantools`, `StanHeaders`, `data.table` and `plotly` to be installed (all avaialbe on CRAN).

