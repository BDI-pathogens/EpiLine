# EpiLine - Estimating epi-curves and distributions from case line list data

This package contains models for estimating epi-curves and individual infection progression distributions simulataneously. The estimators require both total case data by date as well as detailed ("line list") information for a subset of individual.

## Symptom-Report Model
This models the delay between the onset of symptoms and the presenting/testing/reporting of cases to health authorities. Early in outbreaks these delays can are frequently large due to lack of awareness of the symptoms, however, with increased awareness due to public health information the delays will decrease in time. Even if these delays are constant with time, when estimating their distribution it is necessary to consider both censoring effects and the underlying dynamics of the infection to avoid biases. Conversely, when estimating the dynamics of then infection (e.g. r(t) or R(t)) from reported cases, it is necessary to know these distributions. Therefore, if both infection dynamics and reporting delays are varying at the same time, the best way to account for the biases is to simulataneosly estimate both.


## Installation
This is a R package and during the package build the Stan code is compiled. To build this package, clone the repository and then call `R CMD INSTALL --no-multiarch --with-keep.source $FULL_PATH_REPO_DIR`, where `$FULL_PATH_REPO_DIR` is the full path to the directory where the respository was cloned to. The package require `rstan`, `Rcpp`, `rstantools`, `StanHeaders`, `data.table`, `moments`, `R6`, `matrixStats` and `plotly` to be installed (all avaialbe on CRAN).

