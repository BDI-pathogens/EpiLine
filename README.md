# EpiLine - Estimating epi-curves and distributions from case line list data

This package contains models for estimating epi-curves and individual infection progression distributions simulataneously. 
The estimators require both total case data by date as well as detailed ("line list") information for a subset of individual.

## Symptom-Report Model
This models the delay between the onset of symptoms and the presenting/testing/reporting of cases to health authorities. 
Early in outbreaks these delays can are frequently large due to lack of awareness of the symptoms, however, with increased awareness due to public health information the delays will decrease in time. 
Even if these delays are constant with time, when estimating their distribution it is necessary to consider both censoring effects and the underlying dynamics of the infection to avoid biases. 
Conversely, when estimating the dynamics of then infection (e.g. r(t) or R(t)) from reported cases, it is necessary to know these distributions. 
Therefore, if both infection dynamics and reporting delays are varying at the same time, the best way to account for the biases is to simulataneosly estimate both.

### Model Decrtiption
The aim of the model is to understand the interaction between the symptom-report time distribution and the underlying dynamics of the infection rate, therefore we use a very simple for the number of people developing the symptoms each day. 
We model the daily growth rate $r(t)$ with a Gaussian process, so the daily number of people of developing symptoms $S(t)$ is given by

$$
\begin{align}
  r(t) &= r( t - 1 ) + \epsilon(t), \\
  S(t) &= S(t-1) e^{r(t)}, \\
  \epsilon(t) & \sim N(0,\sigma^2_{r_{GP}}),
\end{align}
$$

where $\sigma^2_{r_{GP}}$ is the daily variance of the Gaussian process. 
Note that by making $r(t)$ a Gaussian process instead of $S(t)$ directly a Gaussian process, it means that the prior is that the expected daily chainge in $S(t)$ is the same as the previous day. 
Next we define $f(\tau,t)$, which is probability of someome reporting an infection on day $t$ if they developed symptoms on day $t-\tau$. 
Note that $\tau$ can be negative if a case is found prior to symptoms developing (e.g. it contact-traced and tested positive).
On day $t$
the expected number of case reported is $\mu(t)$ and given by

$$
\mu(t) = \sum_{\tau = -\tau_{\rm post}}^{\tau_{\rm pre}} f(\tau,t) S(t-\tau)
$$

where $\tau_{\rm pre}$ is the maximum number of days pre-reporting the case develops symptoms and 
$\tau_{\rm pre}$ the maximumn number of days post-reporting the case develops symptoms.


### Example Results


## Installation
This is a R package and during the package build the Stan code is compiled. To build this package, clone the repository and then call `R CMD INSTALL --no-multiarch --with-keep.source $FULL_PATH_REPO_DIR`, where `$FULL_PATH_REPO_DIR` is the full path to the directory where the respository was cloned to. The package require `rstan`, `Rcpp`, `rstantools`, `StanHeaders`, `data.table`, `moments`, `R6`, `matrixStats` and `plotly` to be installed (all avaialbe on CRAN).

