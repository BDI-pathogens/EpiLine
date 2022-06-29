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
  r(t) &\sim N( r( t - 1 ), \sigma^2_{r_{GP}}) \\
  S(t) &= S(t-1) e^{r(t)}, 
\end{align}
$$

where $\sigma^2_{r_{GP}}$ is the daily variance of the Gaussian process. 
Note that by making $r(t)$ a Gaussian process instead of $S(t)$ directly a Gaussian process, it means that the prior is that the expected daily chainge in $S(t)$ is the same as the previous day. 
Next we define $f(\tau,t)$, which is probability of someome reporting an infection on day $t+\tau$ if they developed symptoms on day $t$. 
Note that $\tau$ can be negative if a case is found prior to symptoms developing (e.g. it contact-traced and tested positive).
On day $t$ the expected number of case reported is $\mu(t)$ and given by

$$
\mu(t) = \sum_{\tau = -\tau_{\rm post}}^{\tau_{\rm pre}} f(\tau,t-\tau) S(t-\tau)
$$

where $\tau_{\rm pre}$ is the maximum number of days pre-reporting the case develops symptoms and 
$\tau_{\rm pre}$ the maximumn number of days post-reporting the case develops symptoms.
The number of observed reported cases $C(t)$ is modelled as negative binomial variable

$$
C(t) \sim NB(\mu(t),\phi_{OD}),
$$

where $\phi_{OD}$ is the over-dispersion parameter.

The symptom-report time distribution must support both postive and negative values. In addition, empirically it is observed that this distribution can be highly skewed with heavy tails, therefore we model it using the Johnson SU distribution which contains 4 parameters $(\xi, \lambda, \gamma,\delta)$.
To account for the changes in the distribution over time, we model these 4 parameters using Gaussian processes

$$
\begin{align}
  \xi(t) &\sim N( \xi( t - 1 ), \sigma^2_{\xi_{GP}}), \\
  \lambda(t) &\sim N( \lambda( t - 1 ), \sigma^2_{\lambda_{GP}}), \\
  \gamma(t) &\sim N( \gamma( t - 1 ), \sigma^2_{\gamma_{GP}}), \\
  \delta(t) &\sim N( \delta( t - 1 ), \sigma^2_{\delta_{GP}}.) 
\end{align}
$$

These parameters are estimated using line list data of individual cases where the symptoms date report date are known.

Range priors are put on the initial values of all the parameters and the variances of the Gaussian processes. 
For efficiency, we allow for the period between sampling of the Gaussian processes to be greater than day, in which case we linearly interpolate for the intermediate days.
The posterior of the model is sampled using Stan and is contained within this R package.
The function `symptom_report.fit` fits the model to data and `symptom_report.simulator` generates simulated data under the model.

### Example Results

We now demonstrate the model using simulated data which is contained in `examples/linear_r_dist.R`.
In this example, the daily growth rate $r(t)$ declines linearly throughout the simulation from 0.1 to -0.03, and the mean and variance of the symptom-report distribution decrease linearly with time (with constand skewness and kurtosis). 

```
library( EpiLine )

# define the length of the simulatiopn
t_rep          <- 50 # length of time for which data is reported
t_symptom_pre  <- 30 # time before the reporting period to simulate
t_symptom_post <- 5  # time after the reporting period to simulate
t_max          <- t_rep + t_symptom_post + t_symptom_pre

# set up the varaible r(t) and distribution
symptom_0 <- 20                               # initial number of symptomatic people
r         <- 0.1 - 0.13 * ( 1:t_max ) / t_max # r(t) inthe simulation
xi        <- -1 + 6 * ( t_max:1 ) / t_max          # xi parameter in the symptom-report dist
lambda    <- 2 + ( t_max:1 ) / t_max         # lambda parameter in the symptom-report dist

simulation <- symptom_report.simulator(
  t_rep          = t_rep,
  t_symptom_pre  = t_symptom_pre,
  t_symptom_post = t_symptom_post,
  symptom_0   = symptom_0,
  r           = r,
  dist_xi     = xi,
  dist_lambda = lambda
)
```

We next sample from the posterior of the model using Stan (note we're using very short chains for speed so mean and variance of parameter estimates will be of limited accuracy).

```
# data 
reported    <- simulation$reported
ll_report   <- simulation$linelist$report
ll_symptom  <- simulation$linelist$symptom
report_date <- as.Date("2022-04-01")

# fit using model
mcmc_n_samples <- 100
mcmc_n_chains  <- 1
fit <- symptom_report.fit( reported, ll_symptom, ll_report, report_date = report_date, 
                           mcmc_n_samples = mcmc_n_samples, mcmc_n_chains = mcmc_n_chains )
```

Once the fit is complete (this examples takes roughly 20s), we can plot the posterior of the fitted parameters against the simulation parameters.
First we consider the estimate of the number of people developing symptoms on each day (`fit$plot.symptoms( simulation = simulation`).

<img src="https://github.com/BDI-pathogens/EpiLine/blob/main/documentation/linear_symptoms.png" width="700" >

The dotted vertical lines in the chart show the period over which case reports were provided. 
The extended time pre- and post- the reporting window is required because some of the reported cases in this period will have developed symptoms outside of the window,.
Next we plot the estimated $r(t)$ for the entire period (including pre- and post- the reporting window).

<img src="https://github.com/BDI-pathogens/EpiLine/blob/main/documentation/linear_r.png" width="700" >

Note that outside the reportiny window the estimated $r(t)$ flattens out and is different from the simulated one. 
This is because the vast majority of the symptomtic case in these periods will not be contained within the reported data, therefore the prior distribution on $r(t)$ will not be 

Finally we plot the estimated symptom-report distribution at the start and end of the reporting period.

<img src="https://github.com/BDI-pathogens/EpiLine/blob/main/documentation/linear_distribution.png" width="700" >

Note that model successfully estimates the distribution at the start and end of the reporting period.

## Usage on Real Data
To run the model on real data you need to provide the function `symptom_report.fit()` with both the time-series of the number of reported cases and the line list of symptom-report case pairs. Note that it is not necessary to have line list data for all reported cases, so for a reported case where there is no data about the time of symptoms, it should be included in the time-series but not the line list. 

Required arguments of `symptom_report.fit( reported, linelist_symptom, linelist_report )`
1. `reported` - a vector of integers with length of the reporting period containing the total number of cases reported each day. The first entry is for the earliest report date with following entries being ordered for each day of the reporting period (if no reported cases on a day then it should be entered as a 0).
2. `linelist_symptom` - a vector with length of the number of cases for which both symptom and report date are known. The value is an integer of the date that each case developed symptoms relative to the reporting period start date (a value of 1 means symptoms were developed on the first day for which there are reported cases). Note that values are allowed to be positive.
3. `linelist_report` - a vector with length of the number of cases for which both symptom and report date are known, paired with `linelist_symptom`. The value is an integer of the date that each was reported relative to the reporting period start date. Note that all reported cases should be during the reporting periods used in `reported`, i.e. all values are between 1 and the length of `reported`.

The output charts are member functions of the R6 fit object returned by `fit=symptom_report.fit()`. The output is shown in the Example section above for results on simulated data (note with real data you do not provide a simulation object i.e. just call `fit$plot.symptoms()`). 

## Installation
This is a R package and during the package build the Stan code is compiled. To build this package, clone the repository and then call `R CMD INSTALL --no-multiarch --with-keep.source $FULL_PATH_REPO_DIR`, where `$FULL_PATH_REPO_DIR` is the full path to the directory where the respository was cloned to. The package require `rstan`, `Rcpp`, `rstantools`, `StanHeaders`, `data.table`, `moments`, `R6`, `matrixStats` and `plotly` to be installed (all avaialbe on CRAN).

