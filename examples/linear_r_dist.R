# 
# Example:     linear_d_rist
# Description: Example where both r(t) and the symptom-report time distribution 
#              both change linearly with time
#

library( EpiLine )
library( plotly )

set.seed( 1 )

# define the length of the simulation
t_rep          <- 50 # length of time for which data is reported
t_symptom_pre  <- 30 # time before the reporting period to simulate
t_symptom_post <- 5  # time after the reporting period to simulate
t_max          <- t_rep + t_symptom_post + t_symptom_pre

# set up the varaible r(t) and distribution
symptom_0 <- 2                                # initial number of symptomatic people
r         <- 0.1 - 0.13 * ( 1:t_max ) / t_max # r(t) inthe simulation
xi        <- -1 + 6 *( t_max:1 ) / t_max      # xi parameter in the symptom-report dist
lambda    <- 2 + ( t_max:1 ) / t_max          # lambda parameter in the symptom-report dist

simulation <- symptom_report.simulator(
  t_rep          = t_rep,
  t_symptom_pre  = t_symptom_pre,
  t_symptom_post = t_symptom_post,
  symptom_0   = symptom_0,
  r           = r,
  dist_xi     = xi,
  dist_lambda = lambda
)

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

# save plots
plt = fit$plot.symptoms( simulation = simulation )
save_image(plt, "test_linear_symptoms.png")
plt = fit$plot.r( simulation = simulation )
save_image(plt, "test_linear_r.png")
plt = fit$plot.symptom_report.dist( simulation = simulation )
save_image(plt, "test_linear_distribution.png")
plt = fit$plot.symptom_report.quantiles( simulation = simulation, quantiles = c( 0.1,0.5,0.9) )
save_image(plt, "test_linear_distribution_percentiles.png")

fit$plot

# save data
# write.csv()