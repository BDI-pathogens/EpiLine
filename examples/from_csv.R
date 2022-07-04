# 
# Example:     from_csv.R
# Description: Example of a fit where the data is read directly from csv
#
#   file_reported  - is a csv containing 2 columns, the first column the 
#                    date of report and second the number of cases reported 
#                    on that day (first row header: date,reported)
#
#   file_linelist - is a csv containing 2 columns, the first column is the date
#                   a case was reported and the second column is the date on 
#                   which the case had symptoms (first row header: report,symptom)
#
# 

library( EpiLine )

file_reported <- "examples/linear_r_reported.csv"
file_linelist <- "examples/linear_r_linelist.csv"

# fit using model
mcmc_n_samples <- 100
mcmc_n_chains  <- 1
fit <- symptom_report.fit( 
  file_reported  = file_reported,
  file_linelist  = file_linelist,
  mcmc_n_samples = mcmc_n_samples, 
  mcmc_n_chains  = mcmc_n_chains 
)

fit$plot.symptoms() 
fit$plot.r()
fit$plot.symptom_report.dist()
fit$plot.symptom_report.quantiles( quantiles = c( 0.1,0.5,0.9) )



