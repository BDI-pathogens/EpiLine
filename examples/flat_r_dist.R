# 
# Example:     flat_d_rist
# Description: Simple example where both r(t) and the symptom-report time distribution are constant
#

library( EpiLine )

# run the simulation using default paramters
simulation <- symptom_report.simulator()

# data 
reported    <- simulation$reported
ll_report   <- simulation$linelist$report
ll_symptom  <- simulation$linelist$symptom
report_date <- as.Date("2022-04-01")

fit <- symptom_report.fit( reported, ll_symptom, ll_report, report_date = report_date )
fit$plot.symptoms( simulation = simulation ) 
fit$plot.r( simulation = simulation )
fit$plot.symptom_report.dist( simulation = simulation )

