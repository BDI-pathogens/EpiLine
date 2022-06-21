# 
# Example:     flat_d_rist
# Description: Simple example where both r(t) and the symptom-report time distribution are constant
#

library( EpiLine )

# run the simulation using default paramters
simulation <- symptom_report.simulator()

# data 
reported   <- simulation$report
ll_report  <- simulation$linelist_report
ll_symptom <- simulation$linelist_symptom

fit <- symptom_report.fit( reported, ll_symptom, ll_report )
 