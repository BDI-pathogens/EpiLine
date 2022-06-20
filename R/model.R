##################################################################
#  Name: model.symptom_report.stan 
#
#  Description: Gets the Stan model object for the symptom_report model
# 
###################################################################
model.symptom_report.stan <- function() {
  return( stanmodels$symptom_report_model )
}