##################################################################
#  Name: model.symptom_report.stan 
#
#  Description: Gets the Stan model object for the symptom_report model
# 
###################################################################
model.symptom_report.stan <- function( distribution = distribution.johnson_su ) {
  if( distribution == distribution.johnson_su )
    return( stanmodels$symptom_report_model )
  if( distribution == distribution.gamma )
    return( stanmodels$symptom_report_model_gamma )
  
  stop( sprintf( "Distribution %s not implemented", distribution ) )
}