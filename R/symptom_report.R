##################################################################
#  Name: symptom_report.simulator
#
#  Description: Runs a simulation of the symptom_report model and
#  returns the line list of cases by day
# 
#  Arguments:
#  t_rep          - time periods for which data is reported
#  t_symptom_pre  - time before report period to run the simulation
#  t_symptom_post - time after report period to run the simulation
#  symptom_0      - initial number of symptomatic individuals
#  r              - daily r(t) at start of simulation (single value or vector of length t_rep, or t_rep+t_symptom_pre+t_symptom_post)
#  dist_xi        - symptom-report distribution xi parameter (Johnson SU)
#  dist_lambda    - symptom-report distribution lambda parameter (Johnson SU)
#  dist_gamma     - symptom-report distribution gamma parameter (Johnson SU)
#  dist_delta     - symptom-report distribution delta parameter (Johnson SU)
#  report_var_over_mean - the variance over mean for negative binomial distribution of daily cases given expected cases
###################################################################
symptom_report.simulator <- function(
  t_rep          = 55,
  t_symptom_pre  = 20,
  t_symptom_post = 5,
  symptom_0      = 2,
  r              = 0.05,
  dist_xi        = 5,
  dist_lambda    = 3,
  dist_gamma     = -2,
  dist_delta     = 1.5,
  report_var_over_mean = 2
) 
{
  # run the simulation of symptomatic cases for an extended period of time
  # so we can estimate the number of reported cases in the required range
  t_rep_min <- t_symptom_pre
  t_rep_max <- t_rep - t_symptom_post
  t_max     <- t_rep + t_symptom_pre + t_symptom_post 
  if( t_rep < 1 )
    stop( "t_rep must be positive")
  if( t_symptom_pre < 10 )
    stop( "t_symptom_pre must be at least 10")
  if( t_symptom_post < 5 )
    stop( "t_symptom_post must be at least 5")
  
  # check r is acceptable
  if( length( r ) == 1 ) {
    r <- rep( r, t_max )
  } else if( length( r ) == t_rep ) {
    r <- c( rep( r[1], t_symptom_pre ), r, rep( r[1], t_symptom_pre[ t_rep ] ) )
  } else if( length( r ) != t_max )
    stop( "r must be of length 1 or t_rep or t_max")

  # daily symptomic cases is given by 
  symptom <- round( symptom_0 * exp( cumsum( r ) ) )
  
  # check the length of the distribution paramters
  if( length( dist_xi ) == 1 )     dist_xi     <- rep( dist_xi, t_rep )
  if( length( dist_lambda ) == 1 ) dist_lambda <- rep( dist_lambda, t_rep )
  if( length( dist_gamma ) == 1 )  dist_gamma  <- rep( dist_gamma, t_rep )
  if( length( dist_delta ) == 1 )  dist_delta  <- rep( dist_delta, t_rep )
  if( length( dist_xi ) != t_rep | length( dist_lambda ) != t_rep | length( dist_gamma ) != t_rep | length( dist_delta ) != t_rep  )
    stop( "dist_xxxxx paramters must be of length 1 or t_rep" )
  
  # get the expected number reported each day
  report <- rep(0, t_rep)
  for( t in 1:t_rep )
    for( t_off in -4:20 ) {
      tt        <- t + t_rep_min - t_off
      report[t] <- report[t] + symptom[tt] * .djsu( t_off, dist_xi[t], dist_lambda[t], dist_gamma[t], dist_delta[t])
    }
  
  # the number of actual report cases is assumed to be a negative binomial based on this
  if( report_var_over_mean == 1 ) {
    report <- rpois( t_rep, report );
  } else if( report_var_over_mean > 1 ) {
    report <- rnbinom( t_rep, size = report / ( report_var_over_mean - 1 ), mu = report )
  } else 
    stop( "variance over mean must be greater or eqaual to one")
 
  # generate a global linelist
  ll_report   <- rep( 1:t_rep, report )
  ll_gamma    <- rep( dist_gamma, report )
  ll_delta    <- rep( dist_delta, report )
  ll_xi       <- rep( dist_xi, report )
  ll_lambda   <- rep( dist_lambda, report )
  ll_symptom  <- ll_report - round( .rjsu( length( ll_report), ll_xi, ll_lambda, ll_gamma, ll_delta ) )
  linelist    <- data.table( report = ll_report, symptom = ll_symptom )
  
  return( list( symptom = symptom, report = report, linelist = linelist ) )
}