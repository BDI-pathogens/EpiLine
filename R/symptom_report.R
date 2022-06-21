.eps = 1e-6

##################################################################
#  Name: symptom_report.simulator
#
#  Description: Runs a simulation of the symptom_report model and
#  returns the line list of cases by day
# 
#  Arguments:
#  t_rep          - time periods for which data is reported
#  t_symptom_pre  - maximum time before report of onset of symptoms
#  t_symptom_post - maximum time after report of onset of symptoms
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
  dist_xi        = 2,
  dist_lambda    = 2.5,
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

##################################################################
#  Name: symptom_report.fit
#
#  Description: Estimates the parameters in the model from the daily
#  total number of reported cases and the line list of symptom-report pairs
# 
#  Arguments:
#  reported         - vector of number of reported cases by day
#  linelist_symptom - vector of symptom dates (1 per person)
#  linelist_report  - vector of symptom dates (1 per person)
#  t_symptom_pre    - maximum time before report of onset of symptoms
#  t_symptom_post   - maximum time after report of onset of symptoms
###################################################################
symptom_report.fit <- function(
  reported,
  linelist_symptom,
  linelist_report,
  t_symptom_pre  = 20,
  t_symptom_post = 5,
  mcmc_n_samples = 1e2,
  mcmc_n_chains  = 1,
  prior_gamma_min  = -20,
  prior_delta_min  = -20,
  prior_xi_min     = -20,
  prior_lambda_min = 0,
  prior_gamma_max  = 20,
  prior_delta_max  = 20,
  prior_xi_max     = 20,
  prior_lambda_max = 20,
  prior_xi_gp_sd_max      = 0.6,
  prior_lambda_gp_sd_max  = 0.4,
  prior_gamma_gp_sd_max   = 0.4,
  prior_delta_gp_sd_max   = 0.2,
  prior_r_0_min           = -0.1,
  prior_r_0_max           = 0.2,
  prior_r_gp_sd_max       = 0.055,
  prior_phi_od_max        = 1,  
  prior_log_symptoms0_min = log( 0.5 ),
  prior_log_symptoms0_max = log( 50 ),
  hyper_gp_period_r    = 2,
  hyper_gp_period_dist = 2
)
{
    # calculate the times require for the simultion
    t_rep <- length( reported )
    t_max <- t_rep + t_symptom_post + t_symptom_pre
  
    # calculate totals and symptom-report pairs
    if( length( linelist_symptom ) != length( linelist_report ) )
      stop( "linelist of symptom and report dates must be equal" )
    linelist <- data.table( report = linelist_report, symptom = linelist_symptom )
    
    # enforce the max/min times
    linelist[ , symptom := pmin( pmax( symptom, report - t_symptom_pre ), report + t_symptom_post ) ]
    linelist <- linelist[ , .N, by = c("report", "symptom" ) ][ order(report,symptom) ]
    
    # shift the dates by t_symptom_pre since the fitted GPs run for the extra time
    linelist[ , report  := report + t_symptom_pre ]
    linelist[ , symptom := symptom + t_symptom_pre ]
    
    # prepare the data for Stan
    data <- list(
      t_max     = t_max,
      t_rep_min = t_symptom_pre + 1,
      t_rep_max          = t_rep + t_symptom_pre,
      t_rep_symptoms_max = t_symptom_pre - 1,
      t_rep_symptoms_min = -t_symptom_post+1,
      reported    = reported,
      n_ll        = linelist[, .N ],
      ll_report   = linelist[ , report ],
      ll_symptoms = linelist[ , symptom ],
      ll_N        = linelist[ , N], 
      prior_gamma_min  = prior_gamma_min ,
      prior_delta_min  = prior_delta_min,
      prior_xi_min     = prior_xi_min,
      prior_lambda_min = prior_lambda_min,
      prior_gamma_max  = prior_gamma_max,
      prior_delta_max  = prior_delta_max,
      prior_xi_max     = prior_xi_max,
      prior_lambda_max = prior_lambda_max,
      prior_xi_gp_sd_max      = prior_xi_gp_sd_max,
      prior_lambda_gp_sd_max  = prior_lambda_gp_sd_max,
      prior_gamma_gp_sd_max   = prior_gamma_gp_sd_max,
      prior_delta_gp_sd_max   = prior_delta_gp_sd_max,
      prior_log_symptoms0_min = prior_log_symptoms0_min,
      prior_log_symptoms0_max = prior_log_symptoms0_max,
      prior_r_0_min        = prior_r_0_min,
      prior_r_0_max        = prior_r_0_max,
      prior_r_gp_sd_max    = prior_r_gp_sd_max,
      prior_phi_od_max     = prior_phi_od_max,
      hyper_gp_period_r    = hyper_gp_period_r,
      hyper_gp_period_dist = hyper_gp_period_dist
    )
    
    # initialise some of the data in the chain  
    # 1. moment match the Johnson SU parameters for the entire line list
    ll_mean <- mean( linelist[ , rep( report - symptom, N)])
    ll_var  <- var( linelist[ , rep( report - symptom, N)])
    ll_skew <- skewness( linelist[ , rep( report - symptom, N)])
    
    delta  <- runif( 1,1,2)
    t_skew <- function( gamma ) return( .jsu.skewness( 1, 1, gamma, delta ) - ll_skew )
    gamma  <- uniroot( t_skew, lower = -100, upper = 100)$root
    gamma  <- pmin( pmax( gamma, prior_gamma_min + .eps ), prior_gamma_max - .eps )
    
    t_var  <- function( lambda ) return( .jsu.var( 1, lambda, gamma, delta ) - ll_var )
    lambda <- uniroot( t_var, lower = 1e-6, upper = 100)$root
    lambda <- pmin( pmax( lambda, prior_lambda_min + .eps ), prior_lambda_max - .eps )
    
    t_mean <- function( xi )return( .jsu.mean( xi, lambda, gamma, delta ) - ll_mean )
    xi     <- uniroot( t_mean, lower = -100, upper = 100)$root
    xi     <- pmin( pmax( xi, prior_xi_min + .eps ), prior_xi_max - .eps )
    
    # 2.initialize r(0)=0 and for the GPs to be constant  
    init_func <- function(x) {
      return( list(
        r_gp = rep( 0, ceil( data$t_max / data$hyper_gp_period_r ) ),
        r_0 = 0,
        gamma0 = gamma,
        delta0 = delta,
        lambda0 = lambda,
        xi0 = xi,
        xi_gp = rep( 0, ceil( t_rep / data$hyper_gp_period_dist ) ),
        lambda_gp = rep( 0, ceil( t_rep / data$hyper_gp_period_dist ) ),
        gamma_gp = rep( 0, ceil( t_rep / data$hyper_gp_period_dist ) ),
        delta_gp = rep( 0, ceil( t_rep / data$hyper_gp_period_dist ) )
      ) ) } 
    
    # get Stan model and sample
    model <- model.symptom_report.stan()
    init  <- lapply( 1:mcmc_n_chains, init_func )
    raw   <-sampling(
      model,
      chains =mcmc_n_chains,
      cores = min( mcmc_n_chains, 3 ),
      iter = mcmc_n_samples,
      data = data,
      init = init,
      pars = c( "symptoms", "r", "delta", "gamma", "lambda", "xi", "r_gp_sd")
    )
    return( raw )
}