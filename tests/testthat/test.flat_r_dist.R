library( EpiLine )
library( testthat )
library( moments )

# set seed to prevent stochastic failings
set.seed( 1 )

n_tests   <- 2
tol_sd    <- 3
r         <- runif( n_tests, min = -0.1, max = 0.1 )
xi        <- runif( n_tests, min = 0, max = 5 )
lambda    <- runif( n_tests, min = 1, max = 2 )
gamma     <- runif( n_tests, min = -3, max = -1 )
delta     <- runif( n_tests, min = 1.5, max = 2.5 )
symptom_0 <- 20
t_rep     <- 20
t_symptom_post <- 5
t_symptom_pre  <- 30
t_max          <- t_rep + t_symptom_post + t_symptom_pre

for( idx in 1:n_tests ) 
{
  # simulate with the test parameters
  set.seed( 1 )
  simulation  <- symptom_report.simulator( 
    t_rep       = t_rep,
    symptom_0   = symptom_0,
    r           = r[idx],
    dist_xi     = xi[idx],
    dist_lambda = lambda[idx],
    dist_gamma  = gamma[idx],
    dist_delta  = delta[idx],
    t_symptom_post = t_symptom_post,
    t_symptom_pre  = t_symptom_pre 
  )
  
  # quick fit with a small number of samples (so won't be fully converged)
  set.seed( 1 )
  suppressWarnings( invisible( capture.output( fit <- symptom_report.fit( 
    simulation$reported, 
    simulation$linelist$symptom, 
    simulation$linelist$report,
    t_symptom_post = t_symptom_post,
    t_symptom_pre  = t_symptom_pre, 
    mcmc_n_chains  = 1,
    mcmc_n_samples = 1e2
  ) ) ) )
  
  # check the estimation of r(t) at 2 start and end of simulation
  r_1   <- fit$stan_extract$r[ , t_symptom_post ]
  r_max <- fit$stan_extract$r[ , t_symptom_post + t_rep ]
  test_that( "Posterior r(1) mean", expect_lt( abs( mean( r_1 ) - r[idx] ) / sd( r_1 ), tol_sd ) )
  test_that( "Posterior r(t_rep) mean", expect_lt( abs( mean( r_max ) - r[idx] ) / sd( r_max ), tol_sd ) )
  
  # check the moments of the fitted distribution
  xi_1     <- fit$stan_extract$xi[ , t_symptom_pre ]
  lambda_1 <- fit$stan_extract$lambda[ , t_symptom_pre ]
  gamma_1  <- fit$stan_extract$gamma[ , t_symptom_pre ]
  delta_1  <- fit$stan_extract$delta[ , t_symptom_pre ]
  mean_1   <- .jsu.mean( xi_1, lambda_1, gamma_1, delta_1)
  var_1    <- .jsu.var( xi_1, lambda_1, gamma_1, delta_1)
  mean_1_s <- .jsu.mean( xi[idx], lambda[idx], gamma[idx], delta[idx])
  var_1_s  <- .jsu.var( xi[idx], lambda[idx], gamma[idx], delta[idx])

  test_that( "Posterior mean of symptom-report dist", expect_lt( abs( mean( mean_1 ) - mean_1_s ) / sd( mean_1 ), tol_sd ) )
  test_that( "Posterior varaince of symptom-report dist", expect_lt( abs( mean( var_1 ) - var_1_s ) / sd( var_1 ), tol_sd ) )
  
  xi_1     <- fit$stan_extract$xi[ , t_symptom_pre + t_rep  ]
  lambda_1 <- fit$stan_extract$lambda[ , t_symptom_pre + t_rep ]
  gamma_1  <- fit$stan_extract$gamma[ , t_symptom_pre + t_rep ]
  delta_1  <- fit$stan_extract$delta[ , t_symptom_pre + t_rep ]
  mean_1   <- .jsu.mean( xi_1, lambda_1, gamma_1, delta_1)
  var_1    <- .jsu.var( xi_1, lambda_1, gamma_1, delta_1)
  mean_1_s <- .jsu.mean( xi[idx], lambda[idx], gamma[idx], delta[idx])
  var_1_s  <- .jsu.var( xi[idx], lambda[idx], gamma[idx], delta[idx])
  
  test_that( "Posterior mean of symptom-report dist", expect_lt( abs( mean( mean_1 ) - mean_1_s ) / sd( mean_1 ), tol_sd ) )
  test_that( "Posterior varaince of symptom-report dist", expect_lt( abs( mean( var_1 ) - var_1_s ) / sd( var_1 ), tol_sd ) )
}

  
