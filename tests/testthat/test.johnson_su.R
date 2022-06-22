library( EpiLine )
library( testthat )
library( moments )

# set seed to prevent stochastic failings
set.seed( 1 )

# check the sample value of the Johnson SU distribution agree the analytic
# formulas for the distribution moments
n_tests      <- 5
samples      <- 1e6
tol_sample   <- 0.05
tol_analytic <- 1e-4
xi      <- runif( n_tests, min = -5, max = 10 )
lambda  <- runif( n_tests, min = 0.2, max = 5 )
gamma   <- runif( n_tests, min = -5, max = 5 )
delta   <- runif( n_tests, min = 0.2, max = 4 )

for( idx in 1:n_tests ) {
  # check the sample value of the Johnson SU distribution agree the analytic
  # formulas for themoments
  
  sample   <- .rjsu( samples, xi[idx], lambda[idx], gamma[idx], delta[idx] )
  mean     <- .jsu.mean( xi[idx], lambda[idx], gamma[idx], delta[idx] )
  var      <- .jsu.var( xi[idx], lambda[idx], gamma[idx], delta[idx] )
  skewness <- .jsu.skewness( xi[idx], lambda[idx], gamma[idx], delta[idx] )

  test_that( "Estimated JSU sample mean",     expect_lt( abs( mean( sample ) - mean ), tol_sample ) )
  test_that( "Estimated JSU sample variance", expect_lt( abs( var( sample ) - var ), tol_sample ) )
  test_that( "Estimated JSU sample skewness", expect_lt( abs( skewness( sample ) - skewness ), tol_sample ) )

  # check integral of Johnson SU pdf afree with the formulas for the moments
  integrand_1 <- function( x ) return( .djsu( x, xi[idx], lambda[idx], gamma[idx], delta[idx]))
  int_1       <- integrate( integrand_1, lower = -100, upper = 100)$value
  integrand_m <- function( x ) return( x * .djsu( x, xi[idx], lambda[idx], gamma[idx], delta[idx]))
  int_m       <- integrate( integrand_m, lower = -100, upper = 100)$value
  integrand_v <- function( x ) return( x^2 * .djsu( x, xi[idx], lambda[idx], gamma[idx], delta[idx]))
  int_v       <- integrate( integrand_v, lower = -100, upper = 100)$value - int_m^2
  integrand_s <- function( x ) return( x^3 * .djsu( x, xi[idx], lambda[idx], gamma[idx], delta[idx]))
  int_s       <- integrate( integrand_s, lower = -100, upper = 100)$value 
  int_s       <- ( int_s - 3 * int_v * int_m - int_m^3 ) / int_v^1.5
  
  test_that( "Integral of pdf is 1", expect_lt( abs( int_1 - 1 ), tol_analytic ) )
  test_that( "Mean from pdf is analytic formula", expect_lt( abs( int_m / mean - 1), tol_analytic ) )
  test_that( "Variance from pdf is analytic formula", expect_lt( abs( int_v / var - 1), tol_analytic ) )
  test_that( "Skewness from pdf is analytic formula", expect_lt( abs( int_s / skewness - 1), tol_analytic * 20 ) )
  
}

