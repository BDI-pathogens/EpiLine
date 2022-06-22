library( EpiLine )
library( testthat )
library( moments )

# set seed to prevent stochastic failings
set.seed( 1 )

# check the sample value of the Johnson SU distribution agree the analytic
# formulas for the distribution moments
n_tests <- 5
samples <- 1e6
tol     <- 0.05
xi      <- runif( n_tests, min = -5, max = 10 )
lambda  <- runif( n_tests, min = 0.2, max = 5 )
gamma   <- runif( n_tests, min = -5, max = 5 )
delta   <- runif( n_tests, min = 0.2, max = 4 )


for( idx in 1:n_tests ) {
  sample   <- .rjsu( samples, xi[idx], lambda[idx], gamma[idx], delta[idx] )
  mean     <- .jsu.mean( xi[idx], lambda[idx], gamma[idx], delta[idx] )
  var      <- .jsu.var( xi[idx], lambda[idx], gamma[idx], delta[idx] )
  skewness <- .jsu.skewness( xi[idx], lambda[idx], gamma[idx], delta[idx] )

  test_that( "Estimated JSU sample mean",     expect_lt( abs( mean( sample ) - mean ), tol ) )
  test_that( "Estimated JSU sample variance", expect_lt( abs( var( sample ) - var ), tol ) )
  test_that( "Estimated JSU sample skewness", expect_lt( abs( skewness( sample ) - skewness ), tol ) )
}