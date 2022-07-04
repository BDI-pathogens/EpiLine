##################################################################
#  Name: .djsu
#
#  Description: PDF of the Johnson SU distribtuion
##################################################################
.djsu <- function( x, xi, lambda, gamma, delta ) {
  z  <- ( x - xi) / lambda
  z2 <- gamma + delta * asinh( z )
  return( delta/lambda/sqrt(2*pi)/sqrt(1+z^2)*exp(-0.5*z2^2))
}

##################################################################
#  Name: .rjsu
#
#  Description: random variable of Johnson SU distribtuion
##################################################################
.rjsu <- function( n, xi, lambda, gamma, delta ) 
  return( xi + lambda * sinh( ( rnorm(n) -gamma) / delta))

##################################################################
#  Name: .qjsu
#
#  Description: random variable of Johnson SU distribtuion
##################################################################
.qjsu <- function( q, xi, lambda, gamma, delta ) 
  return( xi + lambda * sinh( ( qnorm(q) -gamma) / delta))

##################################################################
#  Name: .jsu.mean
#
#  Description: mean of the Johnson SU distribution
##################################################################
.jsu.mean <- function( xi, lambda, gamma, delta ) {
  omega <- exp(1/delta^2);
  return( xi - lambda *sqrt(omega) * sinh( gamma / delta ) )
}

##################################################################
#  Name: .jsu.var
#
#  Description: variance of the Johnson SU distribution
##################################################################
.jsu.var <- function( xi, lambda, gamma, delta ) {
  omega <- exp(1/delta^2);
  return( 0.5 * (lambda^2 )* (omega - 1)*(omega*cosh(2*gamma/delta) + 1 ) )
}

##################################################################
#  Name: .jsu.skewness
#
#  Description: variance of the Johnson SU distribution
##################################################################
.jsu.skewness <- function( xi, lambda, gamma, delta ) {
  omega <- exp(1/delta^2);
  skew  <- -0.25*(lambda^3)*sqrt(omega)*((omega-1)^2)*(omega*(omega+2)*sinh(3*gamma/delta)+3*sinh(gamma/delta))
  return( skew / .jsu.var( xi, lambda, gamma, delta )^{3/2})
}

