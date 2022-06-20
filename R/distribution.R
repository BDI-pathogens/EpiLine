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
  return( xi + lambda * sinh( ( runif(n) -gamma) / delta))
