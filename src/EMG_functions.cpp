#include <Rcpp.h>  // Includes certain libraries of functions
#include <math.h>
#include <limits>
#include "miscfunctions.h" // Linear interpolation

/*


*/

// Lookup - 01
// Scalar function to generate random draws from convolution of
// exponential and gaussian random variables

double remg_scl( double mu, double sigma, double lambda ) {

  // Initialize output
  double out = NA_REAL;

  // Check for inadmissable values
  if ( ( sigma > 0 ) && ( lambda > 0 ) ) {

    // Convert lambda to rate
    lambda = 1.0/lambda;

    // Sum of independent draws for a gaussian and an
    // exponential variable
    out = R::rnorm( mu, sigma ) + R::rexp( lambda );
  }

  return( out );
}

// Lookup - 02
//' The Exponentially Modified Gaussian Distribution
//'
//' Random generation, density, distribution, and quantile functions for
//' the distribution of the convolution of gaussian and exponential
//' random varibles.
//'
//' @param N the number of draws for random generation.
//' @param x a vector of quantiles.
//' @param mu a vector of the means for the gaussian variable.
//' @param sigma a vector of standard deviations for the gaussian variable
//'   (sigma > 0).
//' @param lambda a vector of rates for the exponential variable
//'   (lambda > 0).
//' @param ln indicates whether the log-likelihood should be returned,
//'   where 1 = True, 0 = False (the default).
//'
//' @section Details:
//' An exponentially modified gaussian distribution describes the sum of
//' independent normal and exponential random variables, possessing
//' a characteristic positive skew due to the exponential variable.
//'
//' For unequal vector lengths, values are recycled.
//'
//' @section References:
//' Forthcoming
//'
//' @examples
//' Forthcoming
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector remg(int N, Rcpp::NumericVector mu,
                              Rcpp::NumericVector sigma,
                              Rcpp::NumericVector lambda) {

  int N_mu = mu.size(); // Number of parameters
  int N_lambda = lambda.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int mu_inc = 0;
  int lambda_inc = 0;
  int sigma_inc = 0;

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector mu_v(N);
  Rcpp::NumericVector lambda_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    mu_v(nv) = mu(mu_inc);
    lambda_v(nv) = lambda(lambda_inc);
    sigma_v(nv) = sigma(sigma_inc);

    mu_inc = mu_inc + 1;
    lambda_inc = lambda_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_mu==mu_inc) mu_inc = 0;
    if (N_lambda==lambda_inc) lambda_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Generate draws
  for (int n = 0; n < N; n++) {
    out(n) = remg_scl( mu_v(n), sigma_v(n), lambda_v(n) );
  }

  return( out );
}

// Lookup - 04
// Scalar function for density of the convolution of an exponential
// and gaussian random variable

double demg_scl( double x, double mu, double sigma,
                 double lambda, int ln = 0 ) {

  // Initialize output
  double out = 0.0;

  // Check for inadmissable variables
  if ( (lambda > 0.0) &&
       (sigma > 0.0) ) {

    double sigma2 = pow( sigma, 2.0 );
    double p1 = lambda/2.0;
    double p2 = 2.0*mu + lambda*sigma2 - 2.0*x;
    double p3 = ( mu + lambda*sigma2 - x )/(pow(2.0,.5)*sigma);

    out = p1*exp(p1*p2)*erfc(p3);

  }

  if ( out < 0.0 ) out = 0.0;
  if ( ln == 1 ) out = log( out );

  return( out );
}

// Lookup - 05
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector demg(Rcpp::NumericVector x,
                         Rcpp::NumericVector mu,
                         Rcpp::NumericVector sigma,
                         Rcpp::NumericVector lambda,
                         int ln = 0 ) {

  int N_x = x.size(); // Number of observations
  int N_mu = mu.size(); // Number of parameters
  int N_lambda = lambda.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int x_inc = 0;
  int mu_inc = 0;
  int lambda_inc = 0;
  int sigma_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_x, N_mu,
                                           N_lambda, N_sigma) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector x_v(N);
  Rcpp::NumericVector mu_v(N);
  Rcpp::NumericVector lambda_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    x_v(nv) = x(x_inc);
    mu_v(nv) = mu(mu_inc);
    lambda_v(nv) = lambda(lambda_inc);
    sigma_v(nv) = sigma(sigma_inc);

    x_inc = x_inc + 1;
    mu_inc = mu_inc + 1;
    lambda_inc = lambda_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_x==x_inc) x_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;
    if (N_lambda==lambda_inc) lambda_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Determine density
  for (int n = 0; n < N; n++) {
    out(n) = demg_scl( x_v(n), mu_v(n), sigma_v(n), lambda_v(n), ln );
  }

  return( out );
}

// Lookup - 06
// Scalar function for CDF of the convolution of an exponential
// and gaussian random variable

double pemg_scl( double x, double mu, double sigma,
                 double lambda ) {

  // Initialize output
  double out = 0.0;

  // Check for inadmissable variables
  if ( (lambda > 0.0) &&
       (sigma > 0.0) ) {

    double u = lambda*(x - mu);
    double v = lambda*sigma;

    double p1 = R::pnorm( u, 0.0, v, 1, 0 );
    double p2 = -u + pow(v,2.0)/2.0;
    double p3 = R::pnorm( u, pow(v,2.0), v, 1, 1 );

    out = p1 - exp( p2 + p3);

  }

  if ( (out < 0.0) || (out > 1.0) ) out = 0.0;

  return( out );
}

// Lookup - 07
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pemg(Rcpp::NumericVector x,
                         Rcpp::NumericVector mu,
                         Rcpp::NumericVector sigma,
                         Rcpp::NumericVector lambda ) {

  int N_x = x.size(); // Number of observations
  int N_mu = mu.size(); // Number of parameters
  int N_lambda = lambda.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int x_inc = 0;
  int mu_inc = 0;
  int lambda_inc = 0;
  int sigma_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_x, N_mu,
                                           N_lambda, N_sigma) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector x_v(N);
  Rcpp::NumericVector mu_v(N);
  Rcpp::NumericVector lambda_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    x_v(nv) = x(x_inc);
    mu_v(nv) = mu(mu_inc);
    lambda_v(nv) = lambda(lambda_inc);
    sigma_v(nv) = sigma(sigma_inc);

    x_inc = x_inc + 1;
    mu_inc = mu_inc + 1;
    lambda_inc = lambda_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_x==x_inc) x_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;
    if (N_lambda==lambda_inc) lambda_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Determine density
  for (int n = 0; n < N; n++) {
    out(n) = pemg_scl( x_v(n), mu_v(n), sigma_v(n), lambda_v(n) );
  }

  return( out );
}


/*
library(seqmodels)

# Define generating parameters
prm = c( mu = .5, sigma = .5, lambda = 4 )

# Simulate values
N = 10000 # Number of draws
# Base R
sim.test = rnorm(N,prm[1],prm[2]) + rexp(N,prm[3])
# C++
sim = remg(N,prm[1],prm[2],prm[3])
x11()
hist( sim, breaks=40, col='grey', border='white', freq=F, main = ' ' )
den = density(sim.test)
lines( den$x, den$y, lty = 2, lwd = 2 )

# Add density
x = seq( -4, 6, length = 1000 )
lines( x, demg( x, prm[1], prm[2], prm[3] ), col ='orange', lwd = 2 )

# CDF
p = 1:length(sim)
p = p/length(sim)
x11();
plot( sort(sim), p, type = 'l', xlab = 'x', ylab = 'CDF', bty = 'l' )
lines( x, pemg(x,prm[1],prm[2],prm[3]), col='orange', lwd = 2 )


*/
