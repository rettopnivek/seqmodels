#include <Rcpp.h> // Includes certain libraries of functions
#include <math.h>
#include <limits>
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Assorted functions for the calculation of the density, distribution,
quantile, and random number generation functions of the exponentially
modified gaussian distribution.

References:
Forthcoming

Index
Lookup - 01:  remg_scl
Lookup - 02:  remg
Lookup - 03:  demg_scl
Lookup - 04:  demg
Lookup - 05:  pemg_scl
Lookup - 06:  pemg
Lookup - 07:  qemg_scl
Lookup - 08:  qemg

### TO DO ###
Add references for the exponentially modified gaussian
Add examples
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

// Lookup - 03
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

// Lookup - 04
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

// Lookup - 05
// Scalar function for CDF of the convolution of an exponential
// and gaussian random variable

double pemg_scl( double x, double mu, double sigma,
                 double lambda ) {

  // Initialize output
  double out = 0.0;

  // Check for inadmissable variables
  if ( (lambda > 0.0) &&
       (sigma > 0.0) &&
       (x != R_NegInf) ) {

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

// Lookup - 06
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

  // Determine CDF
  for (int n = 0; n < N; n++) {
    out(n) = pemg_scl( x_v(n), mu_v(n), sigma_v(n), lambda_v(n) );
  }

  return( out );
}

// Lookup - 07
// A scalar function that calculates the quantile given a cumulative probability
// using linear interpolation.

double qemg_scl( double p, double mu, double sigma, double lambda,
                 double mnRT, double mxRT, int em_stop, double err ) {

  double cur_t = mnRT; // Initialize output

  // Define an initial set of RTs
  std::vector<double> iRT(5);
  for (int i = 0; i < 5; i++) {
    iRT[i] = (i/5.0)*(mxRT-mnRT)+mnRT;
  }

  // Determine the associated CDF values
  std::vector<double> iPrb(5);
  for (int i = 0; i < 5; i++) {
    iPrb[i] = pemg_scl( iRT[i], mu, sigma, lambda );
  }

  // Determine the initial window that the point falls between
  int p0 = minMax( p, iPrb, 0 );
  int p1 = minMax( p, iPrb, 1 );

  double lbp = iPrb[p0]; double lbt = iRT[p0];
  double ubp = iPrb[p1]; double ubt = iRT[p1];

  double prev_t = ubt; double prev_prb = ubp;
  cur_t = linInterp( p, lbp, ubp, lbt, ubt );
  double prb = pemg_scl( cur_t, mu, sigma, lambda );

  double epsilon = 1.0;
  int cnt = 0;

  while ( (cnt < em_stop) & (epsilon > err) ) {
    if (prb < p) {
      lbp = prb;
      lbt = cur_t;
    } else if ( lbp < prb ) {
      ubp = prb;
      ubt = cur_t;
    } else {
      lbp = prb;
      lbt = cur_t;
      ubp = prev_prb;
      ubt = prev_t;
    }
    prev_t = cur_t; prev_prb = prb;
    cur_t = linInterp( p, lbp, ubp, lbt, ubt );
    prb = pemg_scl( cur_t, mu, sigma, lambda );

    cnt = cnt + 1;
    epsilon = std::abs( p - prb );

  }

  if (p==0.0) cur_t = R_NegInf;

  return( cur_t );
}


// Lookup - 08
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qemg ( Rcpp::NumericVector p,
                                Rcpp::NumericVector mu,
                                Rcpp::NumericVector sigma,
                                Rcpp::NumericVector lambda,
                                double mnRT = -1.0,
                                double mxRT = 4.0,
                                int em_stop = 20,
                                double err = 0.001 ) {

  int N_p = p.size(); // Number of observations
  int N_mu = mu.size(); // Number of parameters
  int N_lambda = lambda.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int p_inc = 0;
  int mu_inc = 0;
  int lambda_inc = 0;
  int sigma_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_p, N_mu,
                                           N_lambda, N_sigma) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector p_v(N);
  Rcpp::NumericVector mu_v(N);
  Rcpp::NumericVector lambda_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    p_v(nv) = p(p_inc);
    mu_v(nv) = mu(mu_inc);
    lambda_v(nv) = lambda(lambda_inc);
    sigma_v(nv) = sigma(sigma_inc);

    p_inc = p_inc + 1;
    mu_inc = mu_inc + 1;
    lambda_inc = lambda_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_p==p_inc) p_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;
    if (N_lambda==lambda_inc) lambda_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Determine density
  for (int n = 0; n < N; n++) {
    out(n) = qemg_scl( p_v(n), mu_v(n), sigma_v(n), lambda_v(n),
        mnRT, mxRT, em_stop, err );
  }

  return ( out );
}

/*
library(seqmodels)

# Define generating parameters
prm = c( mu = -1, sigma = 1, lambda = .25 )

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
x = seq( -4, max(sim), length = 1000 )
lines( x, demg( x, prm[1], prm[2], prm[3] ), col ='orange', lwd = 2 )

# CDF
p = 1:length(sim)
p = p/length(sim)
x11();
plot( sort(sim), p, type = 'l', xlab = 'x', ylab = 'CDF', bty = 'l' )
lines( x, pemg(x,prm[1],prm[2],prm[3]), col='orange', lwd = 2 )
prb = seq(.1,.9,.1)
obs = quantile( sim, prb )
points( obs, prb, pch = 15, cex = 1.5 )
pred = qemg( prb, prm[1], prm[2], prm[3], mxRT = max(sim) )
points( pred, prb, pch = 19, cex = .75, col = 'red' )
*/
