#include <Rcpp.h> // Includes certain libraries of functions
#include "emgfunctions.h" // Scalar functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Vectorized functions for the random number generation, density,
distribution, quantile, and moments functions for the exponentially
modified guassian distribution.

Index
Lookup - 01:  remg
Lookup - 02:  demg
Lookup - 03:  pemg
Lookup - 04:  qemg
Lookup - 05:  memg
*/

// Lookup - 01
//' The Exponentially Modified Gaussian Distribution
//'
//' Random generation, density, distribution, quantile, and
//' descriptive moments functions for the the convolution of gaussian
//' and exponential random variables (the ex-gaussian distribution),
//' with location equal to \code{mu}, scale equal to \code{sigma} and
//' rate equal to \code{lambda}.
//'
//' @param x,q vector of quantiles.
//' @param p vector of probabilities.
//' @param n number of observations to be generated.
//' @param mu vector of location parameters for the gaussian
//'   variable.
//' @param sigma vector of scale parameters (i.e., standard
//'   deviations) for the gaussian variable (sigma > 0).
//' @param lambda vector of rate parameters for the exponential
//'   variable (lambda > 0).
//' @param ln logical; if \code{TRUE}, probabilities are given as
//'   log(p).
//' @param lower_tail logical; if \code{TRUE} (default), probabilities
//'   are \eqn{P(X \le x)} otherwise \eqn{P( X > x)}.
//' @param bounds lower and upper limits of the quantiles to explore
//'   for the approximation via linear interpolation.
//' @param em_stop the maximum number of iterations to attempt to
//'   find the quantile via linear interpolation.
//' @param err the number of decimals places to approximate the
//'   cumulative probability during estimation of the quantile function.
//'
//' @section Details:
//' An exponentially modified gaussian distribution describes the sum of
//' independent normal and exponential random variables, possessing
//' a characteristic positive skew due to the exponential variable.
//' The ex-gaussian distribution is therefore useful for providing
//' descriptive fits to response time distributions (e.g., Luce, 1986).
//'
//' A linear interpolation approach is used to approximate the
//' quantile function, estimating the inverse of the cumulative
//' distribution function via an iterative procedure. When
//' the precision of this estimate is set to 8 decimal places,
//' the approximation will be typically accurate to about half of a
//' millisecond.
//'
//' The example section demonstrates how to compute maximum likelihood
//' estimates based on the moments from a set of data.
//'
//' @return
//' \code{demg} gives the density, \code{pemg} gives the
//' distribution function, \code{qemg} gives the quantile function,
//' \code{memg} computes the descriptive moments (mean, variance,
//' standard deviation, skew, and excess kurtosis), and \code{remg}
//' generates random deviates.
//'
//' The length of the result is determined by \code{n} for
//' \code{remg}, and is the maximum of the lengths of the numerical
//' arguments for the other functions.
//'
//' The numerical arguments other than \code{n} are recycled to the
//' length of the result.
//'
//' @section References:
//'
//' Luce, R. D. (1986). Response times: Their role in inferring
//'   elementary mental organization. New York, New York: Oxford
//'   University Press.
//'
//' Terriberry, T. B. (2007). Computing Higher-Order Moments Online.
//'   Retrieved from https://people.xiph.org/~tterribe/notes/homs.html
//'
//' @examples
//' # Density function
//' demg( x = 0.8758, mu = 0.0, sigma = 1.0, lambda = 1.0 )
//' # Distribution function
//' pemg( q = 0.8758, mu = 0.0, sigma = 1.0, lambda = 1.0 )
//' # Quantile function (Accurate to ~4 decimal places)
//' round( qemg( p = .5, mu = 0.0, sigma = 1.0, lambda = 1.0 ), 4 )
//' # Descriptive moments
//' memg( mu = 0.44, sigma = 0.07, lambda = 2.0 )
//'
//' # Simulation
//' sim = remg( n = 1000, mu = 0.44, sigma = 0.07, lambda = 2.0 );
//'
//' # Function to obtain maximum likelihood estimates
//' param_est = function( dat ) {
//'   # Compute 1st, 2nd, and 3rd moments ( Terriberry, 2007).
//'   n = 0; m = 0; m2 = 0; m3 = 0;
//'   for ( k in 1:length( dat ) ) {
//'     n1 = n; n = n + 1
//'     term1 = dat[k] - m; term2 = term1/ n; term3 = term1 * term2 * n1
//'     # Update first moment
//'     m = m + term2
//'     # Update third moment
//'     m3 = m3 + term3 + term2 * (n - 2) - 3 * term2 * m2
//'     # Update second moment
//'     m2 = m2 + term3
//'   }
//'   # Compute standard deviation of sample
//'   s = sqrt( m2 / ( n - 1.0 ) )
//'   # Compute skewness of sample
//'   y = sqrt( n ) * m3 / ( m2^1.5 );
//'   # Estimate parameters
//'   mu_hat = m - s * ( y/2 )^(1/3)
//'   sigma_hat = sqrt( (s^2) * ( 1 - (y/2)^(2/3) ) )
//'   lambda_hat = 1/( s * (y/2)^(1/3) )
//'   return( c( mu = mu_hat, sigma = sigma_hat, lambda = lambda_hat ) )
//' }
//'
//' print( param_est( sim ) )
//'
//' # Plotting
//' layout( matrix( 1:4, 2, 2, byrow = T ) )
//' # Parameters
//' prm = c( m = .44, s = .07, l = 2.0 )
//' # Density
//' obj = quickdist( 'emg', 'PDF', prm )
//' plot( obj ); lines( obj )
//' # CDF
//' obj = quickdist( 'emg', 'CDF', prm )
//' plot( obj ); lines( obj )
//' # Quantiles
//' obj = quickdist( 'emg', 'QF', prm, x = seq( .2, .8, .2 ) )
//' plot( obj ); prb = seq( .2, .8, .2 );
//' abline( h = prb, lty = 2 );
//' lines( obj, type = 'b', pch = 19 )
//' # Hazard function
//' obj = quickdist( 'emg', 'HF', prm )
//' plot( obj ); lines( obj )
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector remg( int n, Rcpp::NumericVector mu,
                          Rcpp::NumericVector sigma,
                          Rcpp::NumericVector lambda ) {
  // Set output vector
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 3 );

  // Fill matrix
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = lambda[ create_index( n, lambda.size() ) ];
  input.column(2) = expanded_input;

  // Generate draws
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(3);
    for ( int c = 0; c < 3; c++ ) prm[c] = input(r,c);
    out(r) = remg_scl( prm );
  }

  return( out );
}

// Lookup - 02
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector demg(Rcpp::NumericVector x,
                         Rcpp::NumericVector mu,
                         Rcpp::NumericVector sigma,
                         Rcpp::NumericVector lambda,
                         bool ln = false ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    x.size(), mu.size(), lambda.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 4 );

  // Fill matrix
  expanded_input = x[ create_index( n, x.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = lambda[ create_index( n, lambda.size() ) ];
  input.column(3) = expanded_input;

  // Determine if log-likelihoods should be returned
  int ln_val = 0; if ( ln ) ln_val = 1;

  // Calculate PDF
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(5);
    for ( int c = 0; c < 4; c++ ) prm[c] = input(r,c);
    prm[4] = ln_val;
    out(r) = demg_scl( prm );
  }

  return( out );
}

// Lookup - 03
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pemg(Rcpp::NumericVector q,
                         Rcpp::NumericVector mu,
                         Rcpp::NumericVector sigma,
                         Rcpp::NumericVector lambda,
                         bool ln = false,
                         bool lower_tail = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    q.size(), mu.size(), lambda.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 4 );

  // Fill matrix
  expanded_input = q[ create_index( n, q.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = lambda[ create_index( n, lambda.size() ) ];
  input.column(3) = expanded_input;

  // Calculate CDF
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(4);
    for ( int c = 0; c < 4; c++ ) prm[c] = input(r,c);
    out(r) = pemg_scl( prm );
  }

  if ( !lower_tail ) out = 1.0 - out;
  if ( ln ) out = log( out );

  return( out );
}

// Lookup - 04
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qemg(Rcpp::NumericVector p,
                         Rcpp::NumericVector mu,
                         Rcpp::NumericVector sigma,
                         Rcpp::NumericVector lambda,
                         Rcpp::NumericVector bounds =
                           Rcpp::NumericVector::create( 0.0, 3.0),
                         double em_stop = 20, double err = 1e-8 ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    p.size(), mu.size(), lambda.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 4 );

  // Fill matrix
  expanded_input = p[ create_index( n, p.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = lambda[ create_index( n, lambda.size() ) ];
  input.column(3) = expanded_input;

  // Calculate quantiles
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(8);
    for ( int c = 0; c < 4; c++ ) prm[c] = input(r,c);
    prm[4] = bounds[0]; prm[5] = bounds[1];
    prm[6] = em_stop; prm[7] = err;
    out(r) = qemg_scl( prm );
  }

  return( out );
}

// Lookup - 05
//' @rdname remg
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame memg( Rcpp::NumericVector mu,
                      Rcpp::NumericVector sigma,
                      Rcpp::NumericVector lambda ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    mu.size(), lambda.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericMatrix out(n,5);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 3 );

  // Fill matrix
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = lambda[ create_index( n, lambda.size() ) ];
  input.column(2) = expanded_input;

  // Calculate moments
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(3);
    for ( int c = 0; c < 3; c++ ) prm[c] = input(r,c);
    std::vector<double> res;
    res = memg_scl( prm );
    for ( int c = 0; c < 5; c++ ) out(r,c) = res[c];
  }

  // Convert to a data-frame
  return Rcpp::DataFrame::create(
    Rcpp::_["mean"]= out.column(0),
    Rcpp::_["var"]= out.column(1),
    Rcpp::_["sd"]= out.column(2),
    Rcpp::_["skew"]= out.column(3),
    Rcpp::_["ex_kurtosis"]= out.column(4) );
}
