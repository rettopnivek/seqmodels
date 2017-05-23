#include <Rcpp.h> // Includes certain libraries of functions
#include "levyfunctions.h"
#include "sigfunctions.h" // Scalar functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Vectorized functions for the random number generation, density,
distribution, quantile, and moments functions for the shifted
inverse gaussian distribution.

Index
Lookup - 01:  rinvgauss
Lookup - 02:  dinvgauss
Lookup - 03:  pinvgauss
Lookup - 04:  qinvgauss
Lookup - 05:  minvgauss
*/

// Lookup - 01
//' The Shifted Inverse Gaussian Distribution
//'
//' Random generation, density, distribution, and quantile functions
//' for the shifted inverse gaussian (or Wald) distribution,
//' parameterized for Brownian motion. \code{kappa} refers to the
//' threshold, \code{xi} refers to the rate of evidence accumulation
//' towards this threshold, \code{tau} is the shift in response times
//' and \code{sigma} refers to the within-trial variability for the
//' rate of evidence accumulation (the coefficient of drift, typically
//' fixed to 1).
//'
//' @param n the number of draws for random generation.
//' @param t a vector of times ( t > 0 ).
//' @param kappa a vector of thresholds determining when a decision
//'   terminates (kappa > 0).
//' @param xi a vector of drift rates, or rates of evidence accumulation
//'   (xi \eqn{\ge} 0).
//' @param tau a vector of shift parameters denoting the lowest
//'   possible time that can be observed (0 \eqn{\ge} tau < t).
//' @param sigma a vector of the within-trial variabilities
//'   (sigma > 0).
//' @param ln logical; if \code{TRUE}, probabilities are given as
//'   log(p).
//' @param lower_tail logical; if \code{TRUE} (default), probabilities
//'   are \eqn{P(X \le x)} otherwise \eqn{P( X > x)}.
//' @param bounds upper limit of the quantiles to explore  for the
//'   approximation via linear interpolation.
//' @param em_stop the maximum number of iterations to attempt to
//'   find the quantile via linear interpolation.
//' @param err the number of decimals places to approximate the
//'   cumulative probability during estimation of the quantile function.
//'
//' @section Details:
//' The inverse gaussian distribution describes the first passage times
//' through a positive threshold kappa for a space and time homogenous
//' Wiener diffusion process.
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
//' \code{dinvgauss} gives the density, \code{pinvgauss} gives the
//' distribution function, \code{qinvgauss} approximates the quantile
//' function, \code{minvgauss} computes the descriptive moments (mean,
//' variance, standard deviation, skew, and excess kurtosis), and
//' \code{rinvgauss} generates random deviates.
//'
//' The length of the result is determined by \code{n} for
//' \code{rinvgauss}, and is the maximum of the length of the
//' numerical arguments for the other functions.
//'
//' The numerical arguments other than \code{n} are recycled to the
//' length of the result.
//'
//' @section References:
//'
//' Dagpunar, J. (1988). Principles of Random Variate Generation.
//'   Oxford: Clarendon Press.
//'
//' Heathcote, A. (2004a). Fitting Wald and ex-Wald distributions to
//'   response time data: An example using functions for the S-PLUS
//'   package. Behavior Research Methods Instruments & Computers, 36,
//'   678 - 694.
//'
//' Heathcote, A. (2004b). rtfit.ssc. Retrieved May 5, 2017 from
//'   Psychonomic Society Web Archive:
//'   http://www.psychonomic.org/ARCHIVE/.
//'
//' @examples
//' # Density
//' dinvgauss( .9758, kappa = 1.0, xi = 1.0, tau = 0.3 )
//' # Distribution function
//' pinvgauss( .9758, kappa = 1.0, xi = 1.0, tau = 0.3 )
//' # Quantile function (Accurate to ~4 decimal places)
//' round( qinvgauss( p = .5, kappa = 1.0, xi = 1.0, tau = 0.3 ), 4 )
//' # Descriptive moments
//' minvgauss( kappa = 1.0, xi = 1.0, tau = 0.3 )
//'
//' # Simulation (No shift)
//' sim = rinvgauss( 1000, kappa = 0.8, xi = 2.0 )
//'
//' # Function to obtain maximum likelihood estimates
//' param_est = function( dat, tau_hat = 0 ) {
//'   # Estimate threshold and drift from first two moments of
//'   # data (Heathcote, 2004):
//'   dat = dat - tau_hat # Apply shift
//'   xi_hat = sqrt( mean( dat )/var( dat ) );
//'   kappa_hat = xi_hat * mean( dat )
//'   return( c( kappa = kappa_hat, xi = xi_hat, tau = tau_hat ) )
//' }
//' print( param_est( sim ) )
//'
//' # Non-zero shift parameter
//' sim = rinvgauss( 1000, kappa = 1.6, xi = 1.5, tau = .3 )
//'
//' # Estimating shift parameter
//'
//' # Function to compute sum of log-likelihoods
//' f = function( tau_hat, dat ) {
//'   prm = param_est( dat, tau_hat = tau_hat )
//'   sll = sum( dinvgauss( dat, prm[1], prm[2], tau = prm[3], ln = T ) )
//'   return( sll )
//' }
//' tau_hat = optimize( f, c( 0.0, min( sim ) ), dat = sim, maximum = T )
//' print( param_est( sim, tau_hat = tau_hat$maximum ) )
//'
//' # Plotting
//' layout( matrix( 1:4, 2, 2, byrow = T ) )
//' # Parameters
//' prm = c( k = 0.8, x = 1.6, t = 0.3, s = 1.0 )
//' # Density
//' obj = quickdist( 'sig', 'PDF', prm )
//' plot( obj ); lines( obj )
//' # CDF
//' obj = quickdist( 'sig', 'CDF', prm )
//' plot( obj ); lines( obj )
//' # Quantiles
//' obj = quickdist( 'sig', 'QF', prm, x = seq( .2, .8, .2 ) )
//' plot( obj ); prb = seq( .2, .8, .2 )
//' abline( h = prb, lty = 2 ); lines( obj, type = 'b', pch = 19 )
//' # Hazard function
//' obj = quickdist( 'sig', 'HF', prm )
//' plot( obj ); lines( obj )
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rinvgauss( int n,
                               Rcpp::NumericVector kappa,
                               Rcpp::NumericVector xi,
                               Rcpp::NumericVector tau =
                                 Rcpp::NumericVector::create(0.0),
                               Rcpp::NumericVector sigma =
                                 Rcpp::NumericVector::create(1.0) ) {

  // Set output vector
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 6 );

  // Fill matrix
  expanded_input = kappa[ create_index( n, kappa.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = Rcpp::rchisq( n, 1 );
  input.column(4) = expanded_input;
  expanded_input = Rcpp::runif( n, 0.0, 1.0 );
  input.column(5) = expanded_input;

  // Generate draws
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(6);
    for ( int c = 0; c < 6; c++ ) prm[c] = input(r,c);
    out(r) = rinvgauss_scl( prm );
  }

  return( out );
}

// Lookup - 02
//' @rdname rinvgauss
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dinvgauss( Rcpp::NumericVector t,
                               Rcpp::NumericVector kappa,
                               Rcpp::NumericVector xi,
                               Rcpp::NumericVector tau =
                                 Rcpp::NumericVector::create(0.0),
                               Rcpp::NumericVector sigma =
                                 Rcpp::NumericVector::create(1.0),
                               bool ln = false ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    t.size(), kappa.size(), xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 5 );

  // Fill matrix
  expanded_input = t[ create_index( n, t.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = kappa[ create_index( n, kappa.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(4) = expanded_input;

  // Determine if log-likelihoods should be returned
  int ln_val = 0; if ( ln ) ln_val = 1;

  // Calculate PDF
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(6);
    for ( int c = 0; c < 5; c++ ) prm[c] = input(r,c);
    prm[5] = ln_val;
    out(r) = dinvgauss_scl( prm );
  }

  return( out );
}

// Lookup - 03
//' @rdname rinvgauss
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pinvgauss( Rcpp::NumericVector t,
                               Rcpp::NumericVector kappa,
                               Rcpp::NumericVector xi,
                               Rcpp::NumericVector tau =
                                 Rcpp::NumericVector::create(0.0),
                               Rcpp::NumericVector sigma =
                                 Rcpp::NumericVector::create(1.0),
                               bool ln = false,
                               bool lower_tail = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    t.size(), kappa.size(), xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 5 );

  // Fill matrix
  expanded_input = t[ create_index( n, t.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = kappa[ create_index( n, kappa.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(4) = expanded_input;

  // Calculate CDF
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(5);
    for ( int c = 0; c < 5; c++ ) prm[c] = input(r,c);
    out(r) = pinvgauss_scl( prm );
  }

  if ( !lower_tail ) out = 1.0 - out;
  if ( ln ) out = log( out );

  return( out );
}

// Lookup - 04
//' @rdname rinvgauss
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qinvgauss( Rcpp::NumericVector p,
                               Rcpp::NumericVector kappa,
                               Rcpp::NumericVector xi,
                               Rcpp::NumericVector tau =
                                 Rcpp::NumericVector::create(0.0),
                               Rcpp::NumericVector sigma =
                                 Rcpp::NumericVector::create(1.0),
                               double bounds = 3.0,
                               double em_stop = 20,
                               double err = 1e-8 ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    p.size(), kappa.size(), xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 5 );

  // Fill matrix
  expanded_input = p[ create_index( n, p.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = kappa[ create_index( n, kappa.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(4) = expanded_input;

  // Calculate quantiles
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(8);
    for ( int c = 0; c < 5; c++ ) prm[c] = input(r,c);
    prm[5] = bounds; prm[6] = em_stop; prm[7] = err;
    out(r) = qinvgauss_scl( prm );
  }

  return( out );
}

// Lookup - 05
//' @rdname rinvgauss
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame minvgauss( Rcpp::NumericVector kappa,
                           Rcpp::NumericVector xi,
                           Rcpp::NumericVector tau =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector sigma =
                             Rcpp::NumericVector::create(1.0) ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    kappa.size(), xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericMatrix out(n,5);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 4 );

  // Fill matrix
  expanded_input = kappa[ create_index( n, kappa.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(3) = expanded_input;

  // Calculate moments
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(4);
    for ( int c = 0; c < 4; c++ ) prm[c] = input(r,c);
    std::vector<double> res;
    res = minvgauss_scl( prm );
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
