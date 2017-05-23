#include <RcppParallel.h>
#include "levyfunctions.h" // For drift of 0
#include "sigfunctions.h" // For inverse gaussian
#include "ewfunctions.h" // Scalar functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Vectorized functions for the random number generation, density,
distribution, and quantile functions for the convolution
of the exponential and Wald distributions.

Index
Lookup - 01:  rexwald
Lookup - 02:  dexwald
*/

// Lookup - 01
//' The ex-Wald distribution
//'
//' Random generation, density, distribution, and quantile functions
//' for the convolution of the exponential and Wald distributions,
//' parameterized for Brownian motion. \code{kappa} refers to the
//' threshold, \code{xi} refers to the rate of evidence accumulation
//' towards this threshold, and \code{tau} is the rate parameter
//' governing the exponentially distribution shift values. Within-trial
//' variability for the rate of evidence accumulation (the coefficient
//' of drift) is fixed to 1.
//'
//' @param n the number of draws for random generation.
//' @param t a vector of times ( t > 0 ).
//' @param kappa a vector of thresholds determining when a decision
//'   terminates (kappa > 0).
//' @param xi a vector of drift rates, or rates of evidence accumulation
//'   (xi \eqn{\ge} 0).
//' @param tau a vector of inverse rate parameters (the mean of the
//'   exponentially distributed residual latency; tau > 0).
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
//' Forthcoming
//'
//' A linear interpolation approach is used to approximate the
//' quantile function, estimating the inverse of the cumulative
//' distribution function via an iterative procedure. When
//' the precision of this estimate is set to 8 decimal places,
//' the approximation will be typically accurate to about half of a
//' millisecond.
//'
//' @return
//' \code{dexwald} gives the density, \code{pexwald} gives the
//' distribution function, \code{qexwald} approximates the quantile
//' function, and \code{rexwald} generates random deviates.
//'
//' The length of the result is determined by \code{n} for
//' \code{rexwald}, and is the maximum of the length of the
//' numerical arguments for the other functions.
//'
//' The numerical arguments other than \code{n} are recycled to the
//' length of the result.
//'
//' @section References:
//' Forthcoming
//'
//' @examples
//' Forthcoming
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rexwald( int n,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector xi,
                             Rcpp::NumericVector tau ) {

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
  expanded_input = Rcpp::rchisq( n, 1 );
  input.column(3) = expanded_input;
  expanded_input = Rcpp::runif( n, 0.0, 1.0 );
  input.column(4) = expanded_input;
  expanded_input = Rcpp::runif( n, 0.0, 1.0 );
  input.column(5) = expanded_input;

  // Generate draws
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(6);
    for ( int c = 0; c < 6; c++ ) prm[c] = input(r,c);
    out(r) = rexwald_scl( prm );
  }

  return( out );
}

// Lookup - 02
//' @rdname rexwald
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dexwald( Rcpp::NumericVector t,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector xi,
                             Rcpp::NumericVector tau,
                             bool ln = false,
                             bool ni = false ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    t.size(), kappa.size(), xi.size(), tau.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 4 );

  // Fill matrix
  expanded_input = t[ create_index( n, t.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = kappa[ create_index( n, kappa.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(3) = expanded_input;

  // Determine if log-likelihoods should be returned
  int ln_val = 0; if ( ln ) ln_val = 1;
  // Determine if numerical integration should be used
  int ni_yes = 0; if ( ni ) ni_yes = 1;

  // Calculate PDF
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(6);
    for ( int c = 0; c < 4; c++ ) prm[c] = input(r,c);
    prm[4] = ln_val;
    prm[5] = ni_yes;
    out(r) = dexwald_scl( prm );
  }

  return( out );
}
