#include <RcppParallel.h>
#include "levyfunctions.h" // Scalar functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Vectorized functions for the density, distribution, quantile,
and random generation functions of the Levy distribution.

Index
Lookup - 01:  dlevy
Lookup - 02;  plevy
Lookup - 03:  qlevy
Lookup - 04:  rlevy
*/

// Lookup - 01
//' The Levy Distribution
//'
//' Density, distribution, random generation, and quantile functions
//' for the Levy distribution, where \code{mu} is a location parameter
//' and \code{sigma} is a scale parameter.
//'
//' @param n the number of draws for random generation.
//' @param x,q a vector of quantiles (must be greater than mu).
//' @param mu a vector of location parameters.
//' @param sigma a vector of scale parameters (sigma > 0).
//' @param ln logical; if \code{TRUE}, probabilities are given as
//'   log(p).
//' @param lower_tail logical; if \code{TRUE} (default), probabilities
//'   are \eqn{P(X \le x)} otherwise \eqn{P( X > x)}.
//'
//' @section Details:
//' A Levy distribution, among other things, can describe the finishing
//' times for a one boundary wiener process when the drift rate is
//' fixed to zero.
//'
//' The mean and variance for the Levy distribution are non-finite.
//'
//' @return
//' \code{dlevy} gives the density, \code{plevy} gives the
//' distribution function, \code{qlevy} gives the quantile function
//' and \code{rlevy} generates random deviates.
//'
//' The length of the result is determined by \code{n} for
//' \code{rlevy}, and is the maximum of the lengths of the numerical
//' arguments for the other functions.
//'
//' The numerical arguments other than \code{n} are recycled to the
//' length of the result.
//'
//' @section References:
//'
//' Applebaum, D. (2010). Lectures on Levy processes and stochastic
//'   calculus, Braunschweig; Lecture 2: Levy processes. Retrieved from
//'   http://www.applebaum.staff.shef.ac.uk/Brauns2notes.pdf.
//'
//' Siegrist, K. (1997). The Levy distribution. Retrieved from
//'   http://www.math.uah.edu/stat/special/Levy.html
//'
//' @examples
//' # Density
//' dlevy( x = 2.199, mu = 0.0, sigma = 1.0 )
//' # Distribution function
//' plevy( q = 2.199, mu = 0.0, sigma = 1.0 )
//' # Quantile function
//' qlevy( p = .5, mu = 0.0, sigma = 1.0 )
//'
//' # Simulations
//' sim = rlevy( n = 1000, mu = 0.0, sigma = 1.0 )
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dlevy ( Rcpp::NumericVector x,
                            Rcpp::NumericVector mu,
                            Rcpp::NumericVector sigma,
                            bool ln = false ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    x.size(), mu.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 3 );

  // Fill matrix
  expanded_input = x[ create_index( n, x.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;

  // Determine if log-likelihoods should be returned
  int ln_val = 0; if ( ln ) ln_val = 1;

  // Calculate the density
  for (int r = 0; r < n; r++) {

    std::vector<double> prm(4);
    for ( int c = 0; c < 3; c++ ) prm[c] = input(r,c);
    prm[3] = ln_val;
    out(r) = dlevy_scl( prm );
  }

  return( out );
}

// Lookup - 02
//' @rdname dlevy
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector plevy ( Rcpp::NumericVector q,
                            Rcpp::NumericVector mu,
                            Rcpp::NumericVector sigma,
                            bool lower_tail = true,
                            bool ln = false ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    q.size(), mu.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 3 );

  // Fill matrix
  expanded_input = q[ create_index( n, q.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;

  // Calculate the density
  for (int r = 0; r < n; r++) {

    std::vector<double> prm(3);
    for ( int c = 0; c < 3; c++ ) prm[c] = input(r,c);
    out(r) = plevy_scl( prm );
  }

  if ( !lower_tail ) out = 1.0 - out;
  if ( ln ) out = log( out );

  return( out );
}

// Lookup - 03
//' @rdname dlevy
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qlevy ( Rcpp::NumericVector p,
                            Rcpp::NumericVector mu,
                            Rcpp::NumericVector sigma ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    p.size(), mu.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 3 );

  // Fill matrix
  expanded_input = p[ create_index( n, p.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;

  // Calculate the density
  for (int r = 0; r < n; r++) {

    std::vector<double> prm(3);
    for ( int c = 0; c < 3; c++ ) prm[c] = input(r,c);
    out(r) = qlevy_scl( prm );
  }

  return( out );
}

// Lookup - 04
//' @rdname dlevy
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rlevy( int n, Rcpp::NumericVector mu,
                           Rcpp::NumericVector sigma ) {

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 3 );

  // Fill matrix
  expanded_input = Rcpp::runif(n,0.0,1.0);
  input.column(0) = expanded_input;
  expanded_input = mu[ create_index( n, mu.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(2) = expanded_input;

  // Calculate the density
  for (int r = 0; r < n; r++) {

    std::vector<double> prm(3);
    for ( int c = 0; c < 3; c++ ) prm[c] = input(r,c);
    out(r) = rlevy_scl( prm );
  }

  return( out );
}

