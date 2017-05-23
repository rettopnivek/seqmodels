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
//' towards this threshold, \code{tau} is the mean for the
//' exponentially distribution shift values, and \code{sigma} is
//' the  within-trial variability for the rate of evidence
//' accumulation (the coefficient of drift, typically fit to 1 for
//' identification purposes).
//'
//' @param n the number of draws for random generation.
//' @param t a vector of times ( t > 0 ).
//' @param kappa a vector of thresholds determining when a decision
//'   terminates (kappa > 0).
//' @param xi a vector of drift rates, or rates of evidence accumulation
//'   (xi \eqn{\ge} 0).
//' @param tau a vector of inverse rate parameters (the mean of the
//'   exponentially distributed residual latency; tau > 0).
//' @param sigma a vector of the within-trial variabilities (the
//'   coefficient of drift; sigma > 0).
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
//' The ex-Wald distribution is the result if one sums two random
//' variables, one distributed according to the inverse Gaussian
//' distribution and the other distributed according to an
//' exponential distribution with rate equal to 1/tau.
//'
//' The density is therefore computed from the convolution of the
//' Wald and exponential distributions. Schwarz (2002) provides a
//' solution for this convolution. If \eqn{\xi^2} is greater than
//' \eqn{2 \sigma^2 / \tau }, a straightforward analytic solution
//' exists (Equation 18). Otherwise, the density is the real
//' part of a complex function. Schwarz provides a solution in
//' this case using the complex error function (Equation 9 and 22).
//' To compute the real and imaginary parts of the complex error
//' function, the density function uses an algorithm from the
//' Faddeeva package (Johnson, 2012), with an additional wrapper
//' taken from the S-PLUS script of Heathcote (2004b).
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
//'
//' Heathcote, A. (2004a). Fitting Wald and ex-Wald distributions to
//'   response time data: An example using functions for the S-PLUS
//'   package. Behavior Research Methods Instruments & Computers, 36,
//'   678 - 694.
//'
//' Heathcote, A. (2004b). rtfit.ssc. Retrieved May 5, 2017 from
//'   Psychonomic Society Web Archive: http://www.psychonomic.org/ARCHIVE/.
//'
//' Johnson, S. G. (2012). Faddeeva Package [Computer software].
//'   Retrieved from
//'   http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package#License
//'
//' Schwarz, W. (2002). On the convolution of the inverse Gaussian and
//'   the exponential random variables. Communications in Statistics,
//'   Theory and Methods, 31, 2113 - 2121.
//'
//' @examples
//' Forthcoming
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rexwald( int n,
                             Rcpp::NumericVector kappa,
                             Rcpp::NumericVector xi,
                             Rcpp::NumericVector tau,
                             Rcpp::NumericVector sigma =
                               Rcpp::NumericVector::create(1.0) ) {

  // Set output vector
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 7 );

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
  expanded_input = Rcpp::runif( n, 0.0, 1.0 );
  input.column(6) = expanded_input;

  // Generate draws
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(7);
    for ( int c = 0; c < 7; c++ ) prm[c] = input(r,c);
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
                             Rcpp::NumericVector sigma =
                               Rcpp::NumericVector::create(1.0),
                             bool ln = false,
                             bool ni = false ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    t.size(), kappa.size(), xi.size(), tau.size(),
    sigma.size() ) );

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
  // Determine if numerical integration should be used
  int ni_yes = 0; if ( ni ) ni_yes = 1;

  // Calculate PDF
  for (int r = 0; r < n; r++) {
    std::vector<double> prm(7);
    for ( int c = 0; c < 5; c++ ) prm[c] = input(r,c);
    prm[5] = ln_val;
    prm[6] = ni_yes;
    out(r) = dexwald_scl( prm );
  }

  return( out );
}
