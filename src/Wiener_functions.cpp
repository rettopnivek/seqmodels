#include <RcppParallel.h>
#include "wpparamverify.h"
#include "navarrofuss.h"
#include "blurtonetal.h"
#include "wienerfunctions.h" // Scalar functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Vectorized functions for the density, distribution, quantile, and
random generation functions for the two boundary wiener process.

Package development:
library(devtools)
library(roxygen2)
1st step: Create folder to store package
2nd step: In Rstudio, set working directory to package folder
          and type 'devtools::use_rcpp()'. This creates the
          basic framework for the package using Rcpp.
3rd step: Select 'Open project' in the file menu and open
          the newly created project file in the package folder.
4th step: Create a C++ file and include the lines:
          //' @useDynLib seqmodels
          //' @importFrom Rcpp sourceCpp

Example of roxygen2 documentation for Rcpp
 //' Title
 //'
 //' Short description.
 //'
 //' @param variable description of variable.
 //'
 //' @return Function's return value(s).
 //' @export
 // [[Rcpp::export]]

Template for documentation of non-exported functions
  Purpose:
  ...
  Arguments:
  ...
  Returns:
  ...

// Printing to R console
Rcpp::Rcout << "Debugging example" << std::endl;

Index
Lookup - 01:  dwienerWorker
Lookup - 02:  dwiener
Lookup - 03:  pwienerWorker
Lookup - 04:  pwiener
Lookup - 05:  qwienerWorker
Lookup - 06:  qwiener
Lookup - 07:  rwienerWorker
Lookup - 08:  rwiener

### TO DO ###
Add moments function
Update .rd page (add plotting examples)
*/

//' @useDynLib seqmodels
//' @importFrom Rcpp sourceCpp

// Lookup - 01
// RcppParallel worker function

struct dwienerWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RVector<double> output;

  // initialize with source and destination
  dwienerWorker( const Rcpp::NumericMatrix input,
                 Rcpp::NumericVector output)
                 : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(9);

      for (int i = 0; i < 9; i++) { prm[i] = input(j,i); }

      output[j] = dwiener_scl( prm );
    }
  }
};

// Lookup - 02
//' Two-boundary Wiener Process for Choice and Response Times
//'
//' Density, distribution, quantile, and random generation functions
//' for a two-boundary wiener process that can be applied to
//' choice and response times (e.g., Luce, 1986; Ratcliff, 1978).
//' \code{alpha} refers to the boundary separation, \code{theta}
//' refers to the proportion governing the start point of accumulation,
//' \code{xi} refers to the rate of evidence accumulation (drift rate),
//' \code{tau} refers to the residual latency (e.g., motor and
//' encoding processes), and \code{sigma} refers to the within-trial
//' variability of evidence accumulation (the coefficient of drift;
//' typically set to 1 or 0.1).
//'
//' @param n the number of draws for random generation.
//' @param rt a vector of responses times ( \code{rt} > 0 ).
//' @param p a vector of probabilities.
//' @param ch a vector of accuracy/choice values ( \code{ch} = {0,1} ).
//' @param alpha a vector of upper boundaries at which the evidence
//'   accumulation terminations.
//' @param theta a vector of proportions determining the starting
//'   point for the evidence accumulation, where the starting point
//'   \eqn{\zeta} = \code{alpha}*\code{theta} ( 0 \eqn{\ge} \code{theta}
//'   \eqn{\ge} 1 ).
//' @param xi a vector of drift rates, the rate of evidence accumulation
//'   ( \code{xi} > 0 ).
//' @param tau a vector of residual latencies for the non-decision
//'   component ( \code{tau} > 0 ).
//' @param sigma a vector giving the coefficients of drift (also known as
//'   within-trial variability; \code{sigma} > 0 ).
//' @param eps the margin of error for the infinite sums being calculated.
//' @param ln logical; if \code{TRUE}, probabilities are given as
//'   log(p).
//' @param joint logical; if \code{FALSE} the conditional density
//'   (normalized to integrate to one) is returned. Otherwise, the
//'   joint density (integrating to the choice probability) is
//'   returned.
//' @param lower_tail logical; if \code{TRUE} (default), probabilities
//'   are \eqn{P(X \le x)} otherwise \eqn{P( X > x)}.
//' @param bounds upper limit of the quantiles to explore
//'   for the approximation via linear interpolation.
//' @param em_stop the maximum number of iterations to attempt to
//'   find the quantile via linear interpolation.
//' @param err the number of decimals places to approximate the
//'   cumulative probability during estimation of the quantile function.
//' @param parYes logical; if \code{TRUE} the code is run in parallel.
//'
//' @section Details:
//' The density function is based on the implementation of Navarro
//' and Fuss (2009). The distribution function is based on the
//' implementation of Blurton et al. (2012).
//'
//' A linear interpolation approach is used to approximate the
//' quantile function and to random deviates by estimating the
//' inverse of the cumulative distribution function via an
//' iterative procedure. When the precision of this estimate is
//' set to 8 decimal places, the approximation will be typically
//' accurate to about half of a millisecond.
//'
//' @return
//' \code{dwiener} gives the density, \code{pwiener} gives the
//' distribution function, \code{qwiener} approximates the quantile
//' function, and \code{rwiener} generates random deviates.
//'
//' The length of the result is determined by \code{n} for \code{rwiener},
//' and is the maximum of the length of the numerical arguments for
//' the other functions.
//'
//' The numerical arguments other than \code{n} are recycled to the
//' length of the result.
//'
//' @section References:
//'
//' Blurton, S. P., Kesselmeier, M., & Gondan, M. (2012). Fast and
//'   accurate calculations for cumulative first-passage time distributions
//'   in Wiener diffusion models. Journal of Mathematical Psychology,
//'   56, 470-475.
//'
//' Luce, R. D. (1986). Response times: Their role in inferring
//'   elementary mental organization. New York, New York: Oxford University
//'   Press.
//'
//' Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate calculations
//'   for first-passage times in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 53, 222-230.
//'
//' Ratcliff, R. (1978). A theory of memory retrieval. Psychological
//'   review, 85, 59 - 108.
//'
//' @examples
//' # Density
//' dwiener( rt = 0.6, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
//'   xi = 1.0, tau = 0.3 )
//' # Distribution function
//' pwiener( rt = 0.6, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
//'   xi = 1.0, tau = 0.3 )
//' # Choice probabilities
//' pwiener( rt = Inf, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
//'   xi = 1.0, tau = 0.3 )
//' # Quantile function (Accurate to ~4 decimal places)
//' round( qwiener( p = .3499, ch = c( 1, 0 ), alpha = 1.6, theta = 0.5,
//'   xi = 1.0, tau = 0.3 ), 4 )
//'
//' # Simulate values
//' sim = rwiener( n = 100, alpha = 0.8, theta = 0.6,
//'   xi = 0.0, tau = 0.3 )
//'
//' # Plotting
//' layout( matrix( 1:4, 2, 2, byrow = T ) )
//' # Parameters
//' prm = c( a = 1.2, z = .4, v = 1.0, t0 = 0.3 )
//' # Density
//' obj = quickdist( 'wp', 'PDF', prm )
//' plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
//' # CDF
//' obj = quickdist( 'wp', 'CDF', prm )
//' plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
//' # Quantiles
//' obj = quickdist( 'wp', 'QF', prm, x = seq( .2, .8, .2 ) )
//' plot( obj ); prb = seq( .2, .8, .2 )
//' abline( h = prb, lty = 2 )
//' # Conditional, not joint
//' lines( obj, type = 'b', pch = 19, weight = 1 )
//' lines( obj, ch = 0, type = 'b', pch = 17, lty = 2, weight = 1 )
//' # Hazard function
//' obj = quickdist( 'wp', 'HF', prm )
//' plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
//'
//' @export
// [[Rcpp::export]]

Rcpp::NumericVector dwiener( Rcpp::NumericVector rt,
                             Rcpp::NumericVector ch,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector theta,
                             Rcpp::NumericVector xi,
                             Rcpp::NumericVector tau,
                             Rcpp::NumericVector sigma =
                               Rcpp::NumericVector::create(1.0),
                             bool ln = false,
                             bool joint = true,
                             double eps = 1e-29,
                             bool parYes = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    rt.size(), ch.size(), alpha.size(), theta.size(),
    xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 9 );

  // Fill matrix
  expanded_input = rt[ create_index( n, rt.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = ch[ create_index( n, ch.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = alpha[ create_index( n, alpha.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = theta[ create_index( n, theta.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(6) = expanded_input;
  double joint_val = 1.0; // Extract indicator for joint density
  if ( !joint ) joint_val = 0.0;
  for (int k = 0; k < n; k++) {
    input(k,7) = eps;
    input(k,8) = joint_val;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate likelihood
  if ( !parYes ) {

    for (int j = 0; j < n; j++) {

      std::vector<double> prm(9);

      for (int i = 0; i < 9; i++)
        prm[i] = input(j,i);

      out(j) = dwiener_scl( prm );
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    dwienerWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  // If indicated, compute log of the density
  if ( ln ) out = log(out);

  return( out );
}

// Lookup - 03
// RcppParallel worker function

struct pwienerWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RVector<double> output;

  // initialize with source and destination
  pwienerWorker( const Rcpp::NumericMatrix input,
                 Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(9);

      for (int i = 0; i < 9; i++) { prm[i] = input(j,i); }

      output[j] = pwiener_scl( prm );
    }
  }
};

// Lookup - 04
//' @rdname dwiener
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pwiener( Rcpp::NumericVector rt,
                             Rcpp::NumericVector ch,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector theta,
                             Rcpp::NumericVector xi,
                             Rcpp::NumericVector tau,
                             Rcpp::NumericVector sigma =
                               Rcpp::NumericVector::create(1.0),
                               bool ln = false,
                               bool joint = true,
                               bool lower_tail = true,
                               double eps = 1e-29,
                               bool parYes = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    rt.size(), ch.size(), alpha.size(), theta.size(),
    xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 9 );

  // Fill matrix
  expanded_input = rt[ create_index( n, rt.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = ch[ create_index( n, ch.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = alpha[ create_index( n, alpha.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = theta[ create_index( n, theta.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(6) = expanded_input;
  double joint_val = 1.0; // Extract indicator for joint distribution
  if ( !joint ) joint_val = 0.0;
  for (int k = 0; k < n; k++) {
    input(k,7) = eps;
    input(k,8) = joint_val;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate distribution function
  if ( !parYes ) {

    for (int j = 0; j < n; j++) {

      std::vector<double> prm(9);

      for (int i = 0; i < 9; i++) { prm[i] = input(j,i); }

      out(j) = pwiener_scl( prm );
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    pwienerWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  if ( !lower_tail ) {
    for ( int l = 0; l < n; l++ ) {
      if ( out( l ) > 0.0 ) out( l ) = 1.0 - out( l );
    }
  }

  if ( ln ) out = log( out );

  return( out );
}

// Lookup - 05
// RcppParallel worker function

struct qwienerWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RVector<double> output;

  // initialize with source and destination
  qwienerWorker( const Rcpp::NumericMatrix input,
                 Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      output[j] = pwiener_scl( prm );
    }
  }
};

// Lookup - 06
//' @rdname dwiener
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qwiener( Rcpp::NumericVector p,
                             Rcpp::NumericVector ch,
                             Rcpp::NumericVector alpha,
                             Rcpp::NumericVector theta,
                             Rcpp::NumericVector xi,
                             Rcpp::NumericVector tau,
                             Rcpp::NumericVector sigma =
                               Rcpp::NumericVector::create(1.0),
                             bool joint = false,
                             double eps = 1e-29,
                             double bounds = 3.0,
                             double em_stop = 20, double err = 1e-8,
                             bool parYes = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    p.size(), ch.size(), alpha.size(), theta.size(),
    xi.size(), tau.size(), sigma.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 12 );

  // Fill matrix
  expanded_input = p[ create_index( n, p.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = ch[ create_index( n, ch.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = alpha[ create_index( n, alpha.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = theta[ create_index( n, theta.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(6) = expanded_input;
  double joint_val = 0.0; // Extract indicator for joint distribution
  if ( joint ) joint_val = 1.0;
  for (int k = 0; k < n; k++) {
    input(k,7) = eps;
    input(k,8) = joint_val;
    input(k,9) = bounds;
    input(k,10) = em_stop;
    input(k,11) = err;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate distribution function
  if ( !parYes ) {

    for (int j = 0; j < n; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      out(j) = qwiener_scl( prm );
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    qwienerWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  return( out );
}

// Lookup - 07
// RcppParallel worker function

struct rwienerWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RMatrix<double> output;

  // initialize with source and destination
  rwienerWorker( const Rcpp::NumericMatrix input,
                 Rcpp::NumericMatrix output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      std::vector<double> rd(2);
      rd = rwiener_scl( prm );

      output(j,0) = rd[0]; output(j,1) = rd[1];
    }
  }
};

// Lookup - 08
//' @rdname dwiener
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame rwiener( int n,
                         Rcpp::NumericVector alpha,
                         Rcpp::NumericVector theta,
                         Rcpp::NumericVector xi,
                         Rcpp::NumericVector tau,
                         Rcpp::NumericVector sigma =
                           Rcpp::NumericVector::create(1.0),
                         double eps = 1e-29,
                         double bounds = 5.0,
                         double em_stop = 30,
                         double err = 1e-16,
                         bool parYes = true ) {

  // Initialize output
  Rcpp::NumericMatrix out(n,2);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 12 );

  // Fill matrix
  expanded_input = Rcpp::runif(n,0.0,1.0);
  input.column(0) = expanded_input;
  expanded_input = Rcpp::runif(n,0.0,1.0);
  input.column(1) = expanded_input;
  expanded_input = alpha[ create_index( n, alpha.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = theta[ create_index( n, theta.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = xi[ create_index( n, xi.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = tau[ create_index( n, tau.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = sigma[ create_index( n, sigma.size() ) ];
  input.column(6) = expanded_input;
  for (int k = 0; k < n; k++) {
    input(k,7) = eps;
    input(k,8) = 0.0;
    input(k,9) = bounds;
    input(k,10) = em_stop;
    input(k,11) = err;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate distribution function
  if ( !parYes ) {

    for (int j = 0; j < n; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      std::vector<double> rd;
      rd = rwiener_scl( prm );

      out(j,0) = rd[0]; out(j,1) = rd[1];
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    rwienerWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  // Convert to a data-frame
  return Rcpp::DataFrame::create(
    Rcpp::_["rt"]= out.column(0),
    Rcpp::_["ch"]= out.column(1));
}
