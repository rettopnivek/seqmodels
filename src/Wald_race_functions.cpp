#include <RcppParallel.h>
#include "levyfunctions.h"
#include "sigfunctions.h"
#include "wrfunctions.h" // Scalar functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Vectorized functions for the random number generation, density,
distribution, and quantile functions for the Wald race model.

Index
Lookup - 01:  rwaldraceWorker
Lookup - 02:  rwaldrace
Lookup - 03:  dwaldraceWorker
Lookup - 04:  dwaldrace
Lookup - 05:  pwaldraceWorker
Lookup - 06:  pwaldrace
Lookup - 07:  dwaldraceWorker
Lookup - 08:  qwaldrace
*/

// Lookup - 01
// RcppParallel worker function

struct rwaldraceWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RMatrix<double> output;

  // initialize with source and destination
  rwaldraceWorker( const Rcpp::NumericMatrix input,
                 Rcpp::NumericMatrix output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(13);

      for (int i = 0; i < 13; i++) { prm[i] = input(j,i); }

      std::vector<double> rd;
      rd = rwaldrace_scl( prm );

      output(j,0) = rd[0]; output(j,1) = rd[1];
    }
  }
};

// Lookup - 02
//' The Wald Race Model
//'
//' Random generation, density, distribution, and quantile functions for
//' a two accumulator version of the Wald (or{ diffusion) race model
//' (Logan et al., 2014). For the racer representing the 1st choice,
//' \code{k1} is the threshold towards which evidence accumulates,
//' \code{x1} is the rate of evidence accumulation, \code{t1} is
//' the residual latency (the non-decision process), and \code{s1}
//' is the within-trial variability (coefficient of drift).
//' The racer representing the second choice has equivalent parameters
//' marked with '0'.
//'
//' @param n the number of draws for random generation.
//' @param rt a vector of responses times ( rt > 0 ).
//' @param p a vector of probabilities.
//' @param ch a vector of accuracy/choice values ( ch = {0,1} ).
//' @param k1 the threshold determining when a decision terminates for
//'   choices coded as 1 ( \code{k1} > 0).
//' @param x1 the average rate of evidence accumulation within a trial
//'   for choices coded as 1 (\code{x1} \eqn{\ge} 0).
//' @param t1 the residual latency for choices coded as 1
//'   (\code{t1} \eqn{\ge} 0).
//' @param k0 the threshold determining when a decision terminates for
//'   choices coded as 0 ( \code{k1} > 0).
//' @param x0 the average rate of evidence accumulation within a trial
//'   for choices coded as 0 (\code{x1} \eqn{\ge} 0).
//' @param t0 the residual latency for choices coded as 0
//'   (\code{t0} \eqn{\ge} 0).
//' @param s1 the within trial variability for choices coded as 1
//'   (\code{s1} > 0).
//' @param s0 the within trial variability for choices coded as 0
//'   (\code{s0} > 0).
//' @param rl logical; if \code{TRUE}, the residual latency shifts
//'   the finishing times for each racer before they are compared.
//'   Otherwise, the residual latency shifts the winning finishing
//'   times.
//' @param ln logical; if \code{TRUE}, probabilities are given as
//'   log(p).
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
//' The Wald (or diffusion) race model assumes that two independent one
//' boundary diffusion processes race each other. Whichever racer reaches
//' its threshold first determines the choice and response time. Because
//' of the independence, the likelihood for the Wald race model is:
//' \deqn{ f(t,y|\alpha)*(1-F(t,y|\beta),}
//' where \eqn{\alpha} and \eqn{\beta} are the sets of parameters for
//' the Wald distribution describing the finishing times for the
//' winning and losing racer respectively, and \eqn{f} and \eqn{F} refer
//' to the density and distribution functions respectively.
//'
//' The distribution function is estimated using numerical integration
//' via Gaussian quadrature.
//'
//' A linear interpolation approach is used to approximate the
//' quantile function, estimating the inverse of the cumulative
//' distribution function via an iterative procedure. When
//' the precision of this estimate is set to 8 decimal places,
//' the approximation will be typically accurate to about half of a
//' millisecond.
//'
//' @return
//' \code{dwaldrace} gives the density, \code{pwaldrace} approximates the
//' distribution function, \code{qwaldrace} approximates the quantile
//' function, and \code{rwaldrace} generates random deviates.
//'
//' The length of the result is determined by \code{n} for
//' \code{rwaldrace}, and is the maximum of the length of the numerical
//' arguments for the other functions.
//'
//' The numerical arguments other than \code{n} are recycled to the
//' length of the result.
//'
//' @section References:
//'
//' Logan, G. D., Van Zandt, T., Verbruggen, F., & Wagenmakers, E. J.
//'   (2014). On the ability to inhibit thought and action: General
//'   and special theories of an act of control. Psychological Review,
//'   121, 66.
//'
//' @examples
//' # Density function
//' dwaldrace( .5, c(1,0), k1 = 1.2, x1 = 2.0, t1 = .2,
//'   k0 = 1.0, x0 = 0.8, t0 = .2 )
//' # Distribution function
//' pwaldrace( .5, c(1,0), k1 = 1.2, x1 = 2.0, t1 = .2,
//'   k0 = 1.0, x0 = 0.8, t0 = .2 )
//' # Choice probabilities
//' pwaldrace( Inf, c(1,0), k1 = 1.2, x1 = 2.0, t1 = .2,
//'   k0 = 1.0, x0 = 0.8, t0 = .2 )
//' # Quantile function (Accurate to ~4 decimal places)
//' qwaldrace( .5, c(1,0), k1 = 1.2, x1 = 2.0, t1 = .2,
//'   k0 = 1.0, x0 = 0.8, t0 = .2 )
//'
//' # Simulate values
//' sim = rwiener( n = 100, alpha = 0.8, theta = 0.6,
//'   xi = 0.0, tau = 0.3 )
//'
//' # Plotting
//' layout( matrix( 1:4, 2, 2, byrow = T ) )
//' # Parameters
//' prm = c( k1 = 1.2, x1 = 2.0, t1 = .2, k0 = 1.0,
//'   x0 = 0.5, t0 = .2 )
//' # Density
//' obj = quickdist( 'wr', 'PDF', prm )
//' plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
//' # CDF
//' obj = quickdist( 'wr', 'CDF', prm )
//' plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
//' # Quantiles
//' obj = quickdist( 'wr', 'QF', prm, x = seq( .2, .8, .2 ) )
//' plot( obj ); prb = seq( .2, .8, .2 )
//' abline( h = prb, lty = 2 )
//' # Conditional, not joint
//' lines( obj, type = 'b', pch = 19, weight = 1 )
//' lines( obj, ch = 0, type = 'b', pch = 17, lty = 2, weight = 1 )
//' # Hazard function
//' obj = quickdist( 'wr', 'HF', prm )
//' plot( obj ); lines( obj ); lines( obj, ch = 0, lty = 2 )
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame rwaldrace( int n,
                           Rcpp::NumericVector k1,
                           Rcpp::NumericVector x1,
                           Rcpp::NumericVector t1,
                           Rcpp::NumericVector k0,
                           Rcpp::NumericVector x0,
                           Rcpp::NumericVector t0,
                           Rcpp::NumericVector s1 =
                             Rcpp::NumericVector::create(1.0),
                           Rcpp::NumericVector s0 =
                             Rcpp::NumericVector::create(1.0),
                           bool rl = false,
                           bool parYes = false ) {

  // Initialize output
  Rcpp::NumericMatrix out(n,2);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 13 );

  // Fill matrix
  expanded_input = k1[ create_index( n, k1.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = x1[ create_index( n, x1.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = t1[ create_index( n, t1.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = s1[ create_index( n, s1.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = k0[ create_index( n, k0.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = x0[ create_index( n, x0.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = t0[ create_index( n, t0.size() ) ];
  input.column(6) = expanded_input;
  expanded_input = s0[ create_index( n, s0.size() ) ];
  input.column(7) = expanded_input;
  // Generate random deviates
  expanded_input = Rcpp::rchisq(n,1);
  input.column(8) = expanded_input;
  expanded_input = Rcpp::runif(n,0.0,1.0);
  input.column(9) = expanded_input;
  expanded_input = Rcpp::rchisq(n,1);
  input.column(10) = expanded_input;
  expanded_input = Rcpp::runif(n,0.0,1.0);
  input.column(11) = expanded_input;
  for ( int k = 0; k < n; k++ ) {
    if (rl) {
      input(k,12) = 1.0;
    } else input(k,12) = 0.0;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate distribution function
  if ( !parYes ) {

    for (int j = 0; j < n; j++) {

      std::vector<double> prm(13);

      for (int i = 0; i < 13; i++) { prm[i] = input(j,i); }

      std::vector<double> rd;
      rd = rwaldrace_scl( prm );

      out(j,0) = rd[0]; out(j,1) = rd[1];
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    rwaldraceWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  // Convert to a data-frame
  return Rcpp::DataFrame::create(
    Rcpp::_["rt"]= out.column(0),
    Rcpp::_["ch"]= out.column(1));
}

// Lookup - 03
// RcppParallel worker function

struct dwaldraceWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RVector<double> output;

  // initialize with source and destination
  dwaldraceWorker( const Rcpp::NumericMatrix input,
                 Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(11);

      for (int i = 0; i < 11; i++)
        prm[i] = input(j,i);

      output[j] = dwaldrace_scl( prm );
    }
  }
};

// Lookup - 04
//' @rdname rwaldrace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dwaldrace( Rcpp::NumericVector rt,
                               Rcpp::NumericVector ch,
                               Rcpp::NumericVector k1,
                               Rcpp::NumericVector x1,
                               Rcpp::NumericVector t1,
                               Rcpp::NumericVector k0,
                               Rcpp::NumericVector x0,
                               Rcpp::NumericVector t0,
                               Rcpp::NumericVector s1 =
                                 Rcpp::NumericVector::create(1.0),
                               Rcpp::NumericVector s0 =
                                 Rcpp::NumericVector::create(1.0),
                               bool rl = false,
                               bool ln = false,
                               bool parYes = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    rt.size(), ch.size(), k1.size(), x1.size(),
    t1.size(), s1.size(), k0.size(), x0.size(),
    t0.size(), s0.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 11 );

  // Fill matrix
  expanded_input = rt[ create_index( n, rt.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = ch[ create_index( n, ch.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = k1[ create_index( n, k1.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = x1[ create_index( n, x1.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = t1[ create_index( n, t1.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = s1[ create_index( n, s1.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = k0[ create_index( n, k0.size() ) ];
  input.column(6) = expanded_input;
  expanded_input = x0[ create_index( n, x0.size() ) ];
  input.column(7) = expanded_input;
  expanded_input = t0[ create_index( n, t0.size() ) ];
  input.column(8) = expanded_input;
  expanded_input = s0[ create_index( n, s0.size() ) ];
  input.column(9) = expanded_input;
  for ( int k = 0; k < n; k++ ) {
    if (rl) {
      input(k,10) = 1.0;
    } else input(k,10) = 0.0;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate likelihood
  if ( !parYes ) {

    for (int j = 0; j < n; j++) {

      std::vector<double> prm(11);

      for (int i = 0; i < 11; i++)
        prm[i] = input(j,i);

      out(j) = dwaldrace_scl( prm );
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    dwaldraceWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  // If indicated, compute the standard density
  if ( !ln ) out = exp(out);

  return( out );
}

// Lookup - 05
// RcppParallel worker function

struct pwaldraceWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RVector<double> output;

  // initialize with source and destination
  pwaldraceWorker( const Rcpp::NumericMatrix input,
                   Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(11);

      for (int i = 0; i < 11; i++)
        prm[i] = input(j,i);

      output[j] = pwaldrace_scl( prm );
    }
  }
};


// Lookup - 06
//' @rdname rwaldrace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pwaldrace( Rcpp::NumericVector rt,
                               Rcpp::NumericVector ch,
                               Rcpp::NumericVector k1,
                               Rcpp::NumericVector x1,
                               Rcpp::NumericVector t1,
                               Rcpp::NumericVector k0,
                               Rcpp::NumericVector x0,
                               Rcpp::NumericVector t0,
                               Rcpp::NumericVector s1 =
                                 Rcpp::NumericVector::create(1.0),
                               Rcpp::NumericVector s0 =
                                 Rcpp::NumericVector::create(1.0),
                               bool rl = false,
                               bool ln = false,
                               bool lower_tail = true,
                               bool parYes = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    rt.size(), ch.size(), k1.size(), x1.size(),
    t1.size(), s1.size(), k0.size(), x0.size(),
    t0.size(), s0.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 11 );

  // Fill matrix
  expanded_input = rt[ create_index( n, rt.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = ch[ create_index( n, ch.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = k1[ create_index( n, k1.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = x1[ create_index( n, x1.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = t1[ create_index( n, t1.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = s1[ create_index( n, s1.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = k0[ create_index( n, k0.size() ) ];
  input.column(6) = expanded_input;
  expanded_input = x0[ create_index( n, x0.size() ) ];
  input.column(7) = expanded_input;
  expanded_input = t0[ create_index( n, t0.size() ) ];
  input.column(8) = expanded_input;
  expanded_input = s0[ create_index( n, s0.size() ) ];
  input.column(9) = expanded_input;
  for ( int k = 0; k < n; k++ ) {
    if (rl) {
      input(k,10) = 1.0;
    } else input(k,10) = 0.0;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Calculate CDF
  if ( !parYes ) {
    for (int j = 0; j < n; j++) {

      std::vector<double> prm(11);

      for (int i = 0; i < 11; i++)
        prm[i] = input(j,i);

      out(j) = pwaldrace_scl( prm );
    }
  } else {
    // Function call operator that works for the specified
    // range (begin/end)
    pwaldraceWorker mt(input, out);

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

// Lookup - 07
// RcppParallel worker function

struct qwaldraceWorker : public RcppParallel::Worker
{
  // Input matrix
  const RcppParallel::RMatrix<double> input;

  // Destination matrix
  RcppParallel::RVector<double> output;

  // initialize with source and destination
  qwaldraceWorker( const Rcpp::NumericMatrix input,
                   Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(14);

      for (int i = 0; i < 14; i++)
        prm[i] = input(j,i);

      output[j] = qwaldrace_scl( prm );
    }
  }
};

// Lookup - 08
//' @rdname rwaldrace
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qwaldrace( Rcpp::NumericVector p,
                               Rcpp::NumericVector ch,
                               Rcpp::NumericVector k1,
                               Rcpp::NumericVector x1,
                               Rcpp::NumericVector t1,
                               Rcpp::NumericVector k0,
                               Rcpp::NumericVector x0,
                               Rcpp::NumericVector t0,
                               Rcpp::NumericVector s1 =
                                 Rcpp::NumericVector::create(1.0),
                               Rcpp::NumericVector s0 =
                                 Rcpp::NumericVector::create(1.0),
                               bool rl = false,
                               double bounds = 3.0,
                               double em_stop = 20,
                               double err = 1e-8,
                               bool parYes = true ) {

  // Determine the longest input vector
  int n = max( Rcpp::NumericVector::create(
    p.size(), ch.size(), k1.size(), x1.size(),
    t1.size(), s1.size(), k0.size(), x0.size(),
    t0.size(), s0.size() ) );

  // Initialize output
  Rcpp::NumericVector out(n);

  // Create matrix whose rows are inputs
  // to the scalar function
  Rcpp::NumericVector expanded_input(n);
  Rcpp::NumericMatrix input( n, 14 );

  // Fill matrix
  expanded_input = p[ create_index( n, p.size() ) ];
  input.column(0) = expanded_input;
  expanded_input = ch[ create_index( n, ch.size() ) ];
  input.column(1) = expanded_input;
  expanded_input = k1[ create_index( n, k1.size() ) ];
  input.column(2) = expanded_input;
  expanded_input = x1[ create_index( n, x1.size() ) ];
  input.column(3) = expanded_input;
  expanded_input = t1[ create_index( n, t1.size() ) ];
  input.column(4) = expanded_input;
  expanded_input = s1[ create_index( n, s1.size() ) ];
  input.column(5) = expanded_input;
  expanded_input = k0[ create_index( n, k0.size() ) ];
  input.column(6) = expanded_input;
  expanded_input = x0[ create_index( n, x0.size() ) ];
  input.column(7) = expanded_input;
  expanded_input = t0[ create_index( n, t0.size() ) ];
  input.column(8) = expanded_input;
  expanded_input = s0[ create_index( n, s0.size() ) ];
  input.column(9) = expanded_input;
  for (int k = 0; k < n; k++) {
    if (rl) {
      input(k,10) = 1.0;
    } else input(k,10) = 0.0;
    input(k,11) = bounds;
    input(k,12) = em_stop;
    input(k,13) = err;
  }

  // If number of observations <= 16, no parallelization
  if (n <= 16) parYes = false;

  // Estimate quantiles
  if ( !parYes ) {
    for (int j = 0; j < n; j++) {

      std::vector<double> prm(14);

      for (int i = 0; i < 14; i++)
        prm[i] = input(j,i);

      out(j) = qwaldrace_scl( prm );
    }
  } else {
    // Function call operator that works for the specified
    // range (begin/end)
    qwaldraceWorker mt(input, out);

    // Call parallelFor to do the work
    RcppParallel::parallelFor(0, n, mt);
  }

  return( out );
}
