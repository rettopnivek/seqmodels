// Include guard to protect against multiple definitions error
#ifndef __WIENERFUNCTIONS__
#define __WIENERFUNCTIONS__

#include <Rcpp.h>
#include "wpparamverify.h"
#include "blurtonetal.h"
#include "navarrofuss.h"
#include "miscfunctions.h"

/*
Purpose:
Scalar functions for the density, distribution, quantile,
and random generation functions of the two boundary
wiener process. The quantile and random generation functions
rely on linear interpolation to approximate the inverse of
the cumulative probabilty function.

Index:
Lookup - 01:  dwiener_scl
Lookup - 02:  pwiener_scl
Lookup - 03:  qwiener_scl
Lookup - 04:  rwiener_scl
*/

// Lookup - 01

inline double dwiener_scl( std::vector<double> par ) {
  /*
  Purpose:
  A scalar version to calculate the full version of the
  density for the Wiener process.
  Arguments:
  par - A vector of parameters, where...
  par[0] = response time
  par[1] = choice (0 or 1)
  par[2] = boundary separation
  par[3] = starting point proportion
  par[4] = drift rate
  par[5] = residual latency
  par[6] = coefficient of drift
  par[7] = precision
  par[8] = Indicator for whether joint or conditional
  Returns:
  The density (joint or conditional) for the two-boundary wiener process.
  */

  // Initialize output
  double out = 0.0;

  // Create index
  std::vector<int> index = create_range( 1, 7 );

  // Check for valid inputs
  if ( wp_param_verify( par, index ) ) {

    double rt = par[0]; double ch = par[1];
    double alpha = par[2]; double theta = par[3];
    double xi = par[4]; double tau = par[5];
    double sigma = par[6]; double eps = par[7];
    int joint = par[8];

    // Decision time
    double dt = rt - tau;

    // Determine starting point
    double zeta = alpha*theta;

    // Rescale paramters based on within-trial variance
    alpha = alpha/sigma;
    xi = xi/sigma;
    zeta = zeta/sigma;

    // Adjust parameters based on absorbing boundary
    if (ch == 1.0) {
      zeta = alpha-zeta;
      xi = -1.0*xi;
    }

    out = wfpt( dt, xi, alpha, zeta, eps );

    // If conditional density is requested
    if ( joint == 0 ) {
      // Probablity of reaching specified boundary
      double adj = pr_absorb( alpha, zeta, xi, sigma );
      out = out/adj; // Re-normalize
    }
  }

  return( out );
}

// Lookup - 02

inline double pwiener_scl( std::vector<double> par ) {
  /*
  Purpose:
  A scalar versions to calculate the full version of the
  distribution function for the wiener process.
  Arguments:
  par - A vector of parameters, where...
  par[0] = response time
  par[1] = choice (0 or 1)
  par[2] = boundary separation
  par[3] = starting point proportion
  par[4] = drift rate
  par[5] = residual latency
  par[6] = coefficient of drift
  par[7] = precision
  par[8] = Indicator for whether joint or conditional
  Returns:
  The value for the distribution function.
  */

  // Initialize output
  double out = 0.0;

  // Create index
  std::vector<int> index = create_range( 1, 9 );
  index[7] = 0; // Skip final inputs
  index[8] = 0;
  // Check for valid inputs
  if ( wp_param_verify( par, index ) ) {

    double rt = par[0]; double ch = par[1];
    double alpha = par[2]; double theta = par[3];
    double xi = par[4]; double tau = par[5];
    double sigma = par[6]; double eps = par[7];
    double joint = par[8];

    // Adjust parameters based on absorbing boundary
    if (ch == 1.0) {
      theta = 1.0 - theta;
      xi = -xi;
    }

    // If time is set to Inf, return the probability of absorption
    // at the relevant boundary
    if (rt == R_PosInf) out = Pu( xi, alpha, theta);

    // Decision time
    double dt = rt - tau;

    // Calculate CDF
    out = F_lower( dt, xi, alpha, theta, sigma, eps );

    // If conditional distribution is requested
    if ( joint == 0.0 ) {
      // Probability of reaching specified boundary
      double adj = Pu( xi, alpha, theta);
      out = out/adj; // Re-normalize
    }
  }

  if (out < 0.0) out = 0.0;

  return( out );
}

// Lookup - 03

inline double qwiener_scl( std::vector<double> par ) {
  /*
  Purpose:
  Scalar function for quantile function of the two boundary
  wiener process.
  Arguments:
  par - A numeric vector of inputs, where...
  par[0]  = A probability
  par[1]  = choice (0 or 1)
  par[2]  = boundary separation
  par[3]  = starting point proportion
  par[4]  = drift rate
  par[5]  = residual latency
  par[6]  = coefficient of drift
  par[7]  = precision
  par[8]  = Indicator for whether joint or conditional
  par[9]  = Maximum quantile to explore
  par[10] = Maximum number of iterations to attempt
  par[11] = Precision for inverse CDF approximation
  Returns:
  An estimated quantile.
  */

  // Initialize output
  double out = NA_REAL;

  // Create index
  std::vector<int> index = create_range( 1, 12 );
  index[0] = 8; // First input is a probability
  // Skip nuisance parameters
  for ( int s = 7; s < 12; s++ ) index[s] = 0;
  // Check for valid inputs
  if ( wp_param_verify( par, index ) ) {

    cdf cur_cdf = pwiener_scl;

    // Extract values governing linear interpolation
    double mn = par[5];
    double mx = par[9];
    int em_stop = par[10];
    double err = par[11];

    out = lin_interp_quantile( par, cur_cdf,
                               mn, mx, em_stop, err,
                               par[5], R_PosInf );
  }

  return( out );
}

// Lookup - 04

inline std::vector<double> rwiener_scl( std::vector<double> par ) {
  /*
  Purpose:
  Scalar function for generating random deviates from a
  two boundary wiener process.
  Arguments:
  par - A numeric vector of inputs, where...
  par[0]  = A random probability
  par[1]  = A random probability
  par[2]  = boundary separation
  par[3]  = starting point proportion
  par[4]  = drift rate
  par[5]  = residual latency
  par[6]  = coefficient of drift
  par[7]  = precision
  par[8] = Indicator for whether joint or conditional
  par[9]  = Maximum quantile to explore
  par[10]  = Maximum number of iterations to attempt
  par[11] = Precision for inverse CDF approximation
  Returns:
  A vector with a response time and a choice.
  */

  // Initialize output
  std::vector<double> out(2);
  out[0] = NA_REAL; out[1] = NA_REAL;

  // Create index
  std::vector<int> index = create_range( 1, 12 );
  index[0] = 8; index[1] = 8; // First two inputs are probabilities
  // Skip nuisance parameters
  for ( int s = 7; s < 12; s++ ) index[s] = 0;
  // Check for valid inputs
  if ( wp_param_verify( par, index ) ) {

    // Determine choice
    double ch = 1.0;
    double p_upper = Pu( par[4], par[2], par[3] );
    if ( p_upper > par[1] ) ch = 0.0;
    par[1] = ch; // Save choice
    out[1] = ch;

    cdf cur_cdf = pwiener_scl;

    // Extract values governing linear interpolation
    double mn = par[5];
    double mx = par[9];
    int em_stop = par[10];
    double err = par[11];

    out[0] = lin_interp_quantile( par, cur_cdf,
                                  mn, mx, em_stop, err,
                                  par[5], R_PosInf );
  }

  return( out );
}

#endif // End include guard
