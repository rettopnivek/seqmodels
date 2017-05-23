// Include guard to protect against multiple definitions error
#ifndef __EWFUNCTIONS__
#define __EWFUNCTIONS__

#include <Rcpp.h> // Includes certain libraries of functions
#include "levyfunctions.h" // For drift of 0
#include "sigfunctions.h" // For inverse gaussian
#include "miscfunctions.h" // Linear interpolation
#include "gsl/include/gsl_integration.h" // Numerical integration
#include "gsl/include/gsl_errno.h" // Error handling

/*
Purpose:
Scalar functions for the random generation, density, distribution,
and quantile functions of the convolution of the inverse
Gaussian (Wald) and exponential distributions.

References:
Abramowitz, M., & Stegun, I. A. (Eds.) (1964). Handbook of
  Mathematical Functions with Formulas, Graphs, and Mathematical
  Tables. Washington D. C.: United States Department of Commerce,
  National Bureau of Standards.
Heathcote, A. (2004a). Fitting Wald and ex-Wald distributions to
  response time data: An example using functions for the S-PLUS
  package. Behavior Research Methods Instruments & Computers, 36,
  678 - 694.
Heathcote, A. (2004b). rtfit.ssc. Retrieved May 5, 2017 from
  Psychonomic Society Web Archive: http://www.psychonomic.org/ARCHIVE/.
Schwarz, W. (2002). On the convolution of the inverse Gaussian and
  the exponential random variables. Communications in Statistics,
  Theory and Methods, 31, 2113 - 2121.
Index
Lookup - 01:  ew_param_verify
Lookup - 02:  uandv
Lookup - 03:  rew
Lookup - 04:  rexwald_scl
Lookup - 05:  dexwald_scl
*/

// Lookup - 01

inline bool ew_param_verify( std::vector<double> prm,
                             std::vector<int> index ) {
  /*
  Purpose:
  Function to verify whether inputs are admissable for the
  ex-Wald distribution.
  Arguments:
  prm   - A vector of parameters
  index - The type of parameter being inputted, where...
  0 = skip
  1 = time
  2 = kappa (Threshold)
  3 = xi (Drift rate)
  4 = tau (Mean of residual latency)
  5 = sigma (Coefficient of drift)
  6 = probability
  Returns:
  TRUE if the inputs are admissable, FALSE otherwise.
  */

  // Initialize output
  bool out;

  // Determine size of vector
  int sz = index.size();

  // Initialize variable for checks
  int chk = 0;
  // Variable for number of checks that need to be passed
  int n_pass = 0;

  for ( int i = 0; i < sz; i++ ) {

    if ( index[i] > 0 ) n_pass += 1;

    // Check whether inputted time is a real number and
    // greater than 0
    if ( index[i] == 1 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether threshold is a real number and whether it
    // is greater than 0
    if ( index[i] == 2 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether drift rate is a real number and whether it
    // is greater than or equal to 0
    if ( index[i] == 3 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] >= 0.0 ) ) chk += 1;
    }

    // Check whether 1/rate parameter is a real number and whether it
    // is greater than 0
    if ( index[i] == 4 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether coefficient of drift is a real number and whether
    // it is greater than 0
    if ( index[i] == 5 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether probability is a real number and whether it
    // is appropriately bounded
    if ( index[i] == 6 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( (prm[i] >= 0.0) | (prm[i] <= 1.0) ) ) chk += 1;
    }

  }

  // Test how many checks were passed
  out = chk == n_pass;

  return( out );
}


// Lookup - 02

inline std::vector<double> uandv( double x, double y,
                                  int firstblock ) {
  /*
  Purpose:
  A function to compute the series approximatin to the real(u)
  and imaginary (v) parts of the complex error function erf( x + iy).
  See Abramowitz and Stegun (1964), 7.1.29. The code is adapted
  from the S-PLUS script of Heathcote (2004).
  Arguments:
  x          - A continuous value
  y          - A continuous value
  firstblock - The number of blocks to use in the series
               approximation
  Returns:
  A vector with the real and imaginary parts of the complex
  error function.
  */

  // Part 1
  double twoxy = 2.0 * x * y;
  double xsq = pow( x, 2.0 );
  double iexpxsqpi = 1.0 / ( M_PI * exp(xsq) );

  // Part 2
  double sin2xy = sin(twoxy);
  double cos2xy = cos(twoxy);

  // Part 3 - 5

  // Vector declarations
  std::vector<double> nmat( firstblock );
  std::vector<double> nsqmat( firstblock );
  std::vector<double> ny( firstblock );
  std::vector<double> twoxcoshny( firstblock );
  std::vector<double> nsinhny( firstblock );
  std::vector<double> nsqfrac( firstblock );
  double sm_1 = 0.0; // Define terms for summations
  double sm_2 = 0.0;

  // Loop
  for ( int i = 0; i < firstblock; i++ ) {
    nmat[i] = i + 1;
    nsqmat[i] = nmat[i] * nmat[i];
    ny[i] = nmat[i] * y;
    twoxcoshny[i] = 2.0 * x * cosh(ny[i]);
    nsinhny[i] = nmat[i] * sinh( ny[i] );
    nsqfrac[i] = ( exp( -1.0 * nsqmat[i]/4.0 ) /
      ( nsqmat[i] + 4.0 * xsq ) );
    // Take dot product
    sm_1 += nsqfrac[i] * ( 2.0 * x - twoxcoshny[i] * cos2xy +
      nsinhny[i] * sin2xy );
    sm_2 += nsqfrac[i] *
      ( twoxcoshny[i] * sin2xy + nsinhny[i] * cos2xy );
  }

  // Part 6
  double u;
  u = 2.0 * R::pnorm( x * sqrt( 2.0 ), 0.0, 1.0, 1, 0 ) - 1.0;
  u += iexpxsqpi*( ( ( 1.0 - cos2xy )/( 2.0 * x ) ) + 2.0*sm_1 );

  // Part 7
  double v;
  v = iexpxsqpi * ( ( sin2xy/( 2.0 * x ) ) + 2.0 * sm_2 );

  // Return real and imaginary parts
  std::vector<double> out(2);
  out[0] = u;
  out[1] = v;

  return(out);
}

// Lookup - 03

inline double rew( double x, double y ) {
  /*
  Purpose:
  Computes the real part of the term in the ex-Wald density
  that involves the complex error function. The code is
  taken from the S-PLUS script of Heathcote (2004).
  Arguments:
  x - A continuous value
  y - A continuous value
  Notes:
  The function assumes the series approximation will converge
  within 20 blocks.
  Returns:
  The real part of the term in the ex-Wald density that involves
  the complex error function.
  */

  std::vector<double> uv = uandv( y, x, 20 );

  double out = exp( pow( y, 2.0 ) - pow( x, 2.0 ) ) *
    ( cos( 2.0 * x * y ) * ( 1.0 - uv[0] ) +
    sin( 2.0 * x * y ) * uv[1] );

  return( out );
}

// Lookup - 04

inline double conv_dig_dexp( double x, void * params) {
  /*
  Purpose:
  A scalar function that can be integrated over in order to
  compute the convolution of the inverse gaussian and the
  exponential distributions.
  Arguments:
  x      - A residual latency
  *param - A pointer to a vector of parameters, where...
  param[0] = Response time
  param[1] = Threshold
  param[2] = Drift rate
  param[3] = 1/Rate
  param[4] = The coefficient of drift
  Returns:
  The convolution of the inverse gaussian and exponential
  distributions.
  */

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;
  out = sig_to_ig( par[0] - x, par[1], par[2], par[4], 0, 1 ) *
    ( exp( -x/par[3] )/par[3] );

  return out;
}

// Lookup - 05

inline double dexwald_ni( std::vector<double> par ) {
  /*
  Purpose:
  ...
  Arguments:
  par - ...
  Returns:
  ...
  */

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  // Initialize result
  double result = 0.0;
  double a = 0.0;

  // Allocate memory
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

  double error;

  gsl_function F;
  F.function = &conv_dig_dexp;
  F.params = &par;

  // GSL QAGSIU parameters:
  // function, lower, absolute error,
  // relative error, subinterval limit,
  // workspace, result, estimated error
  gsl_integration_qagiu (&F, a, 0, 1e-10, 1000,
                         w, &result, &error);

  if (result < 0.0) result = 0.0;

  return( result );
}


// Lookup - 06

inline double rexwald_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function to generate random deviates from
  the ex-Wald distribution.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = A threshold
  prm[1] = A drift rate
  prm[2] = The mean of the exponential component
  prm[3] = The coefficient of drift
  prm[4] = A chi-square random deviate with 1 degree of freedom
  prm[5] = A unit uniform deviate
  prm[6] = A unit uniform deviate
  Returns:
  A random deviate.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 2, 5 );

  // Check for valid inputs
  if ( ew_param_verify( prm, index ) ) {

    // Compute random deviate from inverse gausian
    double ig_rd = sig_to_ig( prm[4], prm[1], prm[2], prm[4],
                              prm[5], 3 );

    // Compute random deviate from exponential
    double e_rd = -log( 1.0 - prm[6] ) * prm[3];

    // Output is sum of the two random deviates
    out = ig_rd + e_rd;
  }

  return( out );
}

// Lookup - 07

inline double dexwald_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function to compute the density of the ex-Wald
  distribution. Code is taken from the S-PLUS script of
  Heathcote (2004).
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = A response time
  prm[1] = A threshold
  prm[2] = A drift rate
  prm[3] = The mean for the exponential component
  prm[4] = The coefficient of drift
  prm[5] = An indicator for the log-likelihood
  prm[6] = An indicator for numerical integration
  Returns:
  The likelihood or log-likelihood.
  */

  // Initialize output
  double out = log( 0.0 );

  // Determine if log-density should be returned
  int ln = prm[5];

  // Define index
  std::vector<int> index = create_range( 1, 5 );

  // Check for valid inputs
  if ( ew_param_verify( prm, index ) ) {

    // Extract parameters
    double rt = prm[0];
    double kappa = prm[1];
    double xi = prm[2];
    double tau = prm[3];
    double sigma = prm[4];

    // If mean for residual latencies is
    // larger than 2 ms
    if ( (tau > 0.002) && (prm[6] == 0.0) ) {

      // Convert mean residual latency to rate
      double gam = 1.0/tau;
      // Square the coefficient of drift
      double s2 = pow( sigma, 2.0 );

      double k = pow( xi, 2.0 ) - 2.0 * gam * s2;

      // Check if xi-squared is greater than
      // 2 x 1/tau x sigma-sqaured
      if ( k > 0.0 ) {

        // Compute equation 18 from Schwarz (2002)
        k = sqrt( k );
        double term1 = gam * exp( -gam * rt + xi * kappa / s2 );
        double term2 = exp( -k * kappa / s2 );
        double Phi1 = R::pnorm(
          ( k * rt - kappa )/( sigma * sqrt( rt ) ),
          0.0, 1.0, 1, 0 );
        double term3 = exp( k * kappa / s2 );
        double Phi2 = R::pnorm(
          -( k * rt + kappa )/( sigma * sqrt( rt ) ),
          0.0, 1.0, 1, 0 );
        out = log( term1) + log( term2 * Phi1 + term3 * Phi2 );

      } else {
        out = xi * kappa - pow( kappa, 2.0 )/( 2.0 * rt ) -
          rt * pow( xi, 2.0)/2.0;
        out += log( rew( sqrt( -rt * k / 2.0 ),
                         kappa/sqrt( 2.0 * rt ) ) / tau );
      }
    } else {
      // Use numerical integration for accuracy
      out = log( dexwald_ni( prm ) );
    }
  }

  // If specified return the density
  if (ln==0.0) out = exp(out);

  return( out );
}

#endif // End include guard
