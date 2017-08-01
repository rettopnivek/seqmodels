// Include guard to protect against multiple definitions error
#ifndef __EWRFUNCTIONS__
#define __EWRFUNCTIONS__

#include <Rcpp.h> // Includes certain libraries of functions
#include "levyfunctions.h" // For drift of 0
#include "sigfunctions.h" // Wald distribution
#include "complex_error_function.h"
#include "ewfunctions.h"
#include "miscfunctions.h" // Linear interpolation
#include "gsl/include/gsl_integration.h" // Numerical integration
#include "gsl/include/gsl_errno.h" // Error handling

/*
Purpose:
Scalar functions for the random generation, density, distribution,
and quantile functions of the ex-Wald race model.

Index
Lookup - 01:  ewr_param_verify
Lookup - 02:  rexwaldrace_scl
Lookup - 03:  dexwaldrace_scl
Lookup - 04:  int_dexwaldrace_scl
Lookup - 05:  pexwaldrace_scl
Lookup - 06:  qexwaldrace_scl
*/

// Lookup - 01

inline bool ewr_param_verify( std::vector<double> prm,
                             std::vector<int> index ) {
  /*
  Purpose:
  Function to verify whether inputs are admissable for the
  ex-Wald race model.
  Arguments:
  prm   - A vector of parameters
  index - The type of parameter being inputted, where...
  0 = skip
  1 = time
  2 = choice
  3 = kappa (Threshold)
  4 = xi (Drift rate)
  5 = tau (Shift)
  6 = sigma (Coefficient of drift)
  7 = probability
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

    // Check whether inputted response time is a real number
    // and greater than 0
    if ( index[i] == 1 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) {
        chk += 1;
        rt_yes = prm[i];
      }
    }

    // Check whether inputted choice is a real number and
    // whether it is 0 or 1
    if ( index[i] == 2 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( ( prm[i] == 0.0 ) || ( prm[i] == 1.0 ) ) ) {
        chk += 1;
        ch_yes = prm[i];
      }
    }

    // Check whether threshold is a real number and whether it
    // is greater than 0
    if ( index[i] == 3 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether drift rate is a real number and whether it
    // is greater than or equal to 0
    if ( index[i] == 4 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] >= 0.0 ) ) chk += 1;
    }

    // Check whether residual latency is a real number and
    // whether it is greater than or equal to 0
    if ( index[i] == 5 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] >= 0.0 ) ) chk += 1;
    }

    // Check whether coefficient of drift is a real number and whether
    // it is greater than 0
    if ( index[i] == 6 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether probability is a real number and whether it
    // is appropriately bounded
    if ( index[i] == 7 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( (prm[i] >= 0.0) | (prm[i] <= 1.0) ) ) chk += 1;
    }

  }

  // Test how many checks were passed
  out = chk == n_pass;

  return( out );
}

// Lookup - 02

inline std::vector<double> rexwaldrace_scl( std::vector<double> prm ) {
  /*
  Purpose:
  A scalar function to generate random deviates from the Wald race
  model.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = Threshold (1)
  prm[1] = Drift rate (1)
  prm[2] = Residual latency (1)
  prm[3] = Coefficient of drift (1)
  prm[4] = Threshold (0)
  prm[5] = Drift rate (0)
  prm[6] = Residual latency (0)
  prm[7] = Coefficient of drift (0)
  prm[8] = A chi-square random deviate with 1 degree of freedom
  prm[9] = A unit uniform deviate
  prm[10] = A unit uniform deviate
  prm[11] = A chi-square random deviate with 1 degree of freedom
  prm[12] = A unit uniform deviate
  prm[13] = A unit uniform deviate
  prm[14] = Indicator for whether residual latency impacts decision
  Returns:
  A random response time and choice.
  */

  // Initialize output
  std::vector<double> out(2);
  out[0] = NA_REAL; out[1] = NA_REAL;

  // Create index
  std::vector<int> index = create_range( 3, 10 );
  for (int i = 4; i < 8; i++ ) index[i] = i - 1;

  // Check for valid inputs
  if ( wr_param_verify( prm, index ) ) {

    // Vector for inverse Gaussian
    std::vector<double> prm_sig(6);

    // Extract parameters for 1st racer
    for ( int i = 0; i < 4; i++ ) prm_sig[i] = prm[i];
    prm_sig[2] = 0.0; prm_sig[4] = prm[8]; prm_sig[5] = prm[9];
    // Simulate response time
    double t1 = rinvgauss_scl(prm_sig);
    // Simulate residual latency
    double tau_1 = -log( 1.0 - prm[10] ) * prm[2];

    // Extract parameters for 2nd racer
    for ( int i = 0; i < 4; i++ ) prm_sig[i] = prm[i+4];
    prm_sig[2] = 0.0; prm_sig[4] = prm[11]; prm_sig[5] = prm[12];
    // Simulate response time
    double t0 = rinvgauss_scl(prm_sig);
    // Simulate residual latency
    double tau_0 = -log( 1.0 - prm[13] ) * prm[6];

    // Extract indicator
    int rl = prm[14];

    // If the residual latency impacts the decision rule
    if ( rl == 1 ) {
      t1 = t1 + tau_1;
      t0 = t0 + tau_0;
    }

    // Determine the winning accumulator
    if (t1 < t0) {
      out[0] = t1;
      if ( rl == 0 ) out[0] += tau_1;
      out[1] = 1;
    }
    if (t1 > t0) {
      out[0] = t0;
      if ( rl == 0 ) out[0] += tau_0;
      out[1] = 0;
    }

  }

  return( out );
}

// Lookup - 03

inline double dexwaldrace_scl( std::vector<double> prm ) {
  /*
  Purpose:
  A scalar function to the density for the ex-Wald race model.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = response time
  prm[1] = choice (0 or 1)
  prm[2] = Threshold (1)
  prm[3] = Drift rate (1)
  prm[4] = Residual latency (1)
  prm[5] = Coefficient of drift (1)
  prm[6] = Threshold (0)
  prm[7] = Drift rate (0)
  prm[8] = Residual latency (0)
  prm[9] = Coefficient of drift (0)
  prm[10] = Indicator for whether residual latency impacts decision
  prm[11] = Indicator on whether to use numerical integration
  Returns:
  The joint log-density for the ex-Wald race model.
  */

  // Initialize output
  double out = log( 0.0 );

  // Create index
  std::vector<int> index = create_range( 1, 10 );
  for (int i = 6; i < 10; i++ ) index[i] = i - 3;

  // Check for valid inputs
  if ( ewr_param_verify( prm, index ) ) {

    // Extract data
    double rt = prm[0]; double ch = prm[1];
    double rl = prm[10];

    std::vector<double> f(7); // PDF for winning racer
    std::vector<double> F(6); // CDF for losing racer
    f[0] = rt; F[0] = rt; f[5] = 1.0;
    f[6] = prm[11]; F[5] = prm[11];

    if ( ch == 1.0 ) {
      for ( int i = 0; i < 4; i++ ) {
        f[i + 1] = prm[ i + 2 ];
        F[i + 1] = prm[ i + 6 ];
      }
    } else {
      for ( int i = 0; i < 4; i++ ) {
        f[i + 1] = prm[ i + 6 ];
        F[i + 1] = prm[ i + 2 ];
      }
    }

    if ( rl == 0.0 ) {
      if ( ch == 1.0 )
        F[3] = prm[4];
      else
        F[3] = prm[8];
    }

    // Compute log-likelihood
    out = dexwald_scl( f ) + log( 1.0 - pexwald_scl( F ) );
  }

  return( out );
}

// Lookup - 04

inline double int_dexwaldrace_scl( double x, void * params) {
  /*
  Purpose:
  A scalar function for the Wald race model density that can be
  numerically integrated using GSL routines.
  Arguments:
  x      - A response time
  *param - A pointer to a vector of parameters, where...
   param[0] = response time
   param[1] = choice (0 or 1)
   param[2] = Threshold (1)
   param[3] = Drift rate (1)
   param[4] = Residual latency (1)
   param[5] = Coefficient of drift (1)
   param[6] = Threshold (0)
   param[7] = Drift rate (0)
   param[8] = Residual latency (0)
   param[9] = Coefficient of drift (0)
   param[10] = Indicator for whether residual latency impacts
               decision
   param[10] = Indicator for whether to use numerical integration
  Returns:
  The joint density for the ex-Wald race model.
  */

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  par[0] = x;

  // Calculate the density for the Wald race model
  out = exp( dexwaldrace_scl( par ) );

  return out;
}

// Lookup - 05

inline double pexwaldrace_scl( std::vector<double> par ) {
  /*
  Purpose:
  A scalar function that computes the distribution function
  for the Wald race model by numerically integrating over
  the density function.
  Arguments:
  par[0] = Response time
  par[1] = Choice (0 or 1)
  par[2] = Threshold (1)
  par[3] = Drift rate (1)
  par[4] = Residual latency (1)
  par[5] = Coefficient of drift (1)
  par[6] = Threshold (0)
  par[7] = Drift rate (0)
  par[8] = Residual latency (0)
  par[9] = Coefficient of drift (0)
  par[10] = Indicator for whether residual latency
            impacts decision
  par[11] = Indicator for whether to use numerical integration
  Returns:
  The cumulative probability.
  */

  double a = 0.0;
  double b = par[0];

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( a >= 0.0 ) {

    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &int_dexwaldrace_scl;
    F.params = &par;

    // If the upper boundary is infinite
    if ( b == R_PosInf ) {
      // GSL QAGSIU parameters:
      // function, lower, absolute error,
      // relative error, subinterval limit,
      // workspace, result, estimated error
      gsl_integration_qagiu (&F, a, 0, 1e-10, 1000,
                             w, &result, &error);
    } else {
      // GSL QAGS parameters:
      // function, lower, upper, absolute error,
      // relative error, subinterval limit,
      // workspace, result, estimated error
      gsl_integration_qags (&F, a, b, 0, 1e-10, 1000,
                            w, &result, &error);
    }

    // Free up memory
    gsl_integration_workspace_free(w);
  }

  // Check for illegal values
  if (result < 0.0) result = 0.0;
  if (result > 1.0) result = 0.0;

  return(result);
}

// Lookup - 06

inline double qexwaldrace_scl( std::vector<double> prm ) {
  /*
  Purpose:
  A scalar function that approximates the quantile function
  for the ex-Wald race model via linear interpolation.
  Arguments:
  prm[0] = A probability
  prm[1] = Choice (0 or 1)
  prm[2] = Threshold (1)
  prm[3] = Drift rate (1)
  prm[4] = Residual latency (1)
  prm[5] = Coefficient of drift (1)
  prm[6] = Threshold (0)
  prm[7] = Drift rate (0)
  prm[8] = Residual latency (0)
  prm[9] = Coefficient of drift (0)
  prm[10] = Indicator for whether residual latency
            impacts decision
  prm[11] = Indicator for whether to use numerical
            integration
  prm[12] = The upper boundary of the quantiles to
            explore
  prm[13] = The number of iterations to attempt
            during the search
  prm[14] = The precision of the estimate of the
            inverse CDF
  Returns:
  An estimated quantile.
  */

  // Initialize output
  double out = NA_REAL;

  // Create index
  std::vector<int> index = create_range( 1, 10 );
  for (int i = 6; i < 10; i++ ) index[i] = i - 3;
  index[0] = 7; // Set to probability

  // Check for valid inputs
  if ( ewr_param_verify( prm, index ) ) {

    // Rescale probabilities for joint distribution
    double p = prm[0];
    prm[0] = R_PosInf;
    double denom = pexwaldrace_scl( prm );
    prm[0] = p * denom;
    out = prm[ 0 ];

    cdf cur_cdf = pexwaldrace_scl;

    // Extract values governing linear interpolation
    double mn = 0.0;
    double mx = prm[12];
    int em_stop = prm[13];
    double err = prm[14];

    // std::vector<double> prm_int(11);
    // for ( int j = 0; j < 11; j++ ) prm_int[j] = prm[j];
    out = lin_interp_quantile( prm, cur_cdf,
                               mn, mx, em_stop, err,
                               0.0, R_PosInf );
  }

  return( out );
}

#endif // End include guard
