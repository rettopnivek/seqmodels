// Include guard to protect against multiple definitions error
#ifndef __LEVYFUNCTIONS__
#define __LEVYFUNCTIONS__

#include <Rcpp.h> // Includes certain libraries of functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Scalar functions for the density, distribution, quantile,
and random generation functions of the Levy distribution.

References:
Applebaum, D. (2010). Lectures on Levy processes and stochastic
  calculus, Braunschweig; Lecture 2: Levy processes. Retrieved from
  http://www.applebaum.staff.shef.ac.uk/Brauns2notes.pdf.
Siegrist, K. (1997). The Levy distribution. Retrieved from
  http://www.math.uah.edu/stat/special/Levy.html

Index
Lookup - 01:  levy_param_verify
Lookup - 02:  dlevy_scl
Lookup - 03:  plevy_scl
Lookup - 04:  qlevy_scl
Lookup - 05:  rlevy_scl
Lookup - 06:  levy_to_ig
*/

// Lookup - 01

inline bool levy_param_verify( std::vector<double> prm,
                               std::vector<int> index ) {
  /*
  Purpose:
  Function to verify whether inputs are admissable for the
  Levy distribution.
  Arguments:
  prm   - A vector of parameters
  index - The type of parameter being inputted, where...
  0 = skip
  1 = datum
  2 = mu
  3 = sigma
  4 = probability
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
  // Variable to store inputted datum
  std::vector<double> dat(2);
  dat[0] = 0.0; dat[1] = NA_REAL;

  for ( int i = 0; i < sz; i++ ) {

    if ( index[i] > 0 ) n_pass += 1;

    if ( index[i] == 1 ) {
      // Store datum to be compared later against the
      // parameter mu
      dat[0] = 1.0; dat[1] = prm[i];
    }

    // Check whether location parameter is a real number
    // and, if needed, whether inputted datum is larger
    if ( index[i] == 2 ) {
      if ( !Rcpp::NumericVector::is_na( prm[i] ) ) {
        chk += 1;
        if ( dat[0] == 1.0 ) {
          if ( ( !Rcpp::NumericVector::is_na( dat[1] ) ) &&
               ( dat[1] > prm[i] ) ) chk += 1;
        }
      }
    }

    // Check whether scale parameter is a real number and whether
    // it is greater than 0
    if ( index[i] == 3 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether probability is a real number and whether it
    // is appropriately bounded
    if ( index[i] == 4 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] >= 0.0 ) && ( prm[i] <= 1.0 ) ) chk += 1;
    }

  }

  // Test how many checks were passed
  out = chk == n_pass;

  return( out );
}

// Lookup - 02

inline double dlevy_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for the density for the Levy distribution.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = x
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  prm[3] = ln (Indicator for log-density)
  Returns:
  The density (or log-density).
  */

  // Initialize output
  double out = log( 0.0 );

  // Extract indicator for log-likelihood
  int ln = prm[3];

  // Define index
  std::vector<int> index = create_range( 1, 3 );

  // Check for valid inputs
  if ( levy_param_verify( prm, index ) ) {

    double t = prm[0];
    double mu = prm[1];
    double sigma = prm[2];

    double p1 = pow( sigma / ( 2.0 * M_PI ), 0.5 );
    double p2 = pow( t - mu, 3.0/2.0 );
    double p3 = exp( -sigma / ( 2.0 * ( t - mu ) ) );

    out = log( p1 ) + log( p3 ) - log( p2 );
  }

  if ( ln == 0 ) out = exp( out );

  return( out );
}

// Lookup - 03

inline double plevy_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for the CDF of the Levy distribution.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = x
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  Returns:
  The CDF.
  */

  // Initialize output
  double out = 0.0;

  // Define index
  std::vector<int> index = create_range( 1, 3 );

  // Check for valid inputs
  if ( levy_param_verify( prm, index ) ) {

    double t = prm[0];
    double mu = prm[1];
    double sigma = prm[2];

    double x = pow( sigma / ( 2.0 * ( t - mu ) ), 0.5 );

    out = erfc( x );
  }

  // Check for inadmissable values
  if ( out > 1.0 ) out = 1.0;
  if ( out < 0.0 ) out = 0.0;

  return( out );
}

// Lookup - 04

inline double qlevy_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for computing quantiles given a cumulative
  probability for the Levy distribution
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = A probability
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  Returns:
  The associated quantile with the cumulative probability.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 1, 3 );
  index[0] = 4; // First value is a probability

  // Check for valid inputs
  if ( levy_param_verify( prm, index ) ) {

    double p = prm[0];
    double mu = prm[1];
    double sigma = prm[2];

    if ( p == 1.0 ) {
      out = R_PosInf;
    } else if ( p == 0.0 ) {
      out = 0.0;
    } else {
      double denom = pow( R::qnorm( 1.0 - p/2, 0.0, 1.0, 1, 0 ), 2.0 );
      out = sigma / denom + mu;
    }
  }

  return( out );
}


// Lookup - 05

inline double rlevy_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for generating random deviates from a
  Levy distribution
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = A probability
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  Returns:
  A random deviate.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 1, 3 );
  index[0] = 4; // First value is a probability

  // Check for valid inputs
  if ( levy_param_verify( prm, index ) ) {

    double u = prm[0];
    double mu = prm[1];
    double sigma = prm[2];

    double denom = pow( R::qnorm( 1.0 - u/2, 0.0, 1.0, 1, 0 ), 2.0 );

    out = sigma / denom + mu;
  }

  return( out );
}

// Lookup - 06

inline double levy_to_ig( double t,
                          double kappa,
                          double sigma, int ver ) {
  /*
  Purpose:
  A convenience function that computes the density,
  distribution, or random generation functions for the
  Levy distribution given parameters from the inverse
  Gaussian distribution.
  Arguments:
  t     - An observed time or unit uniform deviate
  kappa - A threshold
  sigma - The coefficient of drift
  ver   - If 1, computes the density; if 2, computes the
          distribution function; if 3, generates a
          random deviate.
  Returns:
  A continuous value.
  */

  // Initialize output
  double out = NA_REAL;

  // Density function
  if ( ver == 1 ) {
    std::vector<double> prm_levy(4);
    prm_levy[0] = t;
    prm_levy[1] = 0.0;
    prm_levy[2] = pow( kappa, 2.0 )/pow( sigma, 2.0 );
    prm_levy[3] = 1;
    out = dlevy_scl( prm_levy );
  }

  // Distribution function
  if ( ver == 2 ) {
    std::vector<double> prm_levy(3);
    prm_levy[0] = t;
    prm_levy[1] = 0.0;
    prm_levy[2] = pow( kappa, 2.0 )/pow( sigma, 2.0 );
    out = plevy_scl( prm_levy );
  }

  // Random generation
  if ( ver == 3 ) {
    std::vector<double> prm_levy(3);
    prm_levy[0] = t;
    prm_levy[1] = 0.0;
    prm_levy[2] = pow( kappa, 2.0 )/pow( sigma, 2.0 );
    out = rlevy_scl( prm_levy );
  }

  return( out );
}

#endif // End include guard
