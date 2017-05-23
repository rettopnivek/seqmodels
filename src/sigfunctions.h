// Include guard to protect against multiple definitions error
#ifndef __SIGFUNCTIONS__
#define __SIGFUNCTIONS__

#include <Rcpp.h> // Includes certain libraries of functions
#include "levyfunctions.h" // For drift of 0
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Scalar functions for the random generation, density, distribution,
quantile, and moments functions of the shifted inverse Gaussian
(Wald) distribution.

References:
Dagpunar, J. (1988). Principles of Random Variate Generation.
  Oxford: Clarendon Press.
Heathcote, A. (2004a). Fitting Wald and ex-Wald distributions to
  response time data: An example using functions for the S-PLUS
  package. Behavior Research Methods Instruments & Computers, 36,
  678 - 694.
Heathcote, A. (2004b). rtfit.ssc. Retrieved May 5, 2017 from
  Psychonomic Society Web Archive: http://www.psychonomic.org/ARCHIVE/.

Index
Lookup - 01:  sig_param_verify
Lookup - 02:  rinvgauss_scl
Lookup - 03:  dinvgauss_scl
Lookup - 04:  pinvgauss_scl
Lookup - 05:  qinvgauss_scl
Lookup - 06:  minvgauss_scl
Lookup - 07:  sig_to_ig
*/

// Lookup - 01

inline bool sig_param_verify( std::vector<double> prm,
                              std::vector<int> index ) {
  /*
  Purpose:
  Function to verify whether inputs are admissable for the
  shifted inverse Gaussian distribution.
  Arguments:
  prm   - A vector of parameters
  index - The type of parameter being inputted, where...
  0 = skip
  1 = time
  2 = kappa (Threshold)
  3 = xi (Drift rate)
  4 = tau (Shift)
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
  // Variable for whether a time is present
  double t_yes = 0.0;

  for ( int i = 0; i < sz; i++ ) {

    if ( index[i] > 0 ) n_pass += 1;

    // Check whether inputted time is a real number and
    // greater than 0
    if ( index[i] == 1 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) {
        chk += 1;
        t_yes += prm[i];
      }
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

    // Check whether shift parameter is a real number and whether it
    // is greater than or equal to 0 and, if applicable, whether it
    // is less than the inputted response time
    if ( index[i] == 4 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] >= 0.0 ) ) {
        // Check whether shift parameter is less than the
        // given time
        if ( t_yes > 0.0 ) {
          if ( t_yes > prm[i] ) chk += 1;
        } else {
          chk += 1;
        }
      }
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

inline double rinvgauss_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function to generate random draws from the inverse
  Gaussian distribution.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = A threshold
  prm[1] = A drift rate
  prm[2] = A shift parameter
  prm[3] = The coefficient of drift
  prm[4] = A chi-square random deviate with 1 degree of freedom
  prm[5] = A unit uniform deviate
  Returns:
  A random deviate.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 2, 5 );

  // Check for valid inputs
  if ( sig_param_verify( prm, index ) ) {

    // Extract parameters
    double kappa = prm[0];
    double xi = prm[1];
    double tau = prm[2];
    double sigma = prm[3];

    // Rescale base on sigma
    xi = xi / sigma; kappa = kappa / sigma; sigma = 1.0;

    // Generate values from standard inverse Gaussian

    // Special case when drift rate is 0
    // Finishing times follow a Levy distribution
    if ( xi == 0.0 ) {
      out = levy_to_ig( prm[5], kappa, sigma, 3 );
    } else {

      // Extract random deviates
      double y2 = prm[4]; double u = prm[5];

      // Compute random deviates using script adapted from Heathcote
      // (2004b) which was in turn adapted from pp. 79 - 80 of
      // Dagpuner (1988)
      double y2onm = y2/xi;
      double r1 = ( 2.0 * kappa + y2onm -
                    sqrt( y2onm * ( 4.0 * kappa + y2onm ) ) ) /
                    ( 2.0 * xi );
      double r2 = pow( kappa/xi, 2.0 )/r1;

      if ( u < kappa/( kappa + xi*r1 ) ) {
        out = r1;
      }
      else {
        out = r2;
      }

    }

    // Apply shift parameter
    out = out + tau;
  }

  return(out);
}

// Lookup - 03

inline double dinvgauss_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar version for the density of the shifted inverse
  Gaussian distribution. Code is adapted from Smyth et al.
  (2017).
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = An observed time
  prm[1] = A threshold
  prm[2] = A drift rate
  prm[3] = A shift parameter
  prm[4] = The coefficient of drift
  prm[5] = An indicator for the log-likelihood
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
  if ( sig_param_verify( prm, index ) ) {

    // Extract parameters
    double rt = prm[0];
    double kappa = prm[1];
    double xi = prm[2];
    double tau = prm[3];
    double sigma = prm[4];

    // Rescale base on sigma
    xi = xi / sigma; kappa = kappa / sigma; sigma = 1.0;

    // Shift observed time
    double t = rt - tau;

    // Special case when drift rate is 0
    // Finishing times follow a Levy distribution
    if ( xi == 0.0 ) {
      out = levy_to_ig( t, kappa, sigma, 1 );
    } else {

      // Calculate log of density (Heathcote, 2004b)
      out = log( kappa ) - pow( kappa - xi*t, 2.0)/( 2.0 * t );
      out -= log( sqrt( 2.0 * M_PI * pow( t, 3.0 ) ) );
    }
  }

  // If specified return the density
  if (ln==0.0) out = exp(out);

  return( out );
}

// Lookup - 04

inline double pinvgauss_scl( std::vector<double> prm ) {
  /*
   Purpose:
   Scalar version for the distribution function of the
   shifted inverse Gaussian distribution.
   Arguments:
   prm - A vector of parameters, where...
   prm[0] = An observed time
   prm[1] = A threshold
   prm[2] = A drift rate
   prm[3] = A shift parameter
   prm[4] = The coefficient of drift
   Returns:
   The cumulative probability.
   */

  // Initialize output
  double out = 0.0;

  // Define index
  std::vector<int> index = create_range( 1, 5 );

  // Check for valid inputs
  if ( sig_param_verify( prm, index ) ) {

    // Extract parameters
    double rt = prm[0];
    double kappa = prm[1];
    double xi = prm[2];
    double tau = prm[3];
    double sigma = prm[4];

    // Rescale base on sigma
    xi = xi / sigma; kappa = kappa / sigma; sigma = 1.0;

    // For positive infinity
    if ( rt == R_PosInf ) {
      out = 1.0;
    } else {

      // Shift observed time
      double t = rt - tau;

      // Special case when drift rate is 0
      // Finishing times follow a Levy distribution
      if ( xi == 0.0 ) {
        out = levy_to_ig( t, kappa, sigma, 2 );
      } else {

        // Compute distribution function (Heathcote, 2004b)
        double sqrtt = sqrt( t );
        double term1 = ( xi*t - kappa )/sqrtt; // k1
        double term2 = ( xi*t + kappa)/sqrtt; // k2
        double term3 = exp( 2.0 * kappa * xi );
        double Phi1 = R::pnorm( term1, 0.0, 1.0, 1, 0 );
        double Phi2 = R::pnorm( -term2, 0.0, 1.0, 1, 0 );

        // Protection against numerical error
        if ( ( std::isinf( term3 ) ) || ( Phi2 == 0.0 ) ) {

          out = -pow( term1, 2.0 )/2.0 - 0.94/pow( term2, 2.0 );
          out = exp( out ) / ( term2*pow( 2.0 * M_PI, 0.5 ) );
          out += Phi1;

        } else {
          out = Phi1 + term3 * Phi2;
        }

      }
    }


  }

  return( out );
}

// Lookup - 05

inline double qinvgauss_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for quantile function of the shifted
  inverse Gaussian distribution.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = p (Probability)
  prm[1] = A threshold
  prm[2] = A drift rate
  prm[3] = A shift parameter
  prm[4] = The coefficient of drift
  prm[5] = mx (Uppermost quantile to examine)
  prm[6] = em_stop (Max number of possible iterations)
  prm[7] = err (Desired level of accuracy in estimate)
  Returns:
  The quantile associated with the inputted cumulative probability.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 1, 5 );
  index[0] = 6; // For probability

  // Check for inadmissable values
  if ( sig_param_verify( prm, index ) ) {

    cdf cur_cdf = pinvgauss_scl;

    // Extract values governing linear interpolation
    double mx = prm[5];
    int em_stop = prm[6];
    double err = prm[7];

    out = lin_interp_quantile( prm, cur_cdf,
                               0.0, mx, em_stop, err,
                               0.0, R_PosInf );

  }

  return( out );
}

// Lookup - 06

inline std::vector<double> minvgauss_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function to compute the moments (mean, variance,
  standard deviation, skewness, and exces kurtosis) for the
  shifted inverse Gaussian.
  Arguments:
  prm - A vector of parameters, where...
  prm[1] = A threshold
  prm[2] = A drift rate
  prm[3] = A shift parameter
  prm[4] = The coefficient of drift
  Returns:
  A vector with the mean, variance, standard deviation,
  skewness, and excess kurtosis.
  */

  // Initialize output
  std::vector<double> out(5);
  for ( int i = 0; i < 5; i++ ) out[i] = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 2, 5 );

  // Check for inadmissable values
  if ( sig_param_verify( prm, index ) ) {

    // Extract parameters
    double kappa = prm[0];
    double xi = prm[1];
    double tau = prm[2];
    double sigma = prm[3];

    // Rescale base on sigma
    xi = xi / sigma; kappa = kappa / sigma; sigma = 1.0;

    // Convert from Brownian transformation
    double mu = kappa/xi;
    double lambda = pow( kappa, 2.0 )/pow( sigma, 2.0 );

    // Mean
    out[0] = mu + tau;

    // Variance
    out[1] = pow( mu, 3.0 )/lambda;

    // Standard deviation
    out[2] = pow( out[1], 0.5 );

    // Skewness
    out[3] = 3.0 * pow( mu/lambda, 0.5 );

    // Exces kurtosis
    out[4] = 15.0 * mu / lambda;
  }

  return( out );
}

// Lookup - 07

inline double sig_to_ig( double t, double kappa,
                         double xi, double sigma,
                         double u, int ver ) {
  /*
  Purpose:
  A convenience function that computes the density,
  distribution, or random generation functions given
  subset of parameters for the inverse Gaussian.
  Arguments:
  t     - An observed time or a chi-square deviate
  kappa - A threshold
  xi    - A drift rate
  sigma - The coefficient of drift
  u     - A unit uniform deviate
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

    std::vector<double> ig_prm(6);
    ig_prm[0] = t;
    ig_prm[1] = kappa;
    ig_prm[2] = xi;
    ig_prm[3] = 0.0;
    ig_prm[4] = sigma;
    ig_prm[5] = 0.0;

    out = dinvgauss_scl( ig_prm );
  }

  // Distribution function
  if ( ver == 2 ) {

    std::vector<double> ig_prm(5);
    ig_prm[0] = t;
    ig_prm[1] = kappa;
    ig_prm[2] = xi;
    ig_prm[3] = 0.0;
    ig_prm[4] = sigma;

    out = pinvgauss_scl( ig_prm );
  }

  // Random generation
  if ( ver == 3 ) {

    std::vector<double> ig_prm(6);
    ig_prm[0] = kappa;
    ig_prm[1] = xi;
    ig_prm[2] = 0.0;
    ig_prm[3] = sigma;
    ig_prm[4] = t;
    ig_prm[5] = u;
    double ig_rd = rinvgauss_scl( ig_prm );
  }

  return( out );
}

#endif // End include guard
