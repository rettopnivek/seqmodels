// Include guard to protect against multiple definitions error
#ifndef __EMGFUNCTIONS__
#define __EMGFUNCTIONS__

#include <Rcpp.h> // Includes certain libraries of functions
#include "miscfunctions.h" // Linear interpolation

/*
Purpose:
Scalar functions for the calculation of the density, distribution,
quantile, and random number generation functions of the exponentially
modified gaussian distribution.

References:
Luce, R. D. (1986). Response times: Their role in inferring
  elementary mental organization. New York, New York: Oxford
  University Press.
Olivier J., & Norberg M. M. (2010). Positively skewed data: Revisiting
  the Box−Cox power transformation. International Journal of
  Psychological Research, 3(1), 68 − 75.
Ulrich, R., & Miller, J. (1994). Effects of outlier exclusion on
  reaction time analysis. Journal of Experimental Psychology: General,
  123, 34 – 80. doi:10.1037/0096-3445.123.1.34.

Index
Lookup - 01:  emg_param_verify
Lookup - 02:  remg_scl
Lookup - 03:  demg_scl
Lookup - 04:  pemg_scl
Lookup - 05:  qemg_scl
Lookup - 06:  memg_scl
*/

// Lookup - 01

inline bool emg_param_verify( std::vector<double> prm,
                              std::vector<int> index ) {
  /*
  Purpose:
  Function to verify whether inputs are admissable for the
  exponentially modified gaussian distribution.
  Arguments:
  prm   - A vector of parameters
  index - The type of parameter being inputted, where...
  0 = skip
  1 = response time
  2 = location parameter
  3 = scale parameter
  4 = rate parameter
  5 = probability
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
    if ( index[i] == 1 ) {
      if ( !Rcpp::NumericVector::is_na( prm[i] ) ) chk += 1;
    }

    // Check whether location parameter is a real number
    if ( index[i] == 2 ) {
      if ( !Rcpp::NumericVector::is_na( prm[i] ) ) chk += 1;
    }

    // Check whether scale parameter is a real number and whether it
    // is positive
    if ( index[i] == 3 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether rate parameter is a real number and whether it
    // is positive
    if ( index[i] == 4 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether probability is a real number and whether it
    // is appropriately bounded
    if ( index[i] == 5 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( (prm[i] >= 0.0) | (prm[i] <= 1.0) ) ) chk += 1;
    }

  }

  // Test how many checks were passed
  out = chk == n_pass;

  return( out );
}

// Lookup - 02

inline double remg_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function to generate random draws from convolution of
  exponential and gaussian random variables.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = mu (Location)
  prm[1] = sigma (Scale)
  prm[2] = lambda (Rate)
  Returns:
  A random deviate.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 2, 4 );

  // Check for valid inputs
  if ( emg_param_verify( prm, index ) ) {

    double mu = prm[0];
    double sigma = prm[1];
    double lambda = prm[2];

    // Convert lambda to rate
    lambda = 1.0/lambda;

    // Sum of independent draws for a gaussian and an
    // exponential variable
    out = R::rnorm( mu, sigma ) + R::rexp( lambda );
  }

  return( out );
}

// Lookup - 03

inline double demg_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for density of the convolution of an exponential
  and gaussian random variable.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = x (response time)
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  prm[3] = lambda (Rate)
  prm[4] = ln (Indicator for log-density)
  Returns:
  The density (or log-density).
  */

  // Initialize output
  double out = 0.0;

  // Determine if log-density should be returned
  int ln = prm[4];

  // Define index
  std::vector<int> index = create_range( 1, 4 );

  // Check for valid inputs
  if ( emg_param_verify( prm, index ) ) {

    // Extract parameters
    double x = prm[0];
    double mu = prm[1];
    double sigma = prm[2];
    double lambda = prm[3];

    double sigma2 = pow( sigma, 2.0 );
    double p1 = lambda/2.0;
    double p2 = 2.0*mu + lambda*sigma2 - 2.0*x;
    double p3 = ( mu + lambda*sigma2 - x )/(pow(2.0,.5)*sigma);

    out = p1*exp(p1*p2)*erfc(p3);
  }

  if ( out < 0.0 ) out = 0.0;
  if ( ln == 1 ) out = log( out );

  return( out );
}

// Lookup - 04

inline double pemg_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for CDF of the convolution of an exponential
  and gaussian random variable.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = x (response time)
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  prm[3] = lambda (Rate)
  Returns:
  The cumulative probability.
  */

  // Initialize output
  double out = 0.0;

  // Define index
  std::vector<int> index = create_range( 1, 4 );

  // Check for valid inputs
  if ( emg_param_verify( prm, index ) ) {

    // Extract parameters
    double x = prm[0];
    double mu = prm[1];
    double sigma = prm[2];
    double lambda = prm[3];

    // Flip lambda
    double u = lambda*(x - mu);
    double v = lambda*sigma;

    double p1 = R::pnorm( u, 0.0, v, 1, 0 );
    double p2 = -u + pow(v,2.0)/2.0;
    double p3 = R::pnorm( u, pow(v,2.0), v, 1, 1 );

    out = p1 - exp( p2 + p3);
  }

  if ( (out < 0.0) || (out > 1.0) ) out = 0.0;

  return( out );
}

// Lookup - 05

inline double qemg_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function for quantile function of the convolution of an
  exponential and gaussian random variable.
  Arguments:
  prm - A vector of parameters, where...
  prm[0] = p (Probability)
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  prm[3] = lambda (Rate)
  prm[4] = mn (Lowermost quantile to examine)
  prm[5] = mx (Uppermost quantile to examine)
  prm[6] = em_stop (Max number of possible iterations)
  prm[7] = err (Desired level of accuracy in estimate)
  Returns:
  The quantile associated with the inputted cumulative probability.
  */

  // Initialize output
  double out = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 1, 4 );
  index[0] = 5; // For probability

  // Check for inadmissable values
  if ( emg_param_verify( prm, index ) ) {

    cdf cur_cdf = pemg_scl;

    // Extract values governing linear interpolation
    double mn = prm[4];
    double mx = prm[5];
    int em_stop = prm[6];
    double err = prm[7];

    out = lin_interp_quantile( prm, cur_cdf,
                               mn, mx, em_stop, err,
                               R_NegInf, R_PosInf );

  }

  return( out );
}

// Lookup - 06

inline std::vector<double> memg_scl( std::vector<double> prm ) {
  /*
  Purpose:
  Scalar function to compute the moments (mean, variance,
  standard deviation, skewness, and kurtosis) for the
  exponentially modified gaussian distribution.
  Arguments:
  prm - A vector of parameters, where...
  prm[1] = mu (Location)
  prm[2] = sigma (Scale)
  prm[3] = lambda (Rate)
  Returns:
  A vector with the mean, variance, standard deviation,
  skewness, and kurtosis.
  */

  // Initialize output
  std::vector<double> out(5);
  for ( int i = 0; i < 5; i++ ) out[i] = NA_REAL;

  // Define index
  std::vector<int> index = create_range( 2, 4 );

  // Check for inadmissable values
  if ( emg_param_verify( prm, index ) ) {

    // Extract parameters
    double mu = prm[0];
    double sigma = prm[1];
    double lambda = prm[2];

    double inv_pow_2 = 1.0/( pow( sigma, 2.0 ) * pow( lambda, 2.0 ) );
    double inv_pow_3 = 1.0 / ( pow( sigma, 3.0 ) * pow( lambda, 3.0 ) );
    double inv_pow_4 = 1.0 / ( pow( sigma, 4.0 ) * pow( lambda, 4.0 ) );


    // Mean
    out[0] = mu + 1.0/lambda;

    // Variance
    out[1] = pow( sigma, 2.0 ) + 1.0/pow( lambda, 2.0 );

    // Standard deviation
    out[2] = pow( out[1], 0.5 );

    // Skewness
    out[3] = 2.0 * inv_pow_3 * pow( 1.0 + inv_pow_2, -3.0/2.0 );

    // Extreme kurtosis
    double num = 1.0 + 2.0 * inv_pow_2 + 3.0 * inv_pow_4;
    double denom = pow( 1.0 + inv_pow_2, 2.0 );
    out[4] = 3.0 * num / denom - 3.0;
  }

  return( out );
}


#endif // End include guard
