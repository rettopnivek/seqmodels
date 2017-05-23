// Include guard to protect against multiple definitions error
#ifndef __MISCFUNCTIONS__
#define __MISCFUNCTIONS__

#include <Rcpp.h>

/*
Purpose:
Miscellaneous Rcpp code for use in other functions.

Index:
Lookup - 01:  create_index
Lookup - 02:  create_seq
Lookup - 03:  create_range
Lookup - 04:  erfc
Lookup - 05:  minMax
Lookup - 06:  linInterp
Lookup - 07:  cdf
Lookup - 08:  lin_interp_quantile
*/

// Lookup - 01

inline Rcpp::IntegerVector create_index( int N, int mx ) {
  /*
  Purpose:
  Creates an integer vector of size 'N' that can be used to index
  another vector of size 'mx'. When 'mx' is less than 'N', indices
  are repeated in sequential order. If 'mx' is larger than 'N',
  indices are truncated early.
  Arguments:
  N  - The size of the new index vector
  mx - The size of the original vector
  Returns:
  A integer vector of size 'N' that can be used to index a vector
  of size 'mx'.
  */

  Rcpp::IntegerVector out( N );

  int inc = 0;
  for ( int n = 0; n < N; n++ ) {
    out[ n ] = inc;
    inc += 1;
    if ( inc == mx ) inc = 0;
  }

  return( out );
}

// Lookup - 02

inline std::vector<double> create_seq( double sta, double sto,
                                       double cnt ) {
  /*
  Purpose:
  Generates a numeric vector with a specified number of values
  going from a starting value to an end point in regular increments.
  Arguments:
  sta - The starting value
  sto - The end point
  cnt - The number of values in the sequence
  Returns:
  A vector with the sequence of values.
  */

  double inc = ( sto - sta )/( cnt - 1.0 );
  std::vector<double> out( cnt );
  out[0] = sta;

  for ( int i = 1; i < cnt; i++ ) {
    out[i] = out[i-1] + inc;
  }

  return( out );
}

// Lookup - 03

inline std::vector<int> create_range( int sta, int sto ) {
  /*
  Purpose:
  Creates an integer matrix that goes from a starting value
  to an end point in unit increments.
  Arguments:
  sta - The integer to start at
  sto - The integer to end at
  Returns:
  An integer vector starting from 'sta' and ending at 'sto' in
  unit increments.
  */

  // Determine size of vector
  int sz = sto - sta + 1;

  // Initialize output
  std::vector<int> out( sz );

  // Fill vector with values
  out[0] = sta;
  for ( int i = 1; i < sz; i++ ) out[i] = out[i-1] + 1;

  return( out );
}

// Lookup - 04

inline double erfc( double x ) {
  /*
  Purpose:
  Computes the complementary error function.
  Arguments:
  x - A continuous value
  Returns:
  One minus the error function.
  */

  double erf = 2.0 * R::pnorm( x * pow( 2.0, 0.5 ),
                               0.0, 1.0, 1, 0 ) - 1.0;

  return( 1.0 - erf );
}

// Lookup - 05

inline int minMax( double comp, std::vector<double> x,
                   int Direction ) {
  /*
  Purpose:
  A function to determine the smallest value of 'x' larger than
  'comp' or the largest value of 'x' smaller than 'comp'.
  Arguments:
  comp      - A numeric value to compare against the vector 'x'
  x         - A numeric vector of values sorted from smallest to
              largest
  Direction - If 1, the function searches for the smallest value of
              'x' that is greater than 'comp'; If 0, the function
              searches for the largest value of 'x' that is
              less than 'comp'
  Returns:
  An index (from 0 to the size of the vector 'x' minus 1) of the
  smallest value of 'x' larger than 'comp' or the largest value
  of 'x' smaller than 'comp'.
  */

  // Initialize output
  int out = 0;

  // Extract the size of the vector
  int sz = x.size();

  // Loop over the elements in x
  for ( int i = 0; i < sz; i++ ) {
    // Identify values of x larger than comp
    if ( (Direction == 1) && ( x[i] <= comp ) ) out += 1;
    // Identify values of x smaller than comp
    if ( (Direction == 0) && ( x[i] < comp ) ) out += 1;
  }

  // Adjust index to start at 0 and end at sz - 1
  if ( Direction == 1 ) {
    // If no value of x is larger than comp, take the
    // last value
    if ( out >= sz ) out = sz - 1;
  }

  // Adjust index to start at 0 and end at sz - 1
  if ( Direction == 0 ) {
    out -= 1;
    // If every value of x is larger than comp, take the
    // first value
    if ( out < 0 ) out = 0;
  }

  return( out );
}

// Lookup - 06

inline double linInterp( double yN, double y0, double y1,
                         double x0, double x1 ) {
  /*
  Purpose:
  Given a set of starting and end x and y values for a line,
  computes a new value of x given a value of y via linear
  interpolation.
  Arguments:
  yN - A new value on the y-axis lying somewhere on the line
  y0 - The start point for the line on the y-axis
  y1 - The end point for the line on the y-axis
  x0 - The start point for the line on the x-axis
  x1 - The end point for the line on the x-axis
  Returns:
  The associated value on the line for the x-axis based on
  the y-axis value.
  Notes:
  If the value of yN is outside the boundaries of y0 or y1, then
  the returned value is an extrapolation instead.
  */

  double b1 = ( y1 - y0 ) / ( x1 - x0 ); // Slope
  double b0 = y1 - b1*x1; // Intercept
  double num = yN - b0;
  double xN = ( num )/b1;

  return( xN );
}

// Lookup - 07

// Declare a new type for a CDF
typedef double (*cdf)(std::vector<double>);

// Lookup - 08

inline double lin_interp_quantile( std::vector<double> prm,
                                   cdf my_cdf,
                                   double mn,
                                   double mx,
                                   int em_stop,
                                   double err,
                                   double low_b,
                                   double up_b ) {
  /*
  Purpose:
  A function to compute the quantile associated with a given
  cumulative probability via linear interpolation.
  Arguments:
  p       - A cumulative probability
  my_cdf  - A CDF that takes a standard vector of inputs
            (first value must be a cumulative probability)
  prm     - The standard vector of inputs for my_cdf
  mn      - The lowest quantile value to start exploring
  mx      - The highest quantile value to start exploring
  em_stop - The maximum number of iterations for which an
            approximation is attempted
  err     - The desired accuracy to which the cumulative
            probability is approximated
  low_b   - The lower boundary of support
  up_b    - The upper boundary of support
  Returns:
  The estimated quantile associated with the inputted cumulative
  probability.
  */

  // Initialize output
  double out = NA_REAL;

  // Extract cumulative probability for quantile to estimate
  double p = prm[0];

  if ( p == 1.0 ) {
    out = up_b;
  } else if ( p == 0.0 ) {
    out = low_b;
  } else if ( ( p > 0.0 ) && ( p < 1.0 ) ) {

    // Define counter for iterations
    int em;

    // Define vector of cumulative probabilities
    std::vector<double> prb( 5 );
    prb[0] = 0.0; // Initialize values
    prb[4] = 0.0;

    // Check whether lower limit of the quantiles to
    // explore is too high
    if ( low_b < 0.0 ) {

      em = 0;
      prb[0] = 1.0;
      while( (prb[0] > p) && ( em < em_stop ) ) {

        // If maximum quantile is too large, decrease by
        // unit decrements
        if ( em >= 1 ) mn -= 1.0;

        // Compute cumulative probability
        prm[0] = mn;
        prb[0] = my_cdf( prm );
        em += 1;
      }

      // If an lower limit cannot be found return NA
      if ( em >= em_stop ) goto output;
    }

    // Check whether upper limit of the quantiles to
    // explore is too low
    em = 0;
    while( (prb[4] < p) && ( em < em_stop ) ) {

      // If maximum quantile is too small, increase by
      // .5 s increments
      if ( em >= 1 ) mx += .5;

      // Compute cumulative probability
      prm[0] = mx;
      prb[4] = my_cdf( prm );
      em += 1;
    }

    // If an upper limit cannot be found return NA
    if ( em >= em_stop ) goto output;

    // Otherwise, create a sequence of quantiles
    std::vector<double> q_val = create_seq( mn, mx, 5 );

    // Define variable for current quantile value
    double q;

    // Compute cumulative probability for each point
    for ( int i = 1; i < 4; i++ ) {
      prm[0] = q_val[i];
      prb[i] = my_cdf( prm );
    }

    // Lower boundary is based on last value with a cumulative
    // probability smaller than the desired value
    int lb = minMax( p, prb, 0 );

    // Upper boundary is based on the first value with
    // a cumulative probability greater than the desired value
    int ub = minMax( p, prb, 1 );

    // Initialize variables
    double epsilon = 1; // Error in approximation
    em = 1; // Counter for iterations
    // Lower and upper limits for response times
    double q0 = q_val[ lb ]; double q1 = q_val[ ub ];
    // Lower and upper limits for cumulative probability
    double p0 = prb[ lb ]; double p1 = prb[ ub ];

    // Loop until error in approximation is sufficiently
    // small
    while( ( std::abs( epsilon ) > err ) &&
           ( em <= em_stop ) ) {

      // Compute new q via linear interpolation
      double qN = linInterp( p, p0, p1, q0, q1 );

      // Associated cumulative probability
      prm[0] = qN;
      double pN = my_cdf( prm );

      // Take difference between estimate and desired
      epsilon = pN - p;

      // Adjust bounds of linear interpolation based on
      // whether estimate was smaller than desired
      if ( pN < p ) {

        // Set lower bounds to new estimates
        p0 = pN;
        q0 = qN;

        // If possible, scale upper bounds by 50%
        double qP = qN + 0.5 * (q1 - qN);
        prm[0] = qP;
        double pP = my_cdf( prm );
        if ( pP > p ) {
          p1 = pP;
          q1 = qP;
        }

      }

      // Adjust bounds of linear interpolation based on
      // whether estimate was larger than desired
      if ( pN > p ) {

        // Set upper bounds to new estimates
        p1 = pN;
        q1 = qN;

        // If possible, scale lower bounds by 50%
        double qP = q0 + 0.5 * (qN - q0);
        prm[0] = qP;
        double pP = my_cdf( prm );
        if ( pP < p ) {
          p0 = pP;
          q0 = qP;
        }
      }

      q = qN;
      em += 1;
    }

    // Save final output
    out = q;
  }

  // Label for goto statment
  output:
  return( out );
}

#endif // End include guard
