// Include guard to protect against multiple definitions error
#ifndef __BLURTONETAL__
#define __BLURTONETAL__

#include <Rcpp.h>
#include "wpparamverify.h"
#include "miscfunctions.h"

/*
Purpose:
Implementation of the adaptive algorithm of Blurton et al.
(2012) for computing the joint distribution function of the
two boundary wiener process.

References:
Blurton, S. P., Kesselmeier, M., & Gondan, M. (2012).
  Fast and accurate calculations for cumulative first-passage
  time distributions in Wiener diffusion models. Journal of
  Mathematical Psychology, 56, 470-475.

Index:
Lookup - 01:  exp_pnorm
Lookup - 02:  K_large
Lookup - 03:  K_small
Lookup - 04:  Pu
Lookup - 05:  Fl_lower
Lookup - 06:  Fs0_lower
Lookup - 07:  Fs_lower
Lookup - 08:  F_lower
*/

// Lookup - 01

inline double exp_pnorm( double a, double b ) {
  /*
  Purpose:
  Approximates (e^a)*Phi(b), where Phi is the standard normal
  cumulative distribution.
  Arguments:
  a - A continuous value
  b - A continuous value
  Returns:
  An approximation of (e^a)*Phi(b).
  */

  double r = exp(a)*R::pnorm(b,0.0,1.0,1,0);

  if ( Rcpp::NumericVector::is_na(r) && b < -5.5 ) {

    double p1 = 1.0/pow(2.0,0.5);
    double p2 = exp(a - b*b/2.0);
    double p3 = 0.5641882/b/b/b - 1.0/b/pow(M_PI,.5);
    r = p1*p2*p3;

  }

  return(r);
}

// Lookup - 02

inline int K_large( double t, double v, double a,
                    double w, double epsilon ) {
  /*
  Purpose:
  Computes the number of terms required for the large time
  representation.
  Arguments:
  t       - A decision time
  v       - A drift rate
  a       - Boundary separation
  w       - The proportion for the starting time
  epsilon - The degree of precision
  Returns:
  The number of terms needed for the large time representation.
  */

  double sqrtL1 = pow(1.0/t,.5)*(a/M_PI);
  double p1 = -2.0/t*a*a/M_PI/M_PI;
  double p2 = log( (epsilon*M_PI*t/2.0)*(v*v + M_PI*M_PI/a/a) );
  double p3 = v*a*w + v*v*t/2.0;
  double sqrtL2 = pow( std::max( 1.0, p1*p2 + p3 ), 0.5 );

  return( ceil( std::max( sqrtL1, sqrtL2 ) ) );
}

// Lookup - 03

inline int K_small( double t, double v, double a,
                    double w, double epsilon ) {
  /*
  Purpose:
  Computes the number of terms required for the small time
  representation.
  Arguments:
  t       - A decision time
  v       - A drift rate
  a       - Boundary separation
  w       - The proportion for the starting time
  epsilon - The degree of precision
  Returns:
  The number of terms needed for the small time representation.
  */

  int out;

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  if ( std::abs(v) < machine_double ) {
    double p1 = pow(t,.5)/2.0/a;
    double p2 = R::qnorm( std::max( 0.0,
                                    std::min( 1.0,
                                              epsilon/( 2.0-2.0*w ) ) ),
                                              0.0, 1.0, 1, 0 );
    out = ceil( std::max( 0.0, w/2.0 - p1*p2 ) );
  } else {
    // Positive drift
    if ( v > 0 ) {
      epsilon = exp(-2.0*a*w*v)*epsilon;
      v = -1.0*v;
    }

    std::vector<double> S(4);
    S[0] = 0.0;
    S[1] = w - 1.0 + 1.0/2.0/v/a *
      log(epsilon/2.0 * (1-exp(2.0*v*a)));
    S[2] = (0.535 * sqrt(2.0*t) + v*t + a*w)/2.0/a;
    double p1 = epsilon*a/0.3/pow(2.0*M_PI*t,
                                  0.5)*exp(v*v*t/2.0 + v*a*w);
    double p2 = R::qnorm( std::max( 0.0,
                                    std::min( 1.0, p1 ) ),
                                    0.0, 1.0, 1, 0 );
    S[3] = w/2.0 - pow(t,0.5)/2.0/a * p2;
    double mx = S[0];
    for (int i = 1; i < 4; i++) {
      mx = std::max( mx, S[i] );
    }
    out = ceil( mx );
  }

  return( out );
}

// Lookup - 04

inline double Pu( double v, double a, double w ) {
  /*
  Purpose:
  Computes the probability for absorption at the upper barrier.
  Arguments:
  v - A drift rate
  a - Boundary separation
  w - The proportion governing the starting point
  Returns:
  The probability that the accumulation reaches the upper
  boundary.
  */

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  // Define variable for output
  double out = 1.0 - w; // Drift is near zero or w is near 1
  double e = exp( -2.0*v*a*( 1.0-w ) );
  if (e == R_PosInf)
    out = 1.0;
  // Drift isn't near zero and w isn't near 1
  if ( std::abs( e - 1.0 ) >= machine_double )
    out = (1.0 - e) / (exp(2.0*v*a*w) - e); // Standard case

  return( out );
}


// Lookup - 05

inline double Fl_lower( double t, double v, double a,
                        double w, int K ) {
  /*
  Purpose:
  Computes the large time representation of the lower
  sub-distribution.
  Arguments:
  t - A decision time
  v - A drift rate
  a - Boundary separation
  w - The proportion governing the starting point
  K - The number of terms to compute in the infinite sum
  Returns:
  The large time representation.
  */

  // Define variable for output
  double out = 0.0;

  for (int k = K; k > 0; k--) {
    double oldOut = out;
    double p1 = k/(v*v + k*k*M_PI*M_PI/a/a);
    double p2 = exp(-v*a*w - 1.0/2.0*v*v*t -
                    1.0/2.0*k*k*M_PI*M_PI/a/a*t);
    double p3 = sin(M_PI*k*w);
    out = oldOut - p1*p2*p3;
  }
  out = Pu(v, a, w) + 2.0*M_PI/a/a * out;

  return( out );
}

// Lookup - 06

inline double Fs0_lower( double t, double a, double w, int K ) {
  /*
  Purpose:
  A version of Fs_lower for when the drift rate is zero.
  Arguments:
  t - A decision time
  a - Boundary separation
  w - The proportion governing the starting point
  K - The number of terms to compute in the infinite sum
  Returns:
  The small time representation.
  */

  double out = 0.0;

  for (int k = K; k >= 0; k--) {
    double p1 = (-2.0*k - 2.0 + w)*a/pow(t,0.5);
    double p2 = (-2.0*k - w)*a/pow(t,0.5);
    out = R::pnorm( p1, 0.0, 1.0, 1, 0 ) +
      R::pnorm( p2, 0.0, 1.0, 1, 0 );
  }

  return( 2.0*out );
}

// Lookup - 07

inline double Fs_lower( double t, double v, double a,
                        double w, int K ) {
  /*
  Purpose:
  Computes the small time representation of the upper
  sub-distribution.
  Arguments:
  t - A decision time
  v - A drift rate
  a - Boundary separation
  w - The proportion governing the starting point
  K - The number of terms to compute in the infinite sum
  Returns:
  The small time representation.
  */

  // Initialize output
  double out = 0.0;

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  if ( abs(v) < machine_double ) // zero drift case
    out = Fs0_lower( t, a, w, K );

  double S1 = 0.0; double S2 = 0.0;
  double sqt = pow(t,0.5);
  double sgn = 1.0; if ( v < 0.0 ) sgn = -1.0;

  for (int k = K; k > 0; k--) {
    double p1 = -sgn*(2.0*a*k+a*w+v*t)/sqt;
    double p2 = sgn*(2.0*a*k+a*w-v*t)/sqt;
    S1 = S1 + exp_pnorm( 2.0*v*a*k, p1 ) -
      exp_pnorm( -2.0*v*a*k - 2.0*v*a*w, p2 );
    p1 = sgn*(2.0*a*k-a*w-v*t)/sqt;
    p2 = -sgn*(2.0*a*k-a*w+v*t)/sqt;
    S2 = S2 + exp_pnorm( -2.0*v*a*k, p1 ) -
      exp_pnorm( 2.0*v*a*k - 2.0*v*a*w, p2 );
  }
  double p3 = R::pnorm(-sgn*(a*w+v*t)/sqt,0.0,1.0,1,0);
  double p4 = exp_pnorm(-2.0*v*a*w, sgn*(a*w-v*t)/sqt);

  out = Pu(v, a, w) + sgn*( (p3 - p4) + S1 + S2 );

  return( out );
}

// Lookup - 08

inline double F_lower( double t, double v, double a, double w,
                       double sigma, double epsilon ) {
  /*
  Purpose:
  Computes the value for the lower sub-distribution.
  Arguments:
  t       - A decision time
  v       - A drift rate
  a       - Boundary separation
  w       - The proportion governing the starting point
  sigma   - The coefficient of drift
  epsilon - The degree of precision
  Returns:
  The value for the lower sub-distribution.
  */

  a = a/sigma;
  v = v/sigma;
  int K_l = K_large(t, v, a, w, epsilon);
  int K_s = K_small(t, v, a, w, epsilon);
  double out = 0.0;
  if ( K_l < 10*K_s) {
    out = Fl_lower( t, v, a, w, K_l );
  } else {
    out = Fs_lower( t, v, a, w, K_s );
  }

  return( out );
}

#endif // End include guard
