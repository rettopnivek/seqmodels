// Include guard to protect against multiple definitions error
#ifndef __NAVARROFUSS__
#define __NAVARROFUSS__

#include <Rcpp.h>
#include "wpparamverify.h"
#include "miscfunctions.h"

/*
Purpose:
Implementation of the adaptive algorithm of Navarro and
Fuss (2009) for computing the joint density of the
two boundary wiener process.

References:
Navarro, D. J., & Fuss, I. G. (2009). Fast and
  accurate calculations for first-passage times in
  Wiener diffusion models. Journal of Mathematical
  Psychology, 53, 222-230.

Index:
Lookup - 01:  wfpt
Lookup - 02:  pr_absorb
*/

// Lookup - 01

inline double wfpt( double t, double v, double a,
                    double z, double err ) {
  /*
  Purpose:
  Adaptively calculates the joint PDF for the Wiener
  process based on whether the small or large time
  version is more appropriate (Navarro & Fuss, 2009).
  Arguments:
  t   - A decision time
  v   - A drift rate
  a   - Boundary separation
  z   - Starting point
  err - Precision
  Returns:
  The joint density for the two-boundary wiener process.
  */

  // Variable declaration
  double p1;
  std::vector<double> v1(2);

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  // Use normalized time
  double norm_t = t/pow(a,2.0);

  // Convert to relative start point
  double w = z/a;

  // Calculate number of terms in infinite sum needed for larger times

  // The minimum bound
  double minBound = ceil( 1.0/M_PI/pow(norm_t, 0.5) );

  double lambda = M_PI*norm_t*err;

  // if error threshold was set too high set to boundary condition
  int kl = minBound;
  int ks = 2;

  // If the error threshold is set low enough
  if (lambda < 1.0) {

    // Find appropriate number of iterations
    v1[0] = machine_double;
    v1[1] = M_PI*norm_t*err;
    p1 = pow( -2.0*log( std::max(v1[0],v1[1]) )/pow(M_PI,2.0)/norm_t, 0.5);
    v1[0] = minBound; v1[1] = ceil(p1);
    kl = std::max(v1[0],v1[1]); // Ensure the boundary conditions are met

  }

  // Calculate number of terms in infinite sum needed for smaller times
  lambda = 2.0*pow( 2.0*M_PI*norm_t, 0.5 )*err;

  // If the error threshold is set low enough
  if ( lambda < 1.0 ) {

    // Find appropriate number of iterations
    v1[0] = machine_double; v1[1] = 2.0*err*pow(2.0*M_PI*norm_t,0.5);
    p1 = ceil( 2.0 + pow(-2.0*norm_t*log( std::max(v1[0],v1[1]) ), 0.5) );
    v1[0] = pow(norm_t,0.5)+1.0; v1[1] = p1;
    ks = std::max(v1[0],v1[1]); // Ensure the boundary conditions are met

  }

  // Compute the density for the standard case, f( norm_t | 0, 1, w )
  double den = 0.0; // Initialize density

  if (ks < kl) { // if small t is better...

    // Summands
    std::vector<double> n(ks);
    double lb = -floor((ks-1)/2);

    int inc = 0;
    while (inc < ks) {
      n[inc] = lb + inc;
      inc = inc + 1;
    }

    // Scaling term
    double scl = 1.0/pow(a,2.0)*exp(-v*a*w - pow(v,2.0)*t/2.0);

    for (int k = 0; k < ks; k++) {
      v1[0] = machine_double;
      v1[1] = exp( -1.0*pow( w + 2.0*n[k], 2.0 )/2.0/norm_t);
      p1 = ( w + 2.0*n[k] )*std::max(v1[0],v1[1]);
      den = den + p1;
    }
    den = 1.0/pow( 2.0*M_PI*pow(norm_t,3.0), 0.5)*den*scl;

  } else { // If large t is better

    // Scaling term
    double scl = 1.0/pow(a,2.0)*exp(-v*a*w - pow(v,2.0)*t/2.0);

    for (int k = 0; k < kl; k++) {
      double n = k + 1.0;
      p1 = n*exp( -1.0*pow( n, 2.0 )*pow( M_PI, 2.0 )*norm_t/2.0 );
      den = den + p1*sin(n*M_PI*w);
    }

    v1[0] = M_PI*den*scl; v1[1] = 0.0;
    den = std::max(v1[0],v1[1]);

  }

  return( den );
}

// Lookup - 02

inline double pr_absorb( double a, double z, double v,
                         double s ) {
  /*
  Purpose:
  Calculates the probability of the path reaching
  the lower bound (Navarro & Fuss, 2009).
  Arguments:
  a - Boundary separation
  z - Starting point
  v - Drift rate
  s - Coefficient of drift
  Returns:
  Probability that accumulation path will reach the lower
  boundary.
  */

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  double out = 0.0;

  // Calculate bias
  double w = z/a;
  // If drift rate equals 0
  if (std::abs(v) < machine_double) out = 1.0 - w;
  // Otherwise...
  else {
    double p1 = exp(-2.0*v*a/pow(s,2.0));
    double p2 = exp(-2.0*v*z/pow(s,2.0));
    out = std::min(1.0, (p1-p2)/(p1-1.0));
  }

  return out;
}

#endif // End include guard
