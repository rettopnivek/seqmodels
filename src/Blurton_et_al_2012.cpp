#include <RcppParallel.h>
#include <Rcpp.h>  // Includes certain libraries of functions
#include <math.h>
#include <limits>
#include "gsl/include/gsl_integration.h" // Numerical integration
#include "gsl/include/gsl_errno.h" // Error handling

using namespace RcppParallel;

/*
Purpose:
Rcpp code for calculating the distribution function for Ratcliff's
diffusion model. The distribution function for the wiener process
is determined using the formulas from Blurton et al. (2009). Code to
generated random deviates from the diffusion model is based on the
an approximation to the inverse of the distribution function. The
parameter variability is calculated using adaptive quadrature numerical
integration routines from the GSL QUADPACK.

References:
Ratcliff, R., & Tuerlinckx, F. (2002). Estimating parameters of the
  diffusion model: Approaches to dealing with contaminant reaction
  times and parameter variability. Psychonomic Bulletin & Review, 9,
  438-481.
Galassi, M., et al. (2009). GNU scientific library reference manual.
  Network Theory Ltd, 83.


Index:
Lookup - 01: exp_pnorm
Lookup - 02: K_large
Lookup - 03: K_small
Lookup - 04: Pu
Lookup - 05: Fl_lower
Lookup - 06: Fs0_lower
Lookup - 07: Fs_lower
Lookup - 08: F_lower
Lookup - 09: pwiener_scl
Lookup - 10:
*/

//
// Forthcoming
//

// Lookup - 01
// Calculates exp(a) * pnorm(b) using an approximation by
// Kiani et al. (2008)

double exp_pnorm( double a, double b ) {
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
// Number of terms required for large time representation

int K_large( double t, double v, double a,
             double w, double epsilon ) {

  double sqrtL1 = pow(1.0/t,.5)*(a/M_PI);
  double p1 = -2.0/t*a*a/M_PI/M_PI;
  double p2 = log( (epsilon*M_PI*t/2.0)*(v*v + M_PI*M_PI/a/a) );
  double p3 = v*a*w + v*v*t/2.0;
  double sqrtL2 = pow( std::max( 1.0, p1*p2 + p3 ), 0.5 );

  return( ceil( std::max( sqrtL1, sqrtL2 ) ) );
}

// Lookup - 03
// Number of terms required for small time representation

int K_small( double t, double v, double a,
             double w, double epsilon ) {

  int out;

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  if ( std::abs(v) < machine_double ) {
    double p1 = pow(t,.5)/2.0/a;
    double p2 = R::qnorm( std::max( 0.0,
                                    std::min( 1.0, epsilon/( 2.0-2.0*w ) ) ),
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
// Probability for absorption at upper barrier

double Pu( double v, double a, double w ) {

  // Determine machine precision
  double machine_double = std::numeric_limits<float>::denorm_min();

  // Define variable for output
  double out = 1.0 - w; // drift is near zero or w is near 1
  double e = exp(-2.0 * v * a * ( 1.0-w ) );
  if (e == R_PosInf)
    out = 1.0;
  if ( abs( e - 1.0 ) > machine_double ) // drift isn't near zero and w isn't near 1
    out = (1.0 - e) / (exp(2.0*v*a*w) - e); // standard case

  return( out );
}

// Lookup - 05
// Large time representation of lower subdistribution

double Fl_lower(double t, double v, double a,
                double w, int K) {

  // Define variable for output
  double out = 0.0;

  for (int k = K; k > 0; k--) {
    double oldOut = out;
    double p1 = k/(v*v + k*k*M_PI*M_PI/a/a);
    double p2 = exp(-v*a*w - 1.0/2.0*v*v*t - 1.0/2.0*k*k*M_PI*M_PI/a/a*t);
    double p3 = sin(M_PI*k*w);
    out = oldOut - p1*p2*p3;
  }
  out = Pu(v, a, w) + 2.0*M_PI/a/a * out;

  return( out );
}

// Lookup - 06
// Zero drift version

double Fs0_lower( double t, double a, double w, int K ) {

  double out = 0.0;

  for (int k = K; k >= 0; k--) {
    double p1 = (-2.0*k - 2.0 + w)*a/pow(t,0.5);
    double p2 = (-2.0*k - w)*a/pow(t,0.5);
    out = R::pnorm( p1, 0.0, 1.0, 1, 0 ) + R::pnorm( p2, 0.0, 1.0, 1, 0 );
  }

  return( 2.0*out );
}

// Lookup - 07
// Small time representation of the upper subdistribution

double Fs_lower( double t, double v, double a, double w, int K ) {

  //
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
    S1 = S1 + exp_pnorm( 2.0*v*a*k, p1 ) - exp_pnorm( -2.0*v*a*k - 2.0*v*a*w, p2 );
    p1 = sgn*(2.0*a*k-a*w-v*t)/sqt;
    p2 = -sgn*(2.0*a*k-a*w+v*t)/sqt;
    S2 = S2 + exp_pnorm( -2.0*v*a*k, p1 ) - exp_pnorm( 2.0*v*a*k - 2.0*v*a*w, p2 );
  }
  double p3 = R::pnorm(-sgn*(a*w+v*t)/sqt,0.0,1.0,1,0);
  double p4 = exp_pnorm(-2.0*v*a*w, sgn*(a*w-v*t)/sqt);

  out = Pu(v, a, w) + sgn*( (p3 - p4) + S1 + S2 );

  return( out );
}


// Lookup - 08
// Lower subdistribution

double F_lower( double t, double v, double a, double w,
                double sigma, double epsilon ) {
  //
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

// Lookup - 09
// A scalar versions to calculate the full version of the
// distribution function for the Wiener process

double pwiener_scl( std::vector<double> par ) {

  double rt = par[0]; double ch = par[1];
  double alpha = par[2]; double theta = par[3];
  double xi = par[4]; double sigma = par[5];
  double tau = par[6]; double eps = par[7];

  double out = 0.0;
  // Check for inadmissable values
  if ( ( alpha > 0.0) && (theta >= 0.0) &&
       (theta <= 1.0) && (sigma > 0.0) ) {

    // Adjust parameters based on absorbing boundary
    if (ch==1) {
      theta = 1.0 - theta;
      xi = -xi;
    }

    // If time is set to Inf, return the probability of absorption
    // at the relevant boundary
    if (rt == R_PosInf) {
      out = Pu( xi, alpha, theta);
    } else if ( (rt > 0.0) && (tau >= 0.0) &&
      ( rt - tau > 0.0 ) ) {

      // Decision time
      double dt = rt - tau;

      // Calculate CDf
      out = F_lower( dt, xi, alpha, theta, sigma, eps );
    }

  }

  if (out < 0.0) out = 0.0;

  return( out );
}



// Lookup - ??
//' Distribution function for the Wiener diffusion model
//'
//' Calculates the distribution function for the Wiener diffusion
//' model using the implementation of Blurton et al. (2012).
//'
//' @param rt a vector of responses times ( rt > 0 ).
//' @param ch a vector of accuracy/choice values ( ch = {0,1} ).
//' @param alpha a vector of upper boundaries at which the evidence
//'   accumulation terminations.
//' @param theta a vector of proportions determining the starting
//'   point for the evidence accumulation, where the starting point
//'   zeta = alpha*theta ( 0 >= theta >= 1 ).
//' @param xi a vector of drift rates, the rate of evidence accumulation
//'   ( xi > 0 ).
//' @param tau a vector of residual latencies for the non-decision
//'   component ( tau > 0 ).
//' @param sigma a vector of within-trial variabilities, the drift
//'   coefficient ( sigma > 0; default is 1 ).
//' @param eps the margin of error for the infinite sums being calculated
//'   ( default is 1e-29 ).
//' @param joint indicates whether the conditional ( joint = 0 ) or the
//'   joint ( joint = 1) CDF should be returned ( default is 0 ).
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @section References:
//'
//' @return A vector of values for the cumulative distribution.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pdiff( Rcpp::NumericVector rt,
                           Rcpp::NumericVector ch,
                           Rcpp::NumericVector alpha,
                           Rcpp::NumericVector theta,
                           Rcpp::NumericVector xi,
                           Rcpp::NumericVector tau,
                           Rcpp::NumericVector sigma =
                             Rcpp::NumericVector::create(1.0),
                       double eps = 1e-29, int parYes = 1 ) {

  int N_rt = rt.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_alpha = alpha.size(); // Number of parameters
  int N_theta = theta.size();
  int N_xi = xi.size();
  int N_tau = tau.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int rt_inc = 0.0;
  int ch_inc = 0.0;
  int theta_inc = 0.0;
  int xi_inc = 0.0;
  int alpha_inc = 0.0;
  int sigma_inc = 0.0;
  int tau_inc = 0.0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create( N_rt, N_ch, N_alpha, N_theta,
                                      N_xi, N_tau, N_sigma ) );
  // Set output vector
  Rcpp::NumericVector output(N);

  // Create matrix for parameters
  Rcpp::NumericMatrix input(N,8);

  // Input matrix structure
  // input(,0) = rt; input(,1) = ch;
  // input(,2) = alpha; input(,3) = theta;
  // input(,4) = xi; input(,5) = sigma;
  // input(,6) = tau; input(,7) = eps;

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    input(nv,0) = rt(rt_inc);
    input(nv,1) = ch(ch_inc);
    input(nv,2) = alpha(alpha_inc);
    input(nv,3) = theta(theta_inc);
    input(nv,4) = xi(xi_inc);
    input(nv,5) = sigma(sigma_inc);
    input(nv,6) = tau(tau_inc);
    input(nv,7) = eps;

    rt_inc = rt_inc + 1;
    ch_inc = ch_inc + 1;
    alpha_inc = alpha_inc + 1;
    theta_inc = theta_inc + 1;
    xi_inc = xi_inc + 1;
    tau_inc = tau_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_rt==rt_inc) rt_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_theta==theta_inc) theta_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_tau==tau_inc) tau_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Calculate likelihood
  for (int n = 0; n < N; n++) {

    std::vector<double> prm( input.ncol() );
    for (int i = 0; i < input.ncol(); i++ ) {
      prm[i] = input(n,i);
    }

    output(n) = pwiener_scl( prm );
  }

  return( output );
}

/*
//
// Forthcoming
//

// Lookup - ??
// Forthcoming

int minMax( double comp, Rcpp::NumericVector x, int Direction ) {

  Rcpp::LogicalVector chk( x.size() );
  if (Direction==1) { chk = x > comp; } else { chk = x < comp; }

  int out;
  int stp = 0;
  if ( Direction == 1 ) {
    out = 0;
    int inc = chk.size() - 1;
    while ( (stp==0) & (inc >= 0) ) {
      if (chk(inc) == FALSE) { out = inc+1; stp = 1; }
      inc = inc - 1;
    }
    if ( out == chk.size() ) out = out - 1;
  } else {
    out = chk.size() - 1;
    int inc = 0;
    while ( (stp==0) & (inc < chk.size()) ) {
      if (chk(inc) == FALSE) { out = inc-1; stp = 1; }
      inc = inc + 1;
    }
    if ( out == -1 ) out = 0;
  }

  return( out );
}

// Lookup - ??
// Forthcoming

double linInterp( double yN, double y0, double y1,
                  double x0, double x1 ) {

  double b1 = ( y1 - y0 ) / ( x1 - x0 ); // Slope
  double b0 = y1 - b1*x1; // Intercept
  double num = yN - b0;
  double xN = ( num )/b1;

  return( xN );
}

// Lookup - ??
// Forthcoming

double rwiener_choice( std::vector<double> par ) {

  // Initialize output
  double out;

  // Input vector
  std::vector<double> prm(8);
  prm[0] = R_PosInf; prm[1] = 1.0;
  prm[2] = par[0]; prm[3] = par[1];
  prm[4] = par[2]; prm[5] = par[3];
  prm[6] = par[4]; prm[7] = par[5];

  // Calculate probability of a choice

  double prb = pwiener_scl( prm );

  if (prb == 0.0) {
    out = NA_REAL;
  } else {
    out = 1.0;
    double u = R::runif(0.0,1.0);
    if ( u > prb ) out = 0.0;
  }

  return( out );
}

double rwiener_scl( double ch, double alpha, double theta,
                    double xi, double sigma, double tau,
                    double mxRT = 4, double eps = 1e-29,
                    int em_stop = 20, double err = .001 ) {
  // std::

  // Define an initial set of RTs
  std::vector<double> iRT(5);
  for (int i = 1; i < 5; i++) {
    iRT[i] = ( exp(i)/exp(5) )*mxRT;
  }
  // Determine the associated CDF values
  std::vector<double> iPrb(5);
  for (int i = 0; i < 5; i++) {
    iPrb(i) = pwiener_scl( iRT(i), ch, alpha, theta, xi, sigma,
                           0.0, eps, 0 );
  }

  // Draw from a uniform distribution between 0 and 1
  double rw = R::runif(0.0,1.0);
  // Rcout << "rw: " << rw << std::endl; // For debugging

  // Determine the initial window that the point falls between
  int p0 = minMax( rw, iPrb, 0 );
  int p1 = minMax( rw, iPrb, 1 );

  double lbp = iPrb(p0); double lbt = iRT(p0);
  double ubp = iPrb(p1); double ubt = iRT(p1);

  double prev_t = ubt; double prev_prb = ubp;
  double cur_t = linInterp( rw, lbp, ubp, lbt, ubt );
  double prb = pwiener_scl( cur_t, ch, alpha, theta, xi, sigma, 0.0, eps, 0 );

  double epsilon = 1.0;
  int cnt = 0;

  while ( (cnt < em_stop) & (epsilon > err) ) {
    if (prb < rw) {
      lbp = prb;
      lbt = cur_t;
    } else if ( lbp < prb ) {
      ubp = prb;
      ubt = cur_t;
    } else {
      lbp = prb;
      lbt = cur_t;
      ubp = prev_prb;
      ubt = prev_t;
    }
    prev_t = cur_t; prev_prb = prb;
    cur_t = linInterp( rw, lbp, ubp, lbt, ubt );
    prb = pwiener_scl( cur_t, ch, alpha, theta, xi, sigma, 0.0, eps, 0 );

    cnt = cnt + 1;
    epsilon = abs( rw - prb );

  }
  // Rcout << "prb: " << prb << std::endl; // For debugging

  return( cur_t + tau );
}

// Lookup - ??
//' Test
//'
//' Test.
//'
//' @param ch a vector of accuracy/choice values ( ch = {0,1} ).
//'
//' @return Test.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rwiener( int N, Rcpp::NumericVector alpha, Rcpp::NumericVector theta,
                       Rcpp::NumericVector xi, Rcpp::NumericVector tau,
                       Rcpp::NumericVector sigma = Rcpp::NumericVector::create(1.0),
                       double mxRT = 4, double eps = 1e-29, int em_stop = 20,
                       double err = .001 ) {

  int N_alpha = alpha.size(); // Number of parameters
  int N_theta = theta.size();
  int N_xi = xi.size();
  int N_tau = tau.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int theta_inc = 0.0;
  int xi_inc = 0.0;
  int alpha_inc = 0.0;
  int sigma_inc = 0.0;
  int tau_inc = 0.0;

  // Determine the longest input vector

  // Set output vector
  Rcpp::NumericMatrix out(N,2);
  colnames(out) = Rcpp::CharacterVector::create("rt", "ch");

  // Variable declaration
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector theta_v(N);
  Rcpp::NumericVector xi_v(N);
  Rcpp::NumericVector sigma_v(N);
  Rcpp::NumericVector tau_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    alpha_v(nv) = alpha(alpha_inc);
    theta_v(nv) = theta(theta_inc);
    xi_v(nv) = xi(xi_inc);
    sigma_v(nv) = sigma(sigma_inc);
    tau_v(nv) = tau(tau_inc);

    theta_inc = theta_inc + 1;
    xi_inc = xi_inc + 1;
    alpha_inc = alpha_inc + 1;
    sigma_inc = sigma_inc + 1;
    tau_inc = tau_inc + 1;
    if (N_theta==theta_inc) theta_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
    if (N_tau==tau_inc) tau_inc = 0;
  }

  // Random draws using approximation to inverse CDF
  for (int n = 0; n < N; n++) {
    double ch = rwiener_choice( alpha_v(n), theta_v(n), xi_v(n),
                                sigma_v(n), tau_v(n), eps );
    double rt;
    if ( Rcpp::NumericVector::is_na(ch) ) {
      rt = NA_REAL;
    } else {
      rt = rwiener_scl( ch, alpha_v(n), theta_v(n), xi_v(n),
                        sigma_v(n), tau_v(n), mxRT, eps,
                        em_stop, err );
    }
    out(n,0) = rt; out(n,1) = ch;
  }

  return( out );
}


 library(UtilityFunction)

 a = .8; w = .5; v = 1; rl = .2;
 val = seq(.001,1,length=1000)

 d1 = dwiener( val, 1, a, w, v, rl );
 d0 = dwiener( val, 0, a, w, v, rl );

 xl = lowerUpper( .2, val )
 yl = lowerUpper( .2, c(d1,d0) )
 x11(); plot( xl, yl, type='n', ylab = 'Density', xlab = 'RT(s)',
              bty = 'l' );
 lines( val, d1 ); lines( val, d0, lty = 2 )

 tst = rwiener( 5000, a, w, v, rl );
 e1 = density( tst[ tst[,2]==1,1] )
 e0 = density( tst[ tst[,2]==1,1] )
 e1$y = e1$y*mean( tst[,2] )
 e0$y = e0$y*(1-mean( tst[,2] ) )

 lines( e1$x, e1$y, col='blue'); lines( e0$x, e0$y, lty = 2, col='blue')


 linInterp2 = function( yN, y0, y1, x0, x1 ) {
   b1 = ( y1 - y0 ) / ( x1 - x0 ); # Slope
   b0 = y1 - b1*x1; # Intercept
   num = yN - b0;
   xN = ( num )/b1;
   xN
 }

 iRT = c(0, exp(1:5)/exp(5) )*4
 iPrb = pwiener( iRT, 1, a, w, v, 0, joint = 0)

 val = seq(0,4,length=100)
 tst = pwiener( val, 1, a, w, v, 0, joint = 0)
 x11(); plot( c(0,4), c(0,1), type='n', xlab='RT(s)', ylab=('CDF'), bty = 'l')
 lines(val,tst)

 rw = 0.13158

 lbt = iRT[1]; lbp = iPrb[1];
 ubt = iRT[2]; ubp = iPrb[2];
 cur_t = ubt; prb = ubp;

 points( c(lbt,ubt), c(lbp,ubp), col=c('red','blue'), pch = 19 )
 segments( lbt, lbp, ubt, ubp, lty = 2, col = 'pink' )

 prev_t = cur_t; prev_prb = prb;
 cur_t = linInterp2( rw, lbp, ubp, lbt, ubt )
 prb = pwiener( cur_t, 1, a, w, v, 0 )
 points( cur_t, prb, pch=19, col = 'orange' )
 abline( h=rw )

 if (prb < rw) {
   lbp = prb;
   lbt = cur_t;
 } else if ( lbp < prb ) {
   ubp = prb;
   ubt = cur_t;
 } else {
 lbp = prb;
 lbt = cur_t;
  ubp = prev_prb;
  ubt = prev_t;
 }

*/
