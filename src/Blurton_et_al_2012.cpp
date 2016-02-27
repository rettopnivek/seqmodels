#include <RcppParallel.h>
#include <Rcpp.h>  // Includes certain libraries of functions
#include <math.h>
#include <limits>
#include "gsl/include/gsl_integration.h" // Numerical integration
#include "gsl/include/gsl_errno.h" // Error handling
#include "miscfunctions.h" // Linear interpolation

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
Lookup - 01:  exp_pnorm
Lookup - 02:  K_large
Lookup - 03:  K_small
Lookup - 04:  Pu
Lookup - 05:  Fl_lower
Lookup - 06:  Fs0_lower
Lookup - 07:  Fs_lower
Lookup - 08:  F_lower
Lookup - 09:  pwiener_scl
Lookup - 10:
Lookup - 11:
Lookup - 12:
Lookup - 13:
Lookup - 14:
Lookup - 15:
Lookup - 16:
Lookup - 17:
Lookup - 18:
Lookup - 19:
Lookup - 20:
Lookup - 21:
Lookup - 22:

### TO DO ###
Update references
Update index
Incorporate all variability into pdiffWrapper

*/

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
  double e = exp( -2.0*v*a*( 1.0-w ) );
  if (e == R_PosInf)
    out = 1.0;
  if ( std::abs( e - 1.0 ) >= machine_double ) // drift isn't near zero and w isn't near 1
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
  double xi = par[4]; double tau = par[5];
  double sigma = par[6]; double eps = par[10];

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

// Lookup - 10
// A scalar version that calculates the function to
// integrate over when the drift rate varies according
// to a normal distribution

double pw_vxi_scl( double x, void * params) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate normal density values
  double p2 = R::dnorm(x,par[4],par[7],0);

  // Calculate wiener density
  par[4] = x;
  double p1 = pwiener_scl( par );
  out = p1*p2;

  return out;
}

// Lookup - 11
// A scalar version that calculates the function to
// integrate over when the starting point varies according
// to a uniform distribution

double pw_vtheta_scl( double x, void * params) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate uniform density values
  double lb = par[3] - par[8]/2.0;
  double ub = par[3] + par[8]/2.0;

  double p2 = R::dunif(x,lb,ub,0);

  // Calculate wiener density
  par[3] = x;
  double p1 = pwiener_scl( par );
  out = p1*p2;

  return out;
}

// Lookup - 12
// A scalar version that calculates the function to
// integrate over when the residual latency varies according
// to a uniform distribution

double pw_vtau_scl( double x, void * params) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate uniform density values
  double lb = par[5] - par[9]/2.0;
  double ub = par[5] + par[9]/2.0;

  double p2 = R::dunif(x,lb,ub,0);

  // Calculate wiener density
  par[5] = x;
  double p1 = pwiener_scl( par );
  out = p1*p2;

  return out;
}

// Lookup - 13
// Numerical integration routine for a single parameter with variability

double pw_var(std::vector<double> par, double a,double b, int ver ) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( (ver > 1) && (ver < 5) ) {
    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    if ( ver == 2 ) F.function = &pw_vxi_scl;
    if ( ver == 3 ) F.function = &pw_vtheta_scl;
    if ( ver == 4 ) F.function = &pw_vtau_scl;
    F.params = &par;

    // When integrating over eta, the algorithm can be quite slow for
    // small times. Therefore, ff the cdf is very tiny for the regular
    // wiener process, don't integrate.
    int runIntegration = 1;
    if ( ver == 2 ) {
      double tmp = pwiener_scl( par );
      if ( tmp < 1e-5 ) {
        runIntegration = 0;
        result = tmp;
      }
    }
    if ( runIntegration == 1) {
      gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,
                            w, &result, &error);
    }

    // Free up memory
    gsl_integration_workspace_free (w);
  }

  // Check for illegal values
  if (result < 0.0) result = 0.0;

  return(result);
}

// Lookup - 14
// A scalar version that calculates the function to
// integrate over when both drift and starting point vary

double pw_vxi_vtheta_scl( double x, void * params ) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate uniform density values
  double lb = par[3] - par[8]/2.0;
  double ub = par[3] + par[8]/2.0;

  double p2 = R::dunif(x,lb,ub,0);

  // Calculate wiener density with variable drift
  par[3] = x;
  double p1 = pw_var( par, -5.0*par[7] + par[4],
                      5.0*par[7] + par[4], 2 );
  out = p1*p2;

  return out;
}

// Lookup - 15
// A scalar version that calculates the function to
// integrate over when both drift and residual latency vary

double pw_vxi_vtau_scl( double x, void * params ) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate uniform density values
  double lb = par[5] - par[9]/2.0;
  double ub = par[5] + par[9]/2.0;

  double p2 = R::dunif(x,lb,ub,0);

  // Calculate wiener density with variable drift
  par[5] = x;
  double p1 = pw_var( par, -5.0*par[7] + par[4],
                      5.0*par[7] + par[4], 2 );
  out = p1*p2;

  return out;
}

// Lookup - 16
// A scalar version that calculates the function to
// integrate over when both starting point and residual latency vary

double pw_vtheta_vtau_scl( double x, void * params ) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate uniform density values
  double lb = par[5] - par[9]/2.0;
  double ub = par[5] + par[9]/2.0;

  double p2 = R::dunif(x,lb,ub,0);

  // Calculate wiener density with variable starting point
  par[5] = x;

  // Calculate uniform density values
  lb = par[3] - par[8]/2.0;
  ub = par[3] + par[8]/2.0;

  double p1 = pw_var( par, lb,
                      ub, 3 );
  out = p1*p2;

  return out;
}

// Lookup - 17
// Numerical integration routine for two parameters with variability

double pw_var2(std::vector<double> par, double a, double b, int ver ) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( (ver > 4) && (ver < 8) ) {
    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    if ( ver == 5 ) F.function = &pw_vxi_vtheta_scl;
    if ( ver == 6 ) F.function = &pw_vxi_vtau_scl;
    if ( ver == 7 ) F.function = &pw_vtheta_vtau_scl;
    F.params = &par;

    gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,
                          w, &result, &error);

    // Free up memory
    gsl_integration_workspace_free (w);
  }

  // Check for illegal values
  if (result < 0.0) result = 0.0;

  return(result);
}

// Lookup - 18
// A scalar version that calculates the function to
// integrate over when drift, starting point, and residual latency all
// vary

double pw_vxi_vtheta_vtau_scl( double x, void * params ) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate uniform density values
  double lb = par[5] - par[9]/2.0;
  double ub = par[5] + par[9]/2.0;

  double p2 = R::dunif(x,lb,ub,0);

  // Calculate wiener density with variable drift and starting point
  par[5] = x;

  // Calculate uniform density values
  lb = par[3] - par[8]/2.0;
  ub = par[3] + par[8]/2.0;

  double p1 = pw_var2( par, lb,
                       ub, 5 );
  out = p1*p2;

  return out;
}

// Lookup - 19
// Numerical integration routine for three parameters with variability

double pw_var3(std::vector<double> par, double a, double b, int ver ) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( ver == 8 ) {
    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &pw_vxi_vtheta_vtau_scl;
    F.params = &par;

    gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,
                          w, &result, &error);

    // Free up memory
    gsl_integration_workspace_free (w);
  }

  // Check for illegal values
  if (result < 0.0) result = 0.0;

  return(result);
}


// Lookup - 20
// Overall wrapper function for use with parallelization

double piff_wrapper( std::vector<double> par ) {

  // Extract version
  int ver = par[11];
  double out = 0.0;

  // Standard wiener process
  if ( ver == 1 ) {
    out = pwiener_scl( par );
  }
  // Variability in drift
  if (ver == 2) {
    out = pw_var( par, -5.0*par[7] + par[4],
                  5.0*par[7] + par[4], 2 );
  }
  // Variability in starting point
  if (ver == 3) {
    out = pw_var( par, par[3] - par[8]/2.0,
                  par[3] + par[8]/2.0, 3 );
  }
  // Variability in residual latency
  if (ver == 4) {
    out = pw_var( par, par[5] - par[9]/2.0,
                  par[5] + par[9]/2.0, 4 );
  }

  return( out );
}

// Lookup - 21
// RcppParallel worker function

struct pdiffWorker : public Worker
{
  // Input matrix
  const RMatrix<double> input;

  // Destination matrix
  RVector<double> output;

  // initialize with source and destination
  pdiffWorker(const Rcpp::NumericMatrix input,
              Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      output[j] = piff_wrapper( prm );
    }
  }
};

// Lookup - 22
//' @rdname ddiff
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pdiff( Rcpp::NumericVector rt,
                           Rcpp::NumericVector ch,
                           Rcpp::NumericVector alpha,
                           Rcpp::NumericVector theta,
                           Rcpp::NumericVector xi,
                           Rcpp::NumericVector tau,
                           Rcpp::NumericVector eta =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector stheta =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector stau =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector sigma =
                             Rcpp::NumericVector::create(1.0),
                           double eps = 1e-29, int parYes = 1) {

  int N_rt = rt.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_alpha = alpha.size(); // Number of parameters
  int N_theta = theta.size();
  int N_xi = xi.size();
  int N_tau = tau.size();
  int N_eta = eta.size();
  int N_stheta = stheta.size();
  int N_stau = stau.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int rt_inc = 0.0;
  int ch_inc = 0.0;
  int alpha_inc = 0.0;
  int theta_inc = 0.0;
  int xi_inc = 0.0;
  int tau_inc = 0.0;
  int eta_inc = 0.0;
  int stheta_inc = 0.0;
  int stau_inc = 0.0;
  int sigma_inc = 0.0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create( N_rt, N_ch, N_alpha, N_theta,
                                            N_xi, N_tau, N_eta, N_stheta,
                                            N_stau, N_sigma) );

  // If N == 1, no parallelization
  if (N == 1) parYes = 1;

  // Structure of matrix
  // prm(,0) = rt; prm(,1) = ch; prm(,2) = alpha;
  // prm(,3) = theta; prm(,4) = xi; prm(,5) = tau;
  // prm(,6) = sigma; prm(,7) = eta; prm(,8) = stheta
  // prm(,9) = stau; prm(,10) = eps; prm(,11) = ver;

  // Set output vector
  Rcpp::NumericVector output(N);

  // Set input matrix for parameters
  Rcpp::NumericMatrix input(N,12);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    input(nv,0) = rt(rt_inc);
    input(nv,1) = ch(ch_inc);
    input(nv,2) = alpha(alpha_inc);
    input(nv,3) = theta(theta_inc);
    input(nv,4) = xi(xi_inc);
    input(nv,5) = tau(tau_inc);
    input(nv,6) = sigma(sigma_inc);
    input(nv,7) = eta(eta_inc);
    // Fix small values of eta to 0.0
    if ( input(nv,7) < 1e-4 ) { input(nv,7) = 0.0; }
    input(nv,8) = stheta(stheta_inc);
    // Fix small values of stheta to 0.0
    if ( input(nv,8) < 1e-4 ) { input(nv,8) = 0.0; }
    input(nv,9) = stau(stau_inc);
    // Fix small values of stau to 0.0
    if ( input(nv,9) < 1e-4 ) { input(nv,9) = 0.0; }
    input(nv,10) = eps;

    input(nv,11) = 1.0;
    if ( (input(nv,7) > 0.0) &&
         (input(nv,8) == 0.0) &&
         (input(nv,9) == 0.0) ) input(nv,11) = 2.0;
    if ( (input(nv,7) == 0.0) &&
         (input(nv,8) > 0.0) &&
         (input(nv,9) == 0.0) ) input(nv,11) = 3.0;
    if ( (input(nv,7) == 0.0) &&
         (input(nv,8) == 0.0) &&
         (input(nv,9) > 0.0) ) input(nv,11) = 4.0;

    if ( (input(nv,7) > 0.0) &&
         (input(nv,8) > 0.0) &&
         (input(nv,9) == 0.0) ) input(nv,11) = 5.0;
    if ( (input(nv,7) > 0.0) &&
         (input(nv,8) == 0.0) &&
         (input(nv,9) > 0.0) ) input(nv,11) = 6.0;
    if ( (input(nv,7) == 0.0) &&
         (input(nv,8) > 0.0) &&
         (input(nv,9) > 0.0) ) input(nv,11) = 7.0;

    if ( (input(nv,7) > 0.0) &&
         (input(nv,8) > 0.0) &&
         (input(nv,9) > 0.0) ) input(nv,11) = 8.0;

    if ( (input(nv,7) > 10.0) ||
         ( input(nv,3) - (input(nv,8)/2.0 ) <= 0.0 ) ||
         ( input(nv,3) + (input(nv,8)/2.0 ) >= 1.0 ) ||
         ( input(nv,5) - (input(nv,9)/2.0 ) < 0.0 ) )
      input(nv,11) = 0.0;

    rt_inc = rt_inc + 1;
    ch_inc = ch_inc + 1;
    alpha_inc = alpha_inc + 1;
    theta_inc = theta_inc + 1;
    xi_inc = xi_inc + 1;
    tau_inc = tau_inc + 1;
    eta_inc = eta_inc + 1;
    stheta_inc = stheta_inc + 1;
    stau_inc = stau_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_rt==rt_inc) rt_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_theta==theta_inc) theta_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_tau==tau_inc) tau_inc = 0;
    if (N_eta==eta_inc) eta_inc = 0;
    if (N_stheta==stheta_inc) stheta_inc = 0;
    if (N_stau==stau_inc) stau_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Calculate likelihood
  if (parYes == 0) {

    for (int j = 0; j < N; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      output(j) = piff_wrapper( prm );
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    pdiffWorker mt(input, output);

    // Call parallelFor to do the work
    parallelFor(0, N, mt);
  }

  return( output );
}

// Lookup - 23
// Scalar function that generates random choices based on
// absorption at one of the two boundaries for a wiener
// process

double rwiener_choice( std::vector<double> par ) {

  // Initialize output
  double out;

  // Calculate probability of absorption
  // at lower boundary
  double prb = Pu( par[4], par[2], par[3] );

  // Check for inadmissable values and determine choice
  if ( ( par[2] <= 0.0 ) ||
       ( par[3] <= 0.0 ) ||
       ( par[3] >= 1.0 ) ||
       ( par[5] < 0.0 )  ||
       ( par[6] <= 0.0 ) ||
       ( prb <= 0.0 ) ) {
    out = NA_REAL;
  } else {
    out = 0.0;
    double u = R::runif(0.0,1.0);
    if ( u > prb ) out = 1.0;
  }

  return( out );
}

// Lookup - 24
// Scalar function that generates random deviates from the
// two-boundary wiener process by approximating the inverse
// cdf via linear interpolation

double rwiener_scl( std::vector<double> par ) {

  double mxRT = par[11];
  int em_stop = par[12];
  double err = par[13];

  // Initialize response time
  double cur_t = 0.0;

  // Check for inadmissable values
  if ( (par[1] == 0.0 || par[1] == 1.0) ) {

    double maxP = Pu( par[4], par[2], par[3] );
    // Extract choice probability
    if ( par[1] == 1.0 ) {
      maxP = 1.0 - maxP;
    }

    // Define an initial set of RTs
    std::vector<double> iRT(5);
    for (int i = 1; i < 5; i++) {
      iRT[i] = ( exp(i)/exp(5) )*mxRT;
    }

    // Determine the associated CDF values
    std::vector<double> iPrb(5);
    for (int i = 0; i < 5; i++) {
      par[0] = iRT[i];
      iPrb[i] = pwiener_scl( par );
    }

    // Draw from a uniform distribution between 0 and 1
    double p = R::runif(0.0,maxP);
    // Rcpp::Rcout << "p: " << p << std::endl; // For debugging

    // Determine the initial window that the point falls between
    int p0 = minMax( p, iPrb, 0 );
    int p1 = minMax( p, iPrb, 1 );

    double lbp = iPrb[p0]; double lbt = iRT[p0];
    double ubp = iPrb[p1]; double ubt = iRT[p1];

    double prev_t = ubt; double prev_prb = ubp;
    cur_t = linInterp( p, lbp, ubp, lbt, ubt );
    par[0] = cur_t;
    double prb = pwiener_scl( par );
    // Rcpp::Rcout << "cur_t: " << cur_t << std::endl; // For debugging
    // Rcpp::Rcout << "prb: " << prb << std::endl; // For debugging

    double epsilon = 1.0;
    int cnt = 0;

    while ( (cnt < em_stop) & (epsilon > err) ) {
      if (prb < p) {
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
      cur_t = linInterp( p, lbp, ubp, lbt, ubt );
      par[0] = cur_t;
      prb = pwiener_scl( par );

      cnt = cnt + 1;
      epsilon = std::abs( p - prb );

    }

  } else {
    cur_t = NA_REAL;
  }

  return( cur_t + par[5] );
}

// Lookup - 25
// RcppParallel worker function

struct rdiffWorker : public Worker
{
  // Input matrix
  const RMatrix<double> input;

  // Destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  rdiffWorker(const Rcpp::NumericMatrix input,
              Rcpp::NumericMatrix output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(14);

      for (int i = 0; i < 14; i++) { prm[i] = input(j,i); }

      // Determine choice
      double ch_l = rwiener_choice( prm );
      output(j,1) = ch_l;

      // Determine response time
      prm[1] = ch_l;
      double rt_l = rwiener_scl( prm );
      output(j,0) = rt_l;

    }
  }
};


// Lookup - 26
//' @rdname ddiff
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rdiff( int N,
                           Rcpp::NumericVector alpha,
                           Rcpp::NumericVector theta,
                           Rcpp::NumericVector xi,
                           Rcpp::NumericVector tau,
                           Rcpp::NumericVector eta =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector stheta =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector stau =
                             Rcpp::NumericVector::create(0.0),
                           Rcpp::NumericVector sigma =
                             Rcpp::NumericVector::create(1.0),
                           double eps = 1e-29,
                           double mxRT = 4.0,
                           double em_stop = 20,
                           double err = .0001,
                           int parYes = 1) {

  int N_alpha = alpha.size(); // Number of parameters
  int N_theta = theta.size();
  int N_xi = xi.size();
  int N_tau = tau.size();
  int N_eta = eta.size();
  int N_stheta = stheta.size();
  int N_stau = stau.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int alpha_inc = 0.0;
  int theta_inc = 0.0;
  int xi_inc = 0.0;
  int tau_inc = 0.0;
  int eta_inc = 0.0;
  int stheta_inc = 0.0;
  int stau_inc = 0.0;
  int sigma_inc = 0.0;

  // If N == 1, no parallelization
  if (N == 1) parYes = 1;

  // Structure of matrix
  // prm(,0) = rt; prm(,1) = ch; prm(,2) = alpha;
  // prm(,3) = theta; prm(,4) = xi; prm(,5) = tau;
  // prm(,6) = sigma; prm(,7) = eta; prm(,8) = stheta
  // prm(,9) = stau; prm(,10) = eps; prm(,11) = ver;

  // Set output vector
  Rcpp::NumericMatrix output(N,2);

  // Set input matrix for parameters
  Rcpp::NumericMatrix input(N,14);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    // Parameters for trial-to-trial variability

    input(nv,7) = eta(eta_inc);
    // Fix small values of eta to 0.0
    if ( input(nv,7) < 1e-4 ) { input(nv,7) = 0.0; }
    input(nv,8) = stheta(stheta_inc);
    // Fix small values of stheta to 0.0
    if ( input(nv,8) < 1e-4 ) { input(nv,8) = 0.0; }
    input(nv,9) = stau(stau_inc);
    // Fix small values of stau to 0.0
    if ( input(nv,9) < 1e-4 ) { input(nv,9) = 0.0; }

    // Generate random values for parameters
    input(nv,3) = R::runif(
      theta(theta_inc) - input(nv,8)/2.0,
      theta(theta_inc) + input(nv,8)/2.0 );
    input(nv,4) = R::rnorm( xi(xi_inc),
          input(nv,7) );
    input(nv,5) = R::runif(
      tau(tau_inc) - input(nv,9)/2.0,
      tau(tau_inc) + input(nv,9)/2.0 );

    // Parameters that don't vary
    input(nv,0) = 0.0;
    input(nv,1) = 0.0;
    input(nv,2) = alpha(alpha_inc);
    input(nv,6) = sigma(sigma_inc);
    input(nv,10) = eps;
    input(nv,11) = mxRT;
    input(nv,12) = em_stop;
    input(nv,13) = err;

    alpha_inc = alpha_inc + 1;
    theta_inc = theta_inc + 1;
    xi_inc = xi_inc + 1;
    tau_inc = tau_inc + 1;
    eta_inc = eta_inc + 1;
    stheta_inc = stheta_inc + 1;
    stau_inc = stau_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_theta==theta_inc) theta_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_tau==tau_inc) tau_inc = 0;
    if (N_eta==eta_inc) eta_inc = 0;
    if (N_stheta==stheta_inc) stheta_inc = 0;
    if (N_stau==stau_inc) stau_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Calculate likelihood
  if (parYes == 0) {

    for (int j = 0; j < N; j++) {

      std::vector<double> prm(14);

      for (int i = 0; i < 14; i++) { prm[i] = input(j,i); }

      // Determine choice
      double ch_l = rwiener_choice( prm );
      output(j,1) = ch_l;

      // Determine response time
      prm[1] = ch_l;
      double rt_l = rwiener_scl( prm );
      output(j,0) = rt_l;
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    rdiffWorker mt(input, output);

    // Call parallelFor to do the work
    parallelFor(0, N, mt);
  }

  colnames(output) =
    Rcpp::CharacterVector::create("rt", "ch");

  return( output );
}



/*


*/
