#include <RcppParallel.h>
#include <Rcpp.h>  // Includes certain libraries of functions
#include <math.h>
#include <limits>
#include "gsl/include/gsl_integration.h" // Numerical integration
#include "gsl/include/gsl_errno.h" // Error handling

using namespace RcppParallel;

/*
Purpose:
Rcpp code for calculating the density function for Ratcliff's
diffusion model. The density for the wiener process is determined
using the formulas from Navarro and Fuss (2009). The parameter
variability is calculated using adaptive quadrature numerical
integration routines from the GSL QUADPACK.

References:
Navarro, D. J., & Fuss, I. G. (2009). Fast and
  accurate calculations for first-passage times in
  Wiener diffusion models. Journal of Mathematical
  Psychology, 53, 222-230.
Ratcliff, R., & Tuerlinckx, F. (2002). Estimating parameters of the
   diffusion model: Approaches to dealing with contaminant reaction
   times and parameter variability. Psychonomic Bulletin & Review, 9,
   438-481.
Galassi, M., et al. (2009). GNU scientific library reference manual.
  Network Theory Ltd, 83.

Package development:
library(devtools)
library(roxygen2)
1st step: Create folder to store package
2nd step: In Rstudio, set working directory to package folder
          and type 'devtools::use_rcpp()'. This creates the
          basic framework for the package using Rcpp.
3rd step: Select 'Open project' in the file menu and open
          the newly created project file in the package folder.
4th step: Create a C++ file and include the lines:
          //' @useDynLib seqmodels
          //' @importFrom Rcpp sourceCpp

Example of roxygen2 documentation for Rcpp
 //' Title
 //'
 //' Short description.
 //'
 //' @param variable description of variable.
 //'
 //' @return Function's return value(s).
 //' @export
 // [[Rcpp::export]]

// Printing to R console
Rcpp::Rcout << "Debugging example" << std::endl;

Index
Lookup - 01:  wfpt
Lookup - 02:  pr_absorb
Lookup - 03:  dwiener_scl
Lookup - 04:  dw_vxi_scl
Lookup - 05:  dw_vtheta_scl
Lookup - 06:  dw_vtau_scl
Lookup - 07:  dw_var
Lookup - 08:  dw_vxi_vtheta_scl
Lookup - 09:  dw_vxi_vtau_scl
Lookup - 10:  dw_vtheta_vtau_scl
Lookup - 11:  dw_var2
Lookup - 12:  dw_vxi_vtheta_vtau_scl
Lookup - 13:  dw_var3
Lookup - 14:  diff_wrapper
Lookup - 15:  ddiffWorker
Lookup - 16:  ddiff

### TO DO ###
Update .rd page

*/

//' @useDynLib seqmodels
//' @importFrom Rcpp sourceCpp

// Lookup - 01
// Adaptively calcuates the PDF for the Wiener process
// based on whether the small or large time versions is
// most appropriate

double wfpt( double t, double v, double a,
             double z, double err) {

  // t > 0
  // a > 0
  // 0 >= z <= a

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
// Calculates the probability of the path being absorved at
// the lower bound

double pr_absorb(double a, double z, double v,
                 double s) {

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

// Lookup - 03
// A scalar version to calculate the full version of the
// density for the Wiener process

double dwiener_scl( std::vector<double> par ) {

  double rt = par[0]; double ch = par[1];
  double alpha = par[2]; double theta = par[3];
  double xi = par[4]; double tau = par[5];
  double sigma = par[6]; double eps = par[10];

  // Decision time
  double dt = rt - tau;

  // Determine starting point
  double zeta = alpha*theta;

  // Rescale paramters based on within-trial variance
  alpha = alpha/sigma;
  xi = xi/sigma;
  zeta = zeta/sigma;

  // Adjust parameters based on absorbing boundary
  if (ch==1) {
    zeta = alpha-zeta;
    xi = -1.0*xi;
  }

  // Initialize density
  double out = 0.0;
  // Check for inadmissable values
  if ( (rt > 0.0) & (dt > 0.0) & (tau >= 0.0) & (alpha > 0.0) &
       (theta >= 0.0) & (theta <= 1.0) &
       (sigma > 0.0) ) {
    out = wfpt( dt, xi, alpha, zeta, eps );
  }

  return( out );
}

// Lookup - 04
// A scalar version that calculates the function to
// integrate over when the drift rate varies according
// to a normal distribution

double dw_vxi_scl( double x, void * params) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate normal density values
  double p2 = R::dnorm(x,par[4],par[7],0);

  // Calculate wiener density
  par[4] = x;
  double p1 = dwiener_scl( par );
  out = p1*p2;

  return out;
}

// Lookup - 05
// A scalar version that calculates the function to
// integrate over when the starting point varies according
// to a uniform distribution

double dw_vtheta_scl( double x, void * params) {

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
  double p1 = dwiener_scl( par );
  out = p1*p2;

  return out;
}

// Lookup - 06
// A scalar version that calculates the function to
// integrate over when the residual latency varies according
// to a uniform distribution

double dw_vtau_scl( double x, void * params) {

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
  double p1 = dwiener_scl( par );
  out = p1*p2;

  return out;
}

// Lookup - 07
// Numerical integration routine for a single parameter with variability

double dw_var(std::vector<double> par, double a,double b, int ver ) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( (ver > 1) && (ver < 5) ) {
    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    if ( ver == 2 ) F.function = &dw_vxi_scl;
    if ( ver == 3 ) F.function = &dw_vtheta_scl;
    if ( ver == 4 ) F.function = &dw_vtau_scl;
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

// Lookup - 08
// A scalar version that calculates the function to
// integrate over when both drift and starting point vary

double dw_vxi_vtheta_scl( double x, void * params ) {

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
  double p1 = dw_var( par, -5.0*par[7] + par[4],
                      5.0*par[7] + par[4], 2 );
  out = p1*p2;

  return out;
}

// Lookup - 09
// A scalar version that calculates the function to
// integrate over when both drift and residual latency vary

double dw_vxi_vtau_scl( double x, void * params ) {

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
  double p1 = dw_var( par, -5.0*par[7] + par[4],
                      5.0*par[7] + par[4], 2 );
  out = p1*p2;

  return out;
}

// Lookup - 10
// A scalar version that calculates the function to
// integrate over when both starting point and residual latency vary

double dw_vtheta_vtau_scl( double x, void * params ) {

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

  double p1 = dw_var( par, lb,
                      ub, 3 );
  out = p1*p2;

  return out;
}

// Lookup - 11
// Numerical integration routine for two parameters with variability

double dw_var2(std::vector<double> par, double a, double b, int ver ) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( (ver > 4) && (ver < 8) ) {
    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    if ( ver == 5 ) F.function = &dw_vxi_vtheta_scl;
    if ( ver == 6 ) F.function = &dw_vxi_vtau_scl;
    if ( ver == 7 ) F.function = &dw_vtheta_vtau_scl;
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

// Lookup - 12
// A scalar version that calculates the function to
// integrate over when drift, starting point, and residual latency all
// vary

double dw_vxi_vtheta_vtau_scl( double x, void * params ) {

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

  double p1 = dw_var2( par, lb,
                      ub, 5 );
  out = p1*p2;

  return out;
}

// Lookup - 13
// Numerical integration routine for three parameters with variability

double dw_var3(std::vector<double> par, double a, double b, int ver ) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( ver == 8 ) {
    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &dw_vxi_vtheta_vtau_scl;
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

// Lookup - 14
// Overall wrapper function for use with parallelization

double diff_wrapper( std::vector<double> par ) {

  // Extract version
  int ver = par[11];
  double out = 0.0;

  // Standard wiener process
  if ( ver == 1 ) {
    out = dwiener_scl( par );
  }
  // Variability in drift
  if (ver == 2) {
    out = dw_var( par, -5.0*par[7] + par[4],
                  5.0*par[7] + par[4], 2 );
  }
  // Variability in starting point
  if (ver == 3) {
    out = dw_var( par, par[3] - par[8]/2.0,
                  par[3] + par[8]/2.0, 3 );
  }
  // Variability in residual latency
  if (ver == 4) {
    out = dw_var( par, par[5] - par[9]/2.0,
                  par[5] + par[9]/2.0, 4 );
  }

  // Variability in drift and starting point
  if (ver == 5) {
    out = dw_var2( par, par[3] - par[8]/2.0,
                   par[3] + par[8]/2.0, 5 );
  }
  // Variability in drift and residual latency
  if (ver == 6) {
    out = dw_var2( par, par[5] - par[9]/2.0,
                  par[5] + par[9]/2.0, 6 );
  }
  // Variability in starting point and residual latency
  if (ver == 7) {
    out = dw_var2( par, par[5] - par[9]/2.0,
                  par[5] + par[9]/2.0, 7 );
  }
  // Variability in starting point and residual latency
  if (ver == 8) {
    out = dw_var3( par, par[5] - par[9]/2.0,
                   par[5] + par[9]/2.0, 7 );
  }

  return( out );
}

// Lookup - 15
// RcppParallel worker function

struct ddiffWorker : public Worker
{
  // Input matrix
  const RMatrix<double> input;

  // Destination matrix
  RVector<double> output;

  // initialize with source and destination
  ddiffWorker(const Rcpp::NumericMatrix input,
              Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> prm(12);

      for (int i = 0; i < 12; i++) { prm[i] = input(j,i); }

      output[j] = diff_wrapper( prm );
    }
  }
};

// Lookup - 16
//' The Ratcliff diffusion model
//'
//' Density, distribution, random generation, and quantile functions
//' for the diffusion model (e.g. Ratcliff & Tuerlinckx, 2002).
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
//' @param eta a vector of standard deviations for the trial-to-trial
//'   variability in the drift rate (eta > 0). A value results in no
//'   trial-to-trial variability.
//' @param stheta a vector of the ranges for the uniform variability
//'   about the starting point ( 0 <= theta +/- (stheta/2) <= 1). A
//'   value results in no trial-to-trial variability.
//' @param stau a vector of the ranges for the uniform variability
//'   about the residual latency ( 0 <= tau - (stau/2) ). A
//'   value results in no trial-to-trial variability.
//' @param sigma a vector of within-trial variabilities, the drift
//'   coefficient ( sigma > 0 ).
//' @param eps the margin of error for the infinite sums being calculated
//'   ( default is 1e-29 ).
//' @param ln indicates whether the likelihood (ln = 0 ) or the
//'   log-likelihood ( ln = 1 ) should be calculated ( default is 0 ).
//'
//' @section Details:
//' The density function is based on the implementation of Navarro
//' and Fuss (2009). The distribution function is based on the
//' implementation of Blurton et al. (2012). For parameter variability
//' the functions use numerical integration with adaptive quadrature.
//' For unequal vector lengths, values are recycled.
//'
//' @section References:
//' Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate calculations
//'   for first-passage times in Wiener diffusion models. Journal of
//'   Mathematical Psychology, 53(4), 222-230.
//' Ratcliff, R., & Tuerlinckx, F. (2002). Estimating parameters of the
//'   diffusion model: Approaches to dealing with contaminant reaction
//'   times and parameter variability. Psychonomic Bulletin & Review, 9,
//'   438-481.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector ddiff( Rcpp::NumericVector rt,
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
                           double eps = 1e-29, int ln = 0,
                           int parYes = 1) {

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
         ( input(nv,5) - (input(nv,9)/2.0 ) < 0.0 ) ) {
      input(nv,11) = 0.0;
    }

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

      output(j) = diff_wrapper( prm );
    }

  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    ddiffWorker mt(input, output);

    // Call parallelFor to do the work
    parallelFor(0, N, mt);
  }

  if (ln == 1) output = log(output);

  return( output );
}
