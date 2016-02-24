#include <RcppParallel.h>
#include <Rcpp.h>  // Includes certain libraries of functions
#include <math.h>
#include <limits>
#include "gsl/include/gsl_integration.h" // Numerical integration
#include "gsl/include/gsl_errno.h" // Error handling
#include "Misc_functions.h" // Error handling

using namespace RcppParallel;

/*
Purpose:
Assorted functions for the calculation of the probability density
function of the Wald Race model (Logan et al., 2014) as well as
its simulation.

Notes:
Forthcoming

References:
Logan, G. D., Van Zandt, T., Verbruggen, F., & Wagenmakers, E. J. (2014).
  On the ability to inhibit thought and action: General and special theories
  of an act of control. Psychological Review, 121(1), 66-95.
  doi:10.1037/a0035230.
Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating random
  variates using transformations with multiple roots. The American
  Statistician, 30(2), 88-90. doi:10.2307/2683801.

Index
Lookup - 01:  rinvgauss_scl
Lookup - 02:  rinvgauss
Lookup - 03:  dinvgauss_scl
Lookup - 04:  dinvgauss
Lookup - 05:  pinvgauss_scl
Lookup - 06:  pinvgauss
Lookup - 07:  rwaldrace
Lookup - 08:  dwaldrace_scl
Lookup - 09:  dwaldrace
Lookup - 10:  int_dwaldrace_scl
Lookup - 11:  pwaldrace_scl
Lookup - 12:  pwaldraceWorker
Lookup - 13:  pwaldrace
*/

// Lookup - 01
// A scalar version to be used in the rinvgauss function

double rinvgauss_scl(double kappa, double xi, double sigma) {

  // Variable declarations
  double out;

  // Convert into standard parameterization for
  // the wald distribution
  double mu = kappa/xi; double lambda = pow(kappa,2.0)/pow(sigma,2.0);

  // Check for illegal values
  if ( (mu<=0.0) | (lambda<=0.0) | (sigma<=0.0) ) out = NA_REAL;
  else {
    // Generate a random draw
    double v = R::rnorm(0.0,1.0);
    double z = R::runif(0.0,1.0);
    double y = pow(v,2.0);
    double p1 = mu + pow(mu,2.0)*y/(2.0*lambda);
    double p2 = mu/(2.0*lambda);
    double p3 = pow(4.0*mu*lambda*y + pow(mu*y,2.0),0.5);
    double x = p1+p2*p3;
    if (z<=mu/(mu+x)) out = x;
    else out = pow(mu,2.0)/x;
  }

  return(out);
}

// Lookup - 02
//' Random Deviates from the Wald Distribution
//'
//' Generates random draws from a Wald distribution (i.e. an
//' inverse gaussian), parameterized based on Brownian
//' motion.
//'
//' @param N the number of draws
//' @param kappa a vector of thresholds determining when a decision
//'   terminates (kappa > 0).
//' @param xi a vector of drift rates, or rates of evidence accumulation
//'   (xi > 0).
//' @param sigma a vector of the within-trial variabilities
//'   (sigma > 0).
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @section References:
//' Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating
//'   random variates using transformations with multiple roots.
//'   The American Statistician, 30 (2), 88-90. doi:10.2307/2683801.
//'
//' @return Returns a numeric vector consisting of random draws from
//'   the Wald distribution.
//'
//' @examples
//' rinvgauss(8,c(1,2,-1,1),c(1,.5,1,-1),c(1,1,1,1))
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rinvgauss(int N, Rcpp::NumericVector kappa,
                              Rcpp::NumericVector xi,
                              Rcpp::NumericVector sigma) {

  int N_kappa = kappa.size(); // Number of parameters
  int N_xi = xi.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int kappa_inc = 0;
  int xi_inc = 0;
  int sigma_inc = 0;

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector kappa_v(N);
  Rcpp::NumericVector xi_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    kappa_v(nv) = kappa(kappa_inc);
    xi_v(nv) = xi(xi_inc);
    sigma_v(nv) = sigma(sigma_inc);

    kappa_inc = kappa_inc + 1;
    xi_inc = xi_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_kappa==kappa_inc) kappa_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Generate draws
  for (int n = 0; n < N; n++) {
    out(n) = rinvgauss_scl( kappa_v(n), xi_v(n), sigma_v(n) );
  }

  return( out );
}

// Lookup - 03
// A scalar version for the dinvgauss function

double dinvgauss_scl (double t, double kappa, double xi,
                      double sigma, int ln = 0) {

  // Variable declaration
  double p1; double p2;

  // To store output
  double out;

  // Check for illegal values
  if ( (t <= 0.0) | (kappa <= 0.0) | ( xi <= 0.0) )
    out = log(0.0);
  else {
    // Calculate the log-density
    p1 = pow( ( pow(kappa,2.0)/pow(sigma,2.0) )/( 2*M_PI*pow( t, 3.0 ) ), 0.5 );
    p2 = exp( -1.0*( pow(kappa,2.0)/pow(sigma,2.0) )*pow( t-( kappa/xi ), 2.0 )/( 2.0*pow( ( kappa/xi ), 2.0 )*t ) );
    out = log(p1) + log(p2);
  }

  // If specified return the density
  if (ln==0) out = exp(out);

  return( out );
}

// Lookup - 04
//' Density Function for the Wald Distribution
//'
//' Calculates the density function for the Wald distribution (i.e.
//' an inverse gaussian), parameterized based on Brownian
//' motion.
//'
//' @param kappa a vector of thresholds determining when a decision
//'   terminates (kappa > 0).
//' @param xi a vector of drift rates, or rates of evidence accumulation
//'   (xi > 0).
//' @param sigma a vector of the within-trial variabilities
//'   (sigma > 0).
//' @param ln indicates whether the log-likelihood should be returned,
//'   where 1 = True, 0 = False (the default).
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @return Returns a numeric vector consisting of the density (or
//'   log density) for the model.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dinvgauss(Rcpp::NumericVector t,
                              Rcpp::NumericVector kappa,
                              Rcpp::NumericVector xi,
                              Rcpp::NumericVector sigma, int ln = 0) {

  int N_t = t.size(); // Number of observations
  int N_kappa = kappa.size(); // Number of parameters
  int N_xi = xi.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int t_inc = 0;
  int kappa_inc = 0;
  int xi_inc = 0;
  int sigma_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_t, N_kappa,
                                           N_xi, N_sigma) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector t_v(N);
  Rcpp::NumericVector kappa_v(N);
  Rcpp::NumericVector xi_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    t_v(nv) = t(t_inc);
    kappa_v(nv) = kappa(kappa_inc);
    xi_v(nv) = xi(xi_inc);
    sigma_v(nv) = sigma(sigma_inc);

    t_inc = t_inc + 1;
    kappa_inc = kappa_inc + 1;
    xi_inc = xi_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_t==t_inc) t_inc = 0;
    if (N_kappa==kappa_inc) kappa_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Calculate the density
  for (int n = 0; n < N; n++) {
    out(n) = dinvgauss_scl( t_v(n), kappa_v(n), xi_v(n), sigma_v(n), ln);
  }

  return( out );
}

// Lookup-05
// A scalar version for the pinvgauss function

double pinvgauss_scl (double t, double kappa, double xi, double sigma) {

  // Variable declaration
  double p1; double p2;

  // To store output
  double out;

  // Check for illegal values
  if ( (t <= 0.0) | (kappa <= 0.0) | ( xi <= 0.0) | (sigma <= 0.0) )
    out = 0.0;
  else {
    // Calculate the cumulative density
    p1 = pow( ( pow(kappa,2.0)/pow(sigma,2.0) )/t,0.5 )*(t/(kappa/xi)-1.0);
    p2 = -1.0*pow( ( pow(kappa,2.0)/pow(sigma,2.0) )/t,0.5 )*(t/(kappa/xi)+1.0);
    out = R::pnorm(p1,0.0,1.0,1,0) +
      exp(2.0*( pow(kappa,2.0)/pow(sigma,2.0) )/(kappa/xi))*R::pnorm(p2,0.0,1.0,1,0);
  }

  return( out );
}

// Lookup - 06
//' Distribution Function for the Wald Distribution
//'
//' Calculates the distribution function for the Wald distribution (i.e.
//' an inverse gaussian), parameterized based on Brownian
//' motion.
//'
//' @param kappa a vector of thresholds determining when a decision
//'   terminates (kappa > 0).
//' @param xi a vector of drift rates, or rates of evidence accumulation
//'   (xi > 0).
//' @param sigma a vector of the within-trial variabilities
//'   (sigma > 0).
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @return Returns a numeric vector giving P(T <= t | kappa, xi, sigma).
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pinvgauss(Rcpp::NumericVector t,
                              Rcpp::NumericVector kappa,
                              Rcpp::NumericVector xi,
                              Rcpp::NumericVector sigma) {

  int N_t = t.size(); // Number of observations
  int N_kappa = kappa.size(); // Number of parameters
  int N_xi = xi.size();
  int N_sigma = sigma.size();

  // Increment variables for loop
  int t_inc = 0;
  int kappa_inc = 0;
  int xi_inc = 0;
  int sigma_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_t, N_kappa, N_xi, N_sigma) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors in Armadillo for stable parallel processing
  Rcpp::NumericVector t_v(N);
  Rcpp::NumericVector kappa_v(N);
  Rcpp::NumericVector xi_v(N);
  Rcpp::NumericVector sigma_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    t_v(nv) = t(t_inc);
    kappa_v(nv) = kappa(kappa_inc);
    xi_v(nv) = xi(xi_inc);
    sigma_v(nv) = sigma(sigma_inc);

    t_inc = t_inc + 1;
    kappa_inc = kappa_inc + 1;
    xi_inc = xi_inc + 1;
    sigma_inc = sigma_inc + 1;
    if (N_t==t_inc) t_inc = 0;
    if (N_kappa==kappa_inc) kappa_inc = 0;
    if (N_xi==xi_inc) xi_inc = 0;
    if (N_sigma==sigma_inc) sigma_inc = 0;
  }

  // Calculate the distribution function
  for (int n = 0; n < N; n++) {
    out(n) = pinvgauss_scl( t_v(n), kappa_v(n), xi_v(n), sigma_v(n) );
  }

  return( out );
}

// Lookup - 07
//' Random Deviates from the Wald Race Model
//'
//' Simulates a set of response times and choices from a two accumulator
//' version of the Wald race model (Logan et al., 2014).
//'
//' @param N the number of observations to simulate.
//' @param k1 the threshold determining when a decision terminates for
//'   choices == 1 ( k1 > 0).
//' @param xi1 the average rate of evidence accumulation within a trial
//'   for choices == 1 (xi1 > 0).
//' @param tau1 the residual latency for choices == 1 (tau1 >= 0).
//' @param k0 the threshold determining when a decision terminates for
//'   choices == 0 ( k0 > 0).
//' @param xi0 the average rate of evidence accumulation within a trial
//'   for choices == 0 (xi0 > 0).
//' @param tau0 the residual latency for choices == 0 (tau0 >= 0).
//' @param s0 the within trial variability for choices == 0 (s0 > 0).
//' @param s1 the within trial variability for choices == 1 (s1 > 0).
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled. Inadmissible values
//' return NA.
//'
//' @return Returns a matrix with two columns, the first consiting of
//' response times, the second indicate the choice (i.e. the winning
//' accumulator, either 1 or 0).
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rwaldrace (int N, Rcpp::NumericVector k1,
                               Rcpp::NumericVector xi1,
                               Rcpp::NumericVector tau1,
                               Rcpp::NumericVector k0,
                               Rcpp::NumericVector xi0,
                               Rcpp::NumericVector tau0,
                               Rcpp::NumericVector s1 =
                                 Rcpp::NumericVector::create(1.0),
                                 Rcpp::NumericVector s0 =
                                   Rcpp::NumericVector::create(1.0) ) {

  int N_k1 = k1.size(); // Number of parameters
  int N_xi1 = xi1.size();
  int N_tau1 = tau1.size();
  int N_k0 = k0.size();
  int N_xi0 = xi0.size();
  int N_tau0 = tau0.size();
  int N_s1 = s1.size();
  int N_s0 = s0.size();

  // Increment variables for loop
  int k1_inc = 0.0;
  int xi1_inc = 0.0;
  int s1_inc = 0.0;
  int k0_inc = 0.0;
  int xi0_inc = 0.0;
  int s0_inc = 0.0;
  int tau1_inc = 0.0;
  int tau0_inc = 0.0;

  // Set output vector
  Rcpp::NumericMatrix out(N,2);

  // Create vectors for parameters
  Rcpp::NumericVector k1_v(N);
  Rcpp::NumericVector xi1_v(N);
  Rcpp::NumericVector s1_v(N);
  Rcpp::NumericVector k0_v(N);
  Rcpp::NumericVector xi0_v(N);
  Rcpp::NumericVector s0_v(N);
  Rcpp::NumericVector tau1_v(N);
  Rcpp::NumericVector tau0_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    k1_v(nv) = k1(k1_inc);
    xi1_v(nv) = xi1(xi1_inc);
    s1_v(nv) = s1(s1_inc);
    k0_v(nv) = k0(k0_inc);
    xi0_v(nv) = xi0(xi0_inc);
    s0_v(nv) = s0(s0_inc);
    tau1_v(nv) = tau1(tau1_inc);
    tau0_v(nv) = tau0(tau0_inc);

    k1_inc = k1_inc + 1;
    xi1_inc = xi1_inc + 1;
    s1_inc = s1_inc + 1;
    k0_inc = k0_inc + 1;
    xi0_inc = xi0_inc + 1;
    s0_inc = s0_inc + 1;
    tau1_inc = tau1_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_k1==k1_inc) k1_inc = 0;
    if (N_xi1==xi1_inc) xi1_inc = 0;
    if (N_s1==s1_inc) s1_inc = 0;
    if (N_k0==k0_inc) k0_inc = 0;
    if (N_xi0==xi0_inc) xi0_inc = 0;
    if (N_s0==s0_inc) s0_inc = 0;
    if (N_tau1==tau1_inc) tau1_inc = 0;
    if (N_tau0==tau0_inc) tau0_inc = 0;
  }

  // Generate draws
  for (int n = 0; n < N; n++) {

    // Check for inadmissable values and return NA if necessary
    if ( (k1_v(n) <= 0.0) | (xi1_v(n) <= 0.0) |
         (k0_v(n) <= 0.0) | (xi0_v(n) <= 0.0) |
         (s1_v(n) <= 0.0) | (s0_v(n) <= 0.0) |
         (tau1_v(n) < 0.0) | (tau0_v(n) < 0.0) ) {
      out(n,0) = NA_REAL;
      out(n,1) = NA_REAL;
    } else {

      // Simulate values from each accumulator
      double t1 = rinvgauss_scl(k1_v(n),
                                xi1_v(n),s1_v(n));
      double t0 = rinvgauss_scl(k0_v(n),
                                xi0_v(n),s0_v(n));

      double rt = -1.0;
      double ch = -1.0;

      // Determine the winning accumulator
      if (t1 < t0) {
        rt = t1 + tau1_v(n);
        ch = 1;
      }
      if (t1 > t0) {
        rt = t0 + tau0_v(n);
        ch = 0;
      }

      // Save the results
      out(n,0) = rt; out(n,1) = ch;
    }
  }

  return ( out );
}

// Lookup - 08
// A scalar version for the dwaldrace function

double dwaldrace_scl (double rt, double ch, double k1, double xi1,
                      double s1, double tau1, double k0, double xi0,
                      double s0, double tau0, int ln = 0) {

  // If residual latency impacts decisions
  // double dt = rt - ch*tau1 - (1-ch)*tau0;
  // double pt = rt - ch*tau0 - (1-ch)*tau1;

  // If residual latency is separate
  double dt = rt - ch*tau1 - (1-ch)*tau0;
  double pt = rt - ch*tau1 - (1-ch)*tau0;

  double dxi = ch*(xi1) + (1-ch)*xi0;
  double pxi = ch*(xi0) + (1-ch)*xi1;
  double dk = ch*(k1) + (1-ch)*k0;
  double pk = ch*(k0) + (1-ch)*k1;
  double ds = ch*(s1) + (1-ch)*s0;
  double ps = ch*(s0) + (1-ch)*s1;
  double cdf = 1.0;
  double out = log(0.0);

  if ( (k1 > 0.0) & (xi1 > 0.0) &
       (k0 > 0.0) & (xi0 > 0.0) &
       (s1 > 0.0) & (s0 > 0.0) &
       (dt > 0.0) ) {
    if (pt > 0.0) cdf = 1.0-pinvgauss_scl(pt,pk,pxi,ps);
    out = dinvgauss_scl(dt,dk,dxi,ds,1) + log(cdf);
  }
  if (ln==0) out = exp( out );

  return( out );
}

// Lookup - 09
//' Density Function for the Wald Race Model
//'
//' Calculates the joint density for a two choice version of the Wald race
//' model.
//'
//' @param rt a vector of response times (rt > 0).
//' @param ch a vector of choices (ch = {0,1}).
//' @param k1 the threshold determining when a decision terminates for
//'   choices == 1 ( k1 > 0).
//' @param xi1 the average rate of evidence accumulation within a trial
//'   for choices == 1 (xi1 > 0).
//' @param tau1 the residual latency for choices == 1 (tau1 >= 0).
//' @param k0 the threshold determining when a decision terminates for
//'   choices == 0 ( k0 > 0).
//' @param xi0 the average rate of evidence accumulation within a trial
//'   for choices == 0 (xi0 > 0).
//' @param tau0 the residual latency for choices == 0 (tau0 >= 0).
//' @param s0 the within trial variability for choices == 0 (s0 > 0;
//'   default is 1.0).
//' @param s1 the within trial variability for choices == 1 (s1 > 0;
//'   default is 1.0).
//' @param ln indicates whether the log-likelihood should be returned,
//'   where 1 = True, 0 = False (the default).
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @return Returns the joint likelihood (or the log-likelihood) for the
//' two-choice version of the Wald race model.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dwaldrace (Rcpp::NumericVector rt,
                               Rcpp::NumericVector ch,
                               Rcpp::NumericVector k1,
                               Rcpp::NumericVector xi1,
                               Rcpp::NumericVector tau1,
                               Rcpp::NumericVector k0,
                               Rcpp::NumericVector xi0,
                               Rcpp::NumericVector tau0,
                               Rcpp::NumericVector s1 =
                                 Rcpp::NumericVector::create(1.0),
                                 Rcpp::NumericVector s0 =
                                   Rcpp::NumericVector::create(1.0),
                         int ln = 0) {

  int N_rt = rt.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_k1 = k1.size(); // Number of parameters
  int N_xi1 = xi1.size();
  int N_s1 = s1.size();
  int N_tau1 = tau1.size();
  int N_k0 = k0.size();
  int N_xi0 = xi0.size();
  int N_s0 = s0.size();
  int N_tau0 = tau0.size();

  // Increment variables for loop
  int rt_inc = 0.0;
  int ch_inc = 0.0;
  int k1_inc = 0.0;
  int xi1_inc = 0.0;
  int s1_inc = 0.0;
  int tau1_inc = 0.0;
  int k0_inc = 0.0;
  int xi0_inc = 0.0;
  int s0_inc = 0.0;
  int tau0_inc = 0.0;

  // Determine the longest input vector
  int N = max(Rcpp::NumericVector::create(N_rt, N_ch, N_k1,
                                    N_xi1, N_s1, N_tau1, N_k0,
                                    N_xi0, N_s0, N_tau0));

  // Set output vector
  Rcpp::NumericVector out(N);

  // Variable declaration
  Rcpp::NumericVector rt_v(N);
  Rcpp::NumericVector ch_v(N);
  Rcpp::NumericVector k1_v(N);
  Rcpp::NumericVector xi1_v(N);
  Rcpp::NumericVector s1_v(N);
  Rcpp::NumericVector tau1_v(N);
  Rcpp::NumericVector k0_v(N);
  Rcpp::NumericVector xi0_v(N);
  Rcpp::NumericVector s0_v(N);
  Rcpp::NumericVector tau0_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    rt_v(nv) = rt(rt_inc);
    ch_v(nv) = ch(ch_inc);
    k1_v(nv) = k1(k1_inc);
    xi1_v(nv) = xi1(xi1_inc);
    s1_v(nv) = s1(s1_inc);
    tau1_v(nv) = tau1(tau1_inc);
    k0_v(nv) = k0(k0_inc);
    xi0_v(nv) = xi0(xi0_inc);
    s0_v(nv) = s0(s0_inc);
    tau0_v(nv) = tau0(tau0_inc);

    rt_inc = rt_inc + 1;
    ch_inc = ch_inc + 1;
    k1_inc = k1_inc + 1;
    xi1_inc = xi1_inc + 1;
    s1_inc = s1_inc + 1;
    tau1_inc = tau1_inc + 1;
    k0_inc = k0_inc + 1;
    xi0_inc = xi0_inc + 1;
    s0_inc = s0_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_rt==rt_inc) rt_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_k1==k1_inc) k1_inc = 0;
    if (N_xi1==xi1_inc) xi1_inc = 0;
    if (N_s1==s1_inc) s1_inc = 0;
    if (N_tau1==tau1_inc) tau1_inc = 0;
    if (N_k0==k0_inc) k0_inc = 0;
    if (N_xi0==xi0_inc) xi0_inc = 0;
    if (N_s0==s0_inc) s0_inc = 0;
    if (N_tau0==tau0_inc) tau0_inc = 0;
  }

  // Calculate likelihood in parallel
  for (int n = 0; n < N; n++) {
    out(n) = dwaldrace_scl(rt_v(n), ch_v(n), k1_v(n), xi1_v(n),
        s1_v(n), tau1_v(n), k0_v(n), xi0_v(n),
        s0_v(n), tau0_v(n), ln);
  }

  return ( out );
}

// Lookup - 10
// A variant of the scalar version of the density function to
// integrate over.

double int_dwaldrace_scl( double x, void * params) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate the density for the Wald race model
  out = dwaldrace_scl( x, par[0], par[1], par[2],
                       par[3], par[4], par[5], par[6],
                       par[7], par[8], 0 );

  return out;
}

// Lookup - 11
// A function that numerically integrates the density function in
// order to determine the distribution function.

double pwaldrace_scl(std::vector<double> par,
                     double a,double b) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( a >= 0.0 ) {

    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &int_dwaldrace_scl;
    F.params = &par;

    // If the upper boundary is infinite
    if ( b == R_PosInf ) {
      gsl_integration_qagiu (&F, a, 0, 1e-7, 1000,
                             w, &result, &error);
    } else {
      gsl_integration_qags (&F, a, b, 0, 1e-7, 1000,
                            w, &result, &error);
    }

    // Free up memory
    gsl_integration_workspace_free (w);

  }

  // Check for illegal values
  if (result < 0.0) result = 0.0;
  if (result > 1.0) result = 0.0;

  return(result);
}

// Lookup - 12
// RcppParallel worker function

struct pwaldraceWorker : public Worker
{
  // Input matrix
  const RMatrix<double> input;

  // Destination matrix
  RVector<double> output;

  // initialize with source and destination
  pwaldraceWorker(const Rcpp::NumericMatrix input,
              Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      double cur_rt = input(j,0); // Extract RT

      std::vector<double> par(9); // Extract parameters
      for (int i = 1; i < 10; i++) { par[i-1] = input(j,i); }

      output[j] = pwaldrace_scl( par, 0.0, cur_rt );
    }
  }
};

// Lookup - 13
//' Distribution Function for the Wald Race Model
//'
//' Calculates the joint distribution function for a two choice version
//' of the Wald race model.
//'
//' @param rt a vector of response times (rt > 0).
//' @param ch a vector of choices (ch = {0,1}).
//' @param k1 the threshold determining when a decision terminates for
//'   choices == 1 ( k1 > 0).
//' @param xi1 the average rate of evidence accumulation within a trial
//'   for choices == 1 (xi1 > 0).
//' @param tau1 the residual latency for choices == 1 (tau1 >= 0).
//' @param k0 the threshold determining when a decision terminates for
//'   choices == 0 ( k0 > 0).
//' @param xi0 the average rate of evidence accumulation within a trial
//'   for choices == 0 (xi0 > 0).
//' @param tau0 the residual latency for choices == 0 (tau0 >= 0).
//' @param s0 the within trial variability for choices == 0 (s0 > 0;
//'   default is 1.0).
//' @param s1 the within trial variability for choices == 1 (s1 > 0;
//'   default is 1.0).
//' @param parYes if set to 1, the code is run in parallel.
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @return Returns the value(s) for the joint distribution function for
//' the two-choice version of the Wald race model.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pwaldrace ( Rcpp::NumericVector rt,
                                Rcpp::NumericVector ch,
                                Rcpp::NumericVector k1,
                                Rcpp::NumericVector xi1,
                                Rcpp::NumericVector tau1,
                                Rcpp::NumericVector k0,
                                Rcpp::NumericVector xi0,
                                Rcpp::NumericVector tau0,
                                Rcpp::NumericVector s1 =
                                Rcpp::NumericVector::create(1.0),
                                Rcpp::NumericVector s0 =
                                  Rcpp::NumericVector::create(1.0),
                                int parYes = 1 ) {

  int N_rt = rt.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_k1 = k1.size(); // Number of parameters
  int N_xi1 = xi1.size();
  int N_s1 = s1.size();
  int N_tau1 = tau1.size();
  int N_k0 = k0.size();
  int N_xi0 = xi0.size();
  int N_s0 = s0.size();
  int N_tau0 = tau0.size();

  // Increment variables for loop
  int rt_inc = 0.0;
  int ch_inc = 0.0;
  int k1_inc = 0.0;
  int xi1_inc = 0.0;
  int s1_inc = 0.0;
  int tau1_inc = 0.0;
  int k0_inc = 0.0;
  int xi0_inc = 0.0;
  int s0_inc = 0.0;
  int tau0_inc = 0.0;

  // Determine the longest input vector
  int N = max(Rcpp::NumericVector::create(N_rt, N_ch, N_k1,
                                          N_xi1, N_s1, N_tau1, N_k0,
                                          N_xi0, N_s0, N_tau0));

  // Structure of matrix
  // prm(,0) = rt; prm(,1) = ch; prm(,2) = alpha;
  // prm(,3) = theta; prm(,4) = xi; prm(,5) = tau;
  // prm(,6) = sigma; prm(,7) = eta; prm(,8) = stheta
  // prm(,9) = stau; prm(,10) = eps; prm(,11) = ver;

  // Set output vector
  Rcpp::NumericVector output(N);

  // Set input matrix for parameters
  Rcpp::NumericMatrix input(N,10);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    input(nv,0) = rt(rt_inc);
    input(nv,1) = ch(ch_inc);
    input(nv,2) = k1(k1_inc);
    input(nv,3) = xi1(xi1_inc);
    input(nv,4) = s1(s1_inc);
    input(nv,5) = tau1(tau1_inc);
    input(nv,6) = k0(k0_inc);
    input(nv,7) = xi0(xi0_inc);
    input(nv,8) = s0(s0_inc);
    input(nv,9) = tau0(tau0_inc);

    rt_inc = rt_inc + 1;
    ch_inc = ch_inc + 1;
    k1_inc = k1_inc + 1;
    xi1_inc = xi1_inc + 1;
    s1_inc = s1_inc + 1;
    tau1_inc = tau1_inc + 1;
    k0_inc = k0_inc + 1;
    xi0_inc = xi0_inc + 1;
    s0_inc = s0_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_rt==rt_inc) rt_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_k1==k1_inc) k1_inc = 0;
    if (N_xi1==xi1_inc) xi1_inc = 0;
    if (N_s1==s1_inc) s1_inc = 0;
    if (N_tau1==tau1_inc) tau1_inc = 0;
    if (N_k0==k0_inc) k0_inc = 0;
    if (N_xi0==xi0_inc) xi0_inc = 0;
    if (N_s0==s0_inc) s0_inc = 0;
    if (N_tau0==tau0_inc) tau0_inc = 0;
  }

  // Calculate likelihood
  if (parYes == 0) {

    for (int j = 0; j < N; j++) {

      double cur_rt = input(j,0);
      std::vector<double> par(9);
      for (int i = 1; i < 10; i++) { par[i-1] = input(j,i); }

      output(j) = pwaldrace_scl( par, 0.0, cur_rt );
    }
  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    pwaldraceWorker mt(input, output);

    // Call parallelFor to do the work
    parallelFor(0, N, mt);
  }

  return ( output );
}

// Lookup - 14
// A scalar function that calculates the quantile given a cumulative probability
// using linear interpolation.

double qwaldrace_scl( std::vector<double> prm ) {

  double p = prm[0]; // The CDF value to invert
  double mxRT = prm[10];
  int em_stop = prm[11];
  double err = prm[12];
  int joint = prm[13];
  double cur_t = 0.0; // Initialize output

  std::vector<double> par(9);
  for ( int k = 1; k < 10; k++ ) {
    par[k-1] = prm[k];
  }

  // If the quantile is not based on the joint distribution
  if ( joint == 0 ) {
    // Rescale the quantile
    p = p*pwaldrace_scl( par, 0.0, R_PosInf );
  }

  // Return the non-decision component if a quantile of 0 is given
  if ( p == 0.0 ) {
    if ( par[0] == 1.0 ) { cur_t = par[4]; } else { cur_t = par[8]; }
  }

  // Define an initial set of RTs
  std::vector<double> iRT(5);
  for (int i = 1; i < 5; i++) {
    iRT[i] = ( exp(i)/exp(5) )*mxRT;
  }

  // Determine the associated CDF values
  std::vector<double> iPrb(5);
  for (int i = 0; i < 5; i++) {
    iPrb[i] = pwaldrace_scl( par, 0.0, iRT[i] );
  }

  // Determine the initial window that the point falls between
  int p0 = minMax( p, iPrb, 0 );
  int p1 = minMax( p, iPrb, 1 );

  double lbp = iPrb[p0]; double lbt = iRT[p0];
  double ubp = iPrb[p1]; double ubt = iRT[p1];

  double prev_t = ubt; double prev_prb = ubp;
  cur_t = linInterp( p, lbp, ubp, lbt, ubt );
  double prb = pwaldrace_scl( par, 0.0, cur_t );

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
    prb = pwaldrace_scl( par, 0.0, cur_t );

    cnt = cnt + 1;
    epsilon = std::abs( p - prb );

  }

  return( cur_t );
}

// Lookup - 15
// RcppParallel worker function

struct qwaldraceWorker : public Worker
{
  // Input matrix
  const RMatrix<double> input;

  // Destination matrix
  RVector<double> output;

  // initialize with source and destination
  qwaldraceWorker(const Rcpp::NumericMatrix input,
                  Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      std::vector<double> par(14); // Extract parameters
      for (int i = 0; i < 14; i++) { par[i] = input(j,i); }

      output[j] = qwaldrace_scl( par );
    }
  }
};

// Lookup - 16
//' Inverse Distribution Function for the Wald Race Model
//'
//' Calculates the inverse of the joint distribution function for a two choice
//' version of the Wald race model using linear interpolation.
//'
//' @param p a vector of probabilities ( 0 >= p >= 1).
//' @param ch a vector of choices (ch = {0,1}).
//' @param k1 the threshold determining when a decision terminates for
//'   choices == 1 ( k1 > 0).
//' @param xi1 the average rate of evidence accumulation within a trial
//'   for choices == 1 (xi1 > 0).
//' @param tau1 the residual latency for choices == 1 (tau1 >= 0).
//' @param k0 the threshold determining when a decision terminates for
//'   choices == 0 ( k0 > 0).
//' @param xi0 the average rate of evidence accumulation within a trial
//'   for choices == 0 (xi0 > 0).
//' @param tau0 the residual latency for choices == 0 (tau0 >= 0).
//' @param s0 the within trial variability for choices == 0 (s0 > 0;
//'   default is 1.0).
//' @param s1 the within trial variability for choices == 1 (s1 > 0;
//'   default is 1.0).
//' @param mxRT the maximum RT response time value that the algorithm is applied to.
//' @param em_step the maximum number of iterations for the linear interpolation.
//' @param err the desired degree of precision for the linear interpolation.
//' @param joint If 1, indicates that the probabilities are based on the joint
//'   distribution function.
//' @param parYes if set to 1, the code is run in parallel.
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled.
//'
//' @return Returns the quantile(s) for the inverse joint distribution function for
//' the two-choice version of the Wald race model.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qwaldrace ( Rcpp::NumericVector p,
                                Rcpp::NumericVector ch,
                                Rcpp::NumericVector k1,
                                Rcpp::NumericVector xi1,
                                Rcpp::NumericVector tau1,
                                Rcpp::NumericVector k0,
                                Rcpp::NumericVector xi0,
                                Rcpp::NumericVector tau0,
                                Rcpp::NumericVector s1 =
                                  Rcpp::NumericVector::create(1.0),
                                Rcpp::NumericVector s0 =
                                  Rcpp::NumericVector::create(1.0),
                                double mxRT = 4.0,
                                double em_stop = 20.0,
                                double err = 0.001,
                                double joint = 1.0,
                                int parYes = 1 ) {

  int N_p = p.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_k1 = k1.size(); // Number of parameters
  int N_xi1 = xi1.size();
  int N_s1 = s1.size();
  int N_tau1 = tau1.size();
  int N_k0 = k0.size();
  int N_xi0 = xi0.size();
  int N_s0 = s0.size();
  int N_tau0 = tau0.size();

  // Increment variables for loop
  int p_inc = 0.0;
  int ch_inc = 0.0;
  int k1_inc = 0.0;
  int xi1_inc = 0.0;
  int s1_inc = 0.0;
  int tau1_inc = 0.0;
  int k0_inc = 0.0;
  int xi0_inc = 0.0;
  int s0_inc = 0.0;
  int tau0_inc = 0.0;

  // Determine the longest input vector
  int N = max(Rcpp::NumericVector::create(N_p, N_ch, N_k1,
                                          N_xi1, N_s1, N_tau1, N_k0,
                                          N_xi0, N_s0, N_tau0));

  // Structure of matrix
  // prm(,0) = rt; prm(,1) = ch; prm(,2) = alpha;
  // prm(,3) = theta; prm(,4) = xi; prm(,5) = tau;
  // prm(,6) = sigma; prm(,7) = eta; prm(,8) = stheta
  // prm(,9) = stau; prm(,10) = eps; prm(,11) = ver;

  // Set output vector
  Rcpp::NumericVector output(N);

  // Set input matrix for parameters
  Rcpp::NumericMatrix input(N,14);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    input(nv,0) = p(p_inc);
    input(nv,1) = ch(ch_inc);
    input(nv,2) = k1(k1_inc);
    input(nv,3) = xi1(xi1_inc);
    input(nv,4) = s1(s1_inc);
    input(nv,5) = tau1(tau1_inc);
    input(nv,6) = k0(k0_inc);
    input(nv,7) = xi0(xi0_inc);
    input(nv,8) = s0(s0_inc);
    input(nv,9) = tau0(tau0_inc);

    input(nv,10) = mxRT;
    input(nv,11) = em_stop;
    input(nv,12) = err;
    input(nv,13) = joint;

    p_inc = p_inc + 1;
    ch_inc = ch_inc + 1;
    k1_inc = k1_inc + 1;
    xi1_inc = xi1_inc + 1;
    s1_inc = s1_inc + 1;
    tau1_inc = tau1_inc + 1;
    k0_inc = k0_inc + 1;
    xi0_inc = xi0_inc + 1;
    s0_inc = s0_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_p==p_inc) p_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_k1==k1_inc) k1_inc = 0;
    if (N_xi1==xi1_inc) xi1_inc = 0;
    if (N_s1==s1_inc) s1_inc = 0;
    if (N_tau1==tau1_inc) tau1_inc = 0;
    if (N_k0==k0_inc) k0_inc = 0;
    if (N_xi0==xi0_inc) xi0_inc = 0;
    if (N_s0==s0_inc) s0_inc = 0;
    if (N_tau0==tau0_inc) tau0_inc = 0;
  }

  // Calculate likelihood
  if (parYes == 0) {

    for (int j = 0; j < N; j++) {

      std::vector<double> par(14); // Extract parameters
      for (int i = 1; i < 14; i++) { par[i] = input(j,i); }

      output(j) = qwaldrace_scl( par );
    }
  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    qwaldraceWorker mt(input, output);

    // Call parallelFor to do the work
    parallelFor(0, N, mt);
  }

  return ( output );
}

/*
 t = seq(0,2,.001)
prm = c( 1, 1, 4, 1, .2, 1, 1, 1, .2 )
cdf = pwaldrace( t, prm[1], prm[2], prm[3], prm[5], prm[6], prm[7], prm[9] )
x11(); plot( t, cdf, xlab = 'RT (s)', ylab = 'P( T < t )', bty = 'l',
    type = 'l')
prm2 = c( .1, prm, 4.0, 20, .001, 0 )
q = qwaldrace_scl( prm2 )
abline( v = q, col = 'red' )
abline( h = .5, lty = 2 )
abline( v = t[ max( which( cdf < .5 ) ) ], lty = 2 )
qnt = qwaldrace( prb, prm[1], prm[2], prm[3], prm[5], prm[6], prm[7],
                 prm[9], joint = 0 )
abline( v = qnt, col = 'blue' )
*/
