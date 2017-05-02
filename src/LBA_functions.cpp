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
 Assorted functions for the calculation of the density, distribution,
 quantile, and random number generation functions of the Frechet
 distribution and several variants of the LBA model (Brown and
 Heathcote, 2008; Terry et al., 2015).

 References:
 Brown, S. D., & Heathcote, A. (2008). The simplest complete model of
   choice response time: Linear ballistic accumulation. Cognitive
   psychology, 57(3), 153-178.
 Terry, A., Marley, A. A. J., Barnwal, E.-J., Wagenmakers, Heathcote,
   A., & Brown, S. D. (2015). Generalising the drift rate distribution
   for linear ballistic accumulators. Journal of Mathematical
   Psychology, 68, 49 - 58.
 Singman, H., Brown, S., Gretton, M., Heathcote, A., Voss, A.,
   Voss, J., & Terry, A. (2016). rtdists: Response Time
   Distributions.  R package version 0.6-6.

 Index
 Lookup - 01:  rfrechet_scl
 Lookup - 02:  rfrechet
 Lookup - 03:  dfrechet_scl
 Lookup - 04:  dfrechet
 Lookup - 05:  pfrechet_scl
 Lookup - 06:  pfrechet
 Lookup - 07:  qfrechet_scl
 Lookup - 08:  qfrechet
 Lookup - 09:  rlba_1acc_scl
 Lookup - 10:  rlba_1acc
 Lookup - 11:  Z_g
 Lookup - 12:  igamma
 Lookup - 13:  Z_fr
 Lookup - 14:  Z_ln
 Lookup - 15:  plba_1acc_scl
 Lookup - 16:  plba_1acc
 Lookup - 17:  dZ_g
 Lookup - 18:  dZ_fr
 Lookup - 19:  dZ_ln
 Lookup - 20:  dlba_1acc_scl
 Lookup - 21:  dlba_1acc
 Lookup - 22:  qlba_1acc_scl
 Lookup - 23:  qlba_1acc
 Lookup - 24:  rlba
 Lookup - 25:  dlba_scl
 Lookup - 26:  dlba
 Lookup - 27:  ?
 Lookup - 28:  ?
 Lookup - 29:  ?

### TO DO ###
 - Add examples
 - Add functions to dist_plot
 */

// Lookup - 01
// A scalar version of the function to generate random deviates from the
// Frechet distribution

double rfrechet_scl( double alpha, double mu ) {

  // Initialize output
  double x = NA_REAL;

  // Check for valid parameter values
  if ( ( alpha > 0.0 ) & (mu > 0.0 ) ) {

    // Generate a value between 0 and 1
    double u = R::runif( 0.0, 1.0 );

    // Use the inverse probability transform to determine the
    // value of the random deviate from the Frechet distribution
    x = mu*pow( -log( u ), -1.0/alpha );
  }

  return( x );
}

// Lookup - 02
//' The Frechet Distribution
//'
//' Random generation, density, distribution, and quantile functions
//' for the Frechet (inverse Weibull) distribution with a shape and
//' scale parameter.
//'
//' @param N the number of draws for random generation.
//' @param t a vector of quantiles (typically response times).
//' @param alpha a vector of shape parameters (alpha > 0).
//' @param mu a vector of scale parameters (mu > 0).
//' @param ln indicates whether the log-likelihood should be returned,
//'   where 1 = True, 0 = False.
//'
//' @section Details:
//' The Frechet distribution is a special case of the generalized extreme
//' value distribution, the distribution of the maxima of a sequence of
//' i.i.d. random variables with no lower boundary. The current
//' parameterization is intended for applications for a variant of the
//' LBA model (Terry et al., 2015).
//'
//' For unequal vector lengths, values are recycled.
//'
//' @section References:
//' Terry, A., Marley, A. A. J., Barnwal, E.-J., Wagenmakers, Heathcote,
//'   A., & Brown, S. D. (2015). Generalising the drift rate distribution
//'   for linear ballistic accumulators. Journal of Mathematical
//'   Psychology, 68, 49 - 58.
//'
//' @examples
//' # Treatment of illegal values and vectorization
//' set.seed(100)
//' rfrechet( 8, c(1,2,-1,1), c(1,.5,1,-1) ) # Returns NA
//' dfrechet( c(.5, -1), c(1,2,-1,1), c(1,.5,1,-1) ) # Returns 0
//' pfrechet( c(.5, -1), c(1,2,-1,1),c(1,.4,1,-1) ) # Returns 0
//'
//' # Distribution function
//' t = seq(0,2,length=1000)
//' plot( t, pfrechet(t,3,.5), type = 'l', xlab = 'Time', ylab =
//'   'Distribution', bty = 'l', yaxt = 'n' )
//' axis(2,seq(0,1,.5))
//' # Quantile function
//' prb = seq( .1, .9, .2 ) # Probabilities
//' qnt = qfrechet( prb, 3, .5 )
//' segments( rep(0,length(prb)), prb, qnt, prb )
//' segments( qnt,rep(0,length(prb)), qnt, prb )
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rfrechet(int N, Rcpp::NumericVector alpha,
                             Rcpp::NumericVector mu) {

  int N_alpha = alpha.size(); // Number of parameters
  int N_mu = mu.size();

  // Increment variables for loop
  int alpha_inc = 0;
  int mu_inc = 0;

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector mu_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    alpha_v(nv) = alpha(alpha_inc);
    mu_v(nv) = mu(mu_inc);

    alpha_inc = alpha_inc + 1;
    mu_inc = mu_inc + 1;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;

  }

  // Generate draws
  for (int n = 0; n < N; n++) {
    out(n) = rfrechet_scl( alpha_v(n), mu_v(n) );
  }

  return( out );
}

// Lookup - 03
// A scalar version of the density function for the Frechet distribution

double dfrechet_scl( double x, double alpha, double mu, int ln ) {

  // Initialize output
  double out = 0.0;

  // Check for valid parameter values
  if ( ( x >= 0.0) & ( alpha > 0.0 ) & (mu > 0.0 ) ) {

    out = ( alpha/mu )*pow( x/mu, -1.0 - alpha ) *
      exp( -pow( x/mu, -alpha ) );

  }
  if ( ln == 1 ) out = log( out );

  return( out );
}

// Lookup - 04
//' @rdname rfrechet
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dfrechet( Rcpp::NumericVector t,
                              Rcpp::NumericVector alpha,
                              Rcpp::NumericVector mu, int ln = 0) {

  int N_t = t.size(); // Number of observations
  int N_alpha = alpha.size(); // Number of parameters
  int N_mu = mu.size();

  // Increment variables for loop
  int t_inc = 0;
  int alpha_inc = 0;
  int mu_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create( N_t, N_alpha, N_mu ) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector t_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector mu_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    t_v(nv) = t(t_inc);
    alpha_v(nv) = alpha(alpha_inc);
    mu_v(nv) = mu(mu_inc);

    t_inc = t_inc + 1;
    alpha_inc = alpha_inc + 1;
    mu_inc = mu_inc + 1;
    if (N_t==t_inc) t_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;
  }

  // Calculate the density
  for (int n = 0; n < N; n++) {
    out(n) = dfrechet_scl( t_v(n), alpha_v(n), mu_v(n), ln);
  }

  return( out );
}

// Lookup - 05
// A scalar version of the cumulative distribution function for the
// Frechet distribution

double pfrechet_scl( double x, double alpha, double mu ) {

  // Initialize output
  double out = 0.0;

  // Check for valid parameter values
  if ( ( x >= 0.0) & ( alpha > 0.0 ) & (mu > 0.0 ) ) {

    out = exp( -pow( x/mu, -alpha ) );

  }

  // Check values
  if ( ( out < 0.0 ) | ( out > 1.0 ) ) out = 0.0;

  return( out );
}

// Lookup - 06
//' @rdname rfrechet
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector pfrechet( Rcpp::NumericVector t,
                              Rcpp::NumericVector alpha,
                              Rcpp::NumericVector mu ) {

  int N_t = t.size(); // Number of observations
  int N_alpha = alpha.size(); // Number of parameters
  int N_mu = mu.size();

  // Increment variables for loop
  int t_inc = 0;
  int alpha_inc = 0;
  int mu_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_t, N_alpha, N_mu ) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors in Armadillo for stable parallel processing
  Rcpp::NumericVector t_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector mu_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    t_v(nv) = t(t_inc);
    alpha_v(nv) = alpha(alpha_inc);
    mu_v(nv) = mu(mu_inc);

    t_inc = t_inc + 1;
    alpha_inc = alpha_inc + 1;
    mu_inc = mu_inc + 1;
    if (N_t==t_inc) t_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;
  }

  // Calculate the distribution function
  for (int n = 0; n < N; n++) {
    out(n) = pfrechet_scl( t_v(n), alpha_v(n), mu_v(n) );
  }

  return( out );
}

// Lookup - 07
// Scalar version of quantile function for Frechet distribution using
// linear interpolation.

double qfrechet_scl( double p, double alpha, double mu ) {

  // Initialize output
  double out = 0.0;

  if ( ( alpha > 0.0 ) & ( mu > 0.0 ) ) {

    if ( p == 1.0 ) out = R_PosInf;

    if ( ( p > 0.0 ) & ( p < 1.0 ) ) {
      out = mu*pow( -log( p ), -1.0/alpha );
    }

  }

  return( out );
}

// Lookup - 08
//' @rdname rfrechet
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qfrechet( Rcpp::NumericVector p,
                              Rcpp::NumericVector alpha,
                              Rcpp::NumericVector mu ) {

  int N_p = p.size(); // Number of observations
  int N_alpha = alpha.size(); // Number of parameters
  int N_mu = mu.size();

  // Increment variables for loop
  int p_inc = 0;
  int alpha_inc = 0;
  int mu_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create(N_p, N_alpha, N_mu) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors
  Rcpp::NumericVector p_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector mu_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    p_v(nv) = p(p_inc);
    alpha_v(nv) = alpha(alpha_inc);
    mu_v(nv) = mu(mu_inc);

    p_inc = p_inc + 1;
    alpha_inc = alpha_inc + 1;
    mu_inc = mu_inc + 1;
    if (N_p==p_inc) p_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_mu==mu_inc) mu_inc = 0;
  }

  // Calculate the distribution function
  for (int n = 0; n < N; n++) {
    out(n) = qfrechet_scl( p_v(n), alpha_v(n), mu_v(n) );
  }

  return( out );
}

// Lookup - 09
// A scalar version for a single LBA accumulator (and its variants)

double rlba_1acc_scl( double A, double b, double alpha, double beta,
                      int ver ) {

  // Initialize output
  double t = NA_REAL;

  // Check for valid parameter values
  if ( ( A >= 0.0 ) & ( b >= A ) ) {

    // Original LBA model
    if ( (ver == 0) & ( beta > 0.0 ) ) {

      // Restrict to positive drift rates
      double v = -1.0;
      while( v < 0.0 ) {
        v = R::rnorm( alpha, beta );
      }

      double a = R::runif( 0.0, A );
      t = (b-a)/v;

    }

    // LBA model with drift rates sampled from a gamma distribution
    if ( (ver == 1) & ( alpha > 0.0) & ( beta > 0.0 ) ) {

      // Use rate formulation in R
      double v = R::rgamma( alpha, 1.0/beta );
      double a = R::runif( 0.0, A );
      t = (b-a)/v;

    }

    // LBA model with drift rates sampled from a Frechet distribution
    if ( (ver == 2) & ( alpha > 0.0) & ( beta > 0.0 ) ) {

      double v = rfrechet_scl( alpha, beta );
      double a = R::runif( 0.0, A );
      t = (b-a)/v;

    }

    // LBA model with drift rates sampled from a log-normal distribution
    if ( (ver == 3) & ( beta > 0.0 ) ) {

      double v = R::rlnorm( alpha, beta );
      double a = R::runif( 0.0, A );
      t = (b-a)/v;

    }

  }

  return( t );
}

// Lookup - 10
//' A Linear Ballistic Accumulator
//'
//' Random generation, density, and distribution functions
//' for a single linear ballistic accumulator (Brown & Heathcote,
//' 2008; Terry et al., 2015).
//'
//' @param N the number of draws for random generation.
//' @param t a vector of quantiles (typically response times).
//' @param A the upper boundary for the uniformly distributed start
//'   points (A >= 0).
//' @param b the threshold for the accumulator (b >= A).
//' @param alpha a location/shape parameter for the distribution
//'   of drift rates.
//' @param beta a scale parameter for the distribution of drift
//'   rates.
//' @param ln indicates whether the log-likelihood should be returned,
//'   where 1 = True, 0 = False (the default).
//' @param ver indicates which variant of drift rate distribution
//'   for the LBA should be used.
//'
//' @note For a standard LBA with normally distributed drifts,
//' the distribution is truncated so that drift rates will only be
//' positive.
//'
//' @section References:
//' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of
//'   choice response time: Linear ballistic accumulation. Cognitive
//'   psychology, 57(3), 153-178.
//' Terry, A., Marley, A. A. J., Barnwal, E.-J., Wagenmakers, Heathcote,
//'   A., & Brown, S. D. (2015). Generalising the drift rate distribution
//'   for linear ballistic accumulators. Journal of Mathematical
//'   Psychology, 68, 49 - 58.
//' Singman, H., Brown, S., Gretton, M., Heathcote, A., Voss, A.,
//'   Voss, J., & Terry, A. (2016). rtdists: Response Time
//'   Distributions.  R package version 0.6-6.
//'
//' @examples
//' # Forthcoming
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rlba_1acc( int N, Rcpp::NumericVector A,
                               Rcpp::NumericVector b,
                               Rcpp::NumericVector alpha,
                               Rcpp::NumericVector beta,
                               int ver = 0 ) {

  int N_A = A.size(); // Number of parameters
  int N_b = b.size();
  int N_alpha = alpha.size();
  int N_beta = beta.size();

  // Increment variables for loop
  int A_inc = 0;
  int b_inc = 0;
  int alpha_inc = 0;
  int beta_inc = 0;

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector A_v(N);
  Rcpp::NumericVector b_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector beta_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    A_v(nv) = A(A_inc);
    b_v(nv) = b(b_inc);
    alpha_v(nv) = alpha(alpha_inc);
    beta_v(nv) = beta(beta_inc);

    A_inc = A_inc + 1;
    b_inc = b_inc + 1;
    alpha_inc = alpha_inc + 1;
    beta_inc = beta_inc + 1;
    if (N_A==A_inc) A_inc = 0;
    if (N_b==b_inc) b_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_beta==beta_inc) beta_inc = 0;
  }

  // Generate draws
  for (int n = 0; n < N; n++) {
    out(n) = rlba_1acc_scl( A_v(n), b_v(n), alpha_v(n),
        beta_v(n), ver );
  }

  return( out );
}

// Lookup - 11
// The truncated mean for the gamma distribution

double Z_g( double t, double A, double b, double alpha, double beta ) {

  double k = R::gammafn( alpha + 1.0 )/( beta*R::gammafn(alpha) );

  double num = R::pgamma( b/t, alpha + 1.0, 1.0/beta, 1, 0 ) -
    R::pgamma( (b-A)/t, alpha + 1.0, 1.0/beta, 1, 0 );

  // G( b/t, alpha, beta ) - G( (b-A)/t, alpha, beta )
  double denom = R::pgamma( b/t, alpha, 1.0/beta, 1, 0 ) -
    R::pgamma( (b-A)/t, alpha, 1.0/beta, 1, 0 );

  return( k*num/denom );
}

// Lookup - 12
// The incomplete gamma function

double igamma( double x, double a ) {

  double out = R::pgamma( a, x, 1.0, 0, 0 ) * R::gammafn( x );

  return( out );
}

// Lookup - 13
// The truncated mean for the Frechet distribution
// Adapted from code from the 'rtdist' R package

double Z_fr( double t, double A, double b, double alpha, double beta ) {

  double mn = b/t;
  double mx = ( b - A )/t;

  double G_mn = pfrechet_scl( mn, alpha, beta );
  double G_mx = pfrechet_scl( mx, alpha, beta );
  double D = G_mx - G_mn;

  double gam = igamma( 1.0 - ( 1.0/alpha),
                       pow( 1.0/beta * mx, -alpha ) ) -
                         igamma( 1.0 - ( 1.0/alpha ),
                                 pow( 1.0/beta * mn, -alpha ) );

  double out = gam/( 1.0/beta * D );

  return( out );
}

// Lookup - 14
// The truncated mean for the log-normal distribution

double Z_ln( double t, double A, double b, double alpha, double beta ) {

  double s2 = pow( beta, 2.0 );

  double k = exp( alpha + s2/2.0 );

  double num = R::pnorm( log( b/t ), alpha + s2, beta, 1, 0 ) -
    R::pnorm( log( (b-A)/t ), alpha + s2, beta, 1, 0 );

  double denom = R::pnorm( log( b/t ), alpha, beta, 1, 0 ) -
    R::pnorm( log( (b-A)/t ), alpha, beta, 1, 0 );

  double out = k*num/denom;

  return( out );
}

// Lookup - 15
// Scalar function for the CDF of the LBA accumulator

double plba_1acc_scl( double t, double A, double b, double alpha,
                      double beta, int ver ) {

  // Initialize output
  double out = 0.0;

  // Check for valid parameter values
  if ( ( t > 0.0 ) & ( A >= 0.0 ) & ( b >= A ) ) {

    // Standard LBA
    if ( ( ver == 0 ) & ( beta > 0.0 ) ) {

      // For small threshold values
      if ( A < 1e-10 ) {

        double p3 = b/t;
        out = 1.0 - R::pnorm( p3, alpha, beta, 1, 0 );

      } else {

        double p1 = (b - A - t*alpha)/(t*beta);
        double p2 = (b - t*alpha)/(t*beta);
        double p3 = (t*beta)/A;
        double p4 = (b - A - t*alpha)/A;
        double p5 = (b - t*alpha)/A;

        double l1 = 1.0 + p4*R::pnorm(p1,0.0,1.0,1,0) -
          p5*R::pnorm(p2,0.0,1.0,1,0);
        double l2 = p3*R::dnorm(p1,0.0,1.0,0) -
          p3*R::dnorm(p2,0.0,1.0,0);

        // Joint distribution - integrates to proportion of
        // finite response times (i.e. with positive drift rates)
        out = l1 + l2;
      }

      // Rescale to distribution function conditional on
      // positive drift rates
      out = out/( 1.0 - R::pnorm( 0.0, alpha, beta, 1, 0 ) );
    }

    // LBA with gamma-distributed drift rates
    if ( ( ver == 1 ) & ( alpha > 0.0) & ( beta > 0.0 ) ) {

      // Determine truncated mean
      double  Z_t = Z_g(t, A, b, alpha, beta);

      // Calculate distribution function
      double p1 = ( ( t*Z_t - b )/A ) *
        R::pgamma( b/t, alpha, 1.0/beta, 1, 0 );
      double p2 = ( ( b - A - t*Z_t )/A ) *
        R::pgamma( ( b - A )/t, alpha, 1.0/beta, 1, 0 );

      out = 1.0 + p1 + p2;
    }

    // LBA with Frechet-distributed drift rates
    if ( ( ver == 2 ) & ( alpha > 0.0 ) & ( beta > 0.0 ) ) {

      // Determine truncated mean
      double  Z_t = Z_fr(t, A, b, alpha, beta);

      // Calculate distribution function
      double p1 = ( ( t*Z_t - b )/A ) *
        pfrechet_scl( b/t, alpha, beta );
      double p2 = ( ( b - A - t*Z_t )/A ) *
        pfrechet_scl( ( b - A )/t, alpha, beta );

      out = 1.0 + p1 + p2;
    }

    // LBA with Frechet-distributed drift rates
    if ( ( ver == 3 ) & ( beta > 0.0 ) ) {

      // Determine truncated mean
      double  Z_t = Z_ln(t, A, b, alpha, beta);

      // Calculate distribution function
      double p1 = ( ( t*Z_t - b )/A ) *
        R::plnorm( b/t, alpha, beta, 1, 0 );
      double p2 = ( ( b - A - t*Z_t )/A ) *
        R::plnorm( ( b - A )/t, alpha, beta, 1, 0 );

      out = 1.0 + p1 + p2;
    }

  }

  // Check for incorrect values
  if ( ( out > 1.0 ) | ( out < 0.0 ) ) out = 0.0;

  return( out );
}

// Lookup - 16
//' @rdname rlba_1acc
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector plba_1acc( Rcpp::NumericVector t,
                               Rcpp::NumericVector A,
                               Rcpp::NumericVector b,
                               Rcpp::NumericVector alpha,
                               Rcpp::NumericVector beta,
                               int ver = 0 ) {

  int N_t = t.size(); // Number of observations
  int N_A = A.size(); // Number of parameters
  int N_b = b.size();
  int N_alpha = alpha.size();
  int N_beta = beta.size();

  // Increment variables for loop
  int t_inc = 0;
  int A_inc = 0;
  int b_inc = 0;
  int alpha_inc = 0;
  int beta_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create( N_t, N_A, N_b,
                                            N_alpha, N_beta ) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector t_v(N);
  Rcpp::NumericVector A_v(N);
  Rcpp::NumericVector b_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector beta_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    t_v(nv) = t(t_inc);
    A_v(nv) = A(A_inc);
    b_v(nv) = b(b_inc);
    alpha_v(nv) = alpha(alpha_inc);
    beta_v(nv) = beta(beta_inc);

    t_inc = t_inc + 1;
    A_inc = A_inc + 1;
    b_inc = b_inc + 1;
    alpha_inc = alpha_inc + 1;
    beta_inc = beta_inc + 1;
    if (N_t==t_inc) t_inc = 0;
    if (N_A==A_inc) A_inc = 0;
    if (N_b==b_inc) b_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_beta==beta_inc) beta_inc = 0;
  }

  // Calculate distribution function
  for (int n = 0; n < N; n++) {
    out(n) = plba_1acc_scl( t_v(n), A_v(n), b_v(n), alpha_v(n),
        beta_v(n), ver );
  }

  return( out );
}

// Lookup - 17
// The derivative of the truncated mean for the gamma distribution

double dZ_g( double t, double A, double b, double alpha, double beta ) {

  double k = R::gammafn( alpha + 1.0 )/( beta * R::gammafn( alpha ) );
  double f = R::pgamma( b/t, alpha + 1.0, 1.0/beta, 1, 0 ) -
    R::pgamma( (b-A)/t, alpha + 1.0, 1.0/beta, 1, 0 );
  double fd = -1.0 * R::dgamma( b/t, alpha + 1.0, 1.0/beta, 0 ) *
    b * pow(t,-2.0) + R::dgamma( (b-A)/t, alpha + 1.0, 1.0/beta, 0 ) *
    ( b - A ) * pow(t,-2.0);
  double h = R::pgamma( b/t, alpha, 1.0/beta, 1, 0 ) -
    R::pgamma( (b-A)/t, alpha, 1.0/beta, 1, 0 );
  double hd = -1.0 * R::dgamma( b/t, alpha, 1.0/beta, 0 ) *
    b * pow(t,-2.0) + R::dgamma( (b-A)/t, alpha, 1.0/beta, 0 ) *
    ( b - A ) * pow(t,-2.0);

  double out = k*( ( fd * h - f * hd )/pow( h, 2.0 ) );

  return( out );
}

// Lookup - 18
// The derivative of the truncated mean for the frechet distribution
// Adapted from code from the 'rtdist' package

double dZ_fr( double t, double A, double b, double alpha, double beta ) {

  double mn = b/t;
  double mx = ( b - A )/t;

  double G_mn = pfrechet_scl( mn, alpha, beta );
  double G_mx = pfrechet_scl( mx, alpha, beta );
  double D = G_mx - G_mn;

  double g_mn = dfrechet_scl( mn, alpha, beta, 0 );
  double g_mx = dfrechet_scl( mx, alpha, beta, 0 );
  double diffG1 = -b/pow( t, 2.0 ) * g_mn;
  double diffG2 = -(b-A)/pow( t, 2.0 ) * g_mx;
  double diffD = diffG1 - diffG2;

  double gam = igamma( 1.0 - ( 1.0/alpha),
                       pow( 1.0/beta * mx, -alpha ) ) -
                         igamma( 1.0 - ( 1.0/alpha ),
                                 pow( 1.0/beta * mn, -alpha ) );

  double subterm = pow( t, -alpha + 2.0 );
  double term1_1 = pow( 1.0/beta * b, -alpha + 1.0 )/subterm;
  double term1_2 = exp( -pow( 1.0/beta * mn, -alpha ) );
  double term1 = -alpha * term1_1 * term1_2;

  double term2_1 = pow( 1.0/beta * ( b - A ), -alpha + 1.0 )/subterm;
  double term2_2 = exp( -pow( 1.0/beta * mx, -alpha ) );
  double term2 = -alpha * term2_1 * term2_2;

  double diffgam = term1 - term2;

  double out = -pow( 1.0/beta, -1.0 ) * (
    -pow( D, -2.0 ) * diffD  * gam +
      ( diffgam * pow( D, -1.0 ) )
  );

  return( out );
}

// Lookup - 19
// The derivative for the truncated mean of a log-normal distribution

double dZ_ln( double t, double A, double b, double alpha, double beta ) {

  // Variable declaration
  double subterm1;
  double subterm2;

  double mn = log( (b-A)/t );
  double mx = log( b/t );

  subterm1 = R::pnorm( ( mx - alpha - pow( beta, 2.0) )/beta,
                       0.0, 1.0, 1, 0 );
  subterm2 = R::pnorm( ( mn - alpha - pow( beta, 2.0 ) )/beta,
                       0.0, 1.0, 1, 0 );
  double u = subterm1 - subterm2;

  subterm1 = R::pnorm( ( mx - alpha )/beta,
                       0.0, 1.0, 1, 0 );
  subterm2 = R::pnorm( ( mn - alpha )/beta,
                       0.0, 1.0, 1, 0 );
  double v = subterm1 - subterm2;

  double k = -1.0/( beta*t );
  subterm1 = k * R::dnorm( ( mx - alpha - pow( beta, 2.0 ) )/beta,
                           0.0, 1.0, 0 );
  subterm2 = k * R::dnorm( ( mn - alpha - pow( beta, 2.0 ) )/beta,
                           0.0, 1.0, 0 );
  double udash = subterm1 - subterm2;

  subterm1 = k * R::dnorm( ( mx - alpha )/beta,
                           0.0, 1.0, 0 );
  subterm2 = k * R::dnorm( ( mn - alpha )/beta,
                           0.0, 1.0, 0 );
  double vdash = subterm1 - subterm2;

  double cnst = exp( alpha + pow( beta, 2.0 )/2.0 );

  double out = ( ( udash*v - vdash*u)/pow( v, 2.0 ) )*cnst;

  return( out );
}

// Lookup - 20
// Scalar function for the PDF of the LBA accumulator

double dlba_1acc_scl( double t, double A, double b, double alpha,
                      double beta, int ver, int ln ) {

  // Initialize output
  double out = 0.0;

  // Check for valid parameter values
  if ( ( t > 0.0 ) & ( A >= 0.0 ) & ( b >= A ) ) {

    // Standard LBA
    if ( ( ver == 0 ) & ( beta > 0.0 ) ) {

      if ( A < 1e-10 ) {

        // For small threshold values
        double phi2 = R::dnorm( b/t, alpha, beta, 0 );
        out = ( b/pow(t,2.0) )*phi2;

      } else {

        // Variable declaration
        double ts = t*beta;
        double tv = t*alpha;
        double p1 = b - tv;
        double p2 = p1/ts;
        double p3 = (p1-A)/ts;

        double Phi2 = R::pnorm(p2,0.0,1.0,1,0);
        double Phi3 = R::pnorm(p3,0.0,1.0,1,0);
        double phi2 = R::dnorm(p2,0.0,1.0,0);
        double phi3 = R::dnorm(p3,0.0,1.0,0);

        // Calculate the likelihood
        out = ( alpha*(Phi2-Phi3) + beta*( phi3-phi2) )/A;

        // Rescale to distribution function conditional on
        // positive drift rates
        out = out/( 1.0 - R::pnorm( 0.0, alpha, beta, 1, 0 ) );

      }
    }

    // LBA with gamma-distributed drift rates
    if ( ( ver == 1 ) & ( alpha > 0.0) & ( beta > 0.0 ) ) {

      // Determine truncated mean
      double Z_t = Z_g( t, A, b, alpha, beta );
      double Z_td = dZ_g( t, A, b, alpha, beta );

      double g_bt = R::dgamma( b/t, alpha, 1.0/beta, 0 );
      double g_bAt = R::dgamma( (b-A)/t, alpha, 1.0/beta, 0 );
      double G_bt = R::pgamma( b/t, alpha, 1.0/beta, 1, 0 );
      double G_bAt = R::pgamma( (b-A)/t, alpha, 1.0/beta, 1, 0 );

      double p1 = -( g_bt * b )/pow( t, 2.0 ) * ( t * Z_t - b )/A;
      double p2 = G_bt * ( Z_t + t * Z_td )/A;
      double p3 = -g_bAt * ( b - A ) * pow( t, -2.0 ) *
        ( b - A - t * Z_t ) / A;
      double p4 = G_bAt * ( -Z_t - t * Z_td )/A;

      out = p1 + p2 + p3 + p4;
    }

    // LBA with frechet-distributed drift rates
    if ( ( ver == 2 ) & ( alpha > 0.0) & ( beta > 0.0 ) ) {

      // Determine truncated mean
      double  Z_t = Z_fr( t, A, b, alpha, beta );
      double Z_td = dZ_fr( t, A, b, alpha, beta );

      double g_bt = dfrechet_scl( b/t, alpha, beta, 0 );
      double g_bAt = dfrechet_scl( (b-A)/t, alpha, beta, 0 );
      double G_bt = pfrechet_scl( b/t, alpha, beta );
      double G_bAt = pfrechet_scl( (b-A)/t, alpha, beta );

      double p1 = -( g_bt * b )/pow( t, 2.0 ) * ( t * Z_t - b )/A;
      double p2 = G_bt * ( Z_t + t * Z_td )/A;
      double p3 = -g_bAt * ( b - A ) * pow( t, -2.0 ) *
        ( b - A - t * Z_t ) / A;
      double p4 = G_bAt * ( -Z_t - t * Z_td )/A;

      out = p1 + p2 + p3 + p4;
    }

    // LBA with lognormal-distributed drift rates
    if ( ( ver == 3 ) & ( beta > 0.0 ) ) {

      // Determine truncated mean
      double  Z_t = Z_ln( t, A, b, alpha, beta );
      double Z_td = dZ_ln( t, A, b, alpha, beta );

      double g_bt = R::dlnorm( b/t, alpha, beta, 0 );
      double g_bAt = R::dlnorm( (b-A)/t, alpha, beta, 0 );
      double G_bt = R::plnorm( b/t, alpha, beta, 1, 0 );
      double G_bAt = R::plnorm( (b-A)/t, alpha, beta, 1, 0 );

      double p1 = -( g_bt * b )/pow( t, 2.0 ) * ( t * Z_t - b )/A;
      double p2 = G_bt * ( Z_t + t * Z_td )/A;
      double p3 = -g_bAt * ( b - A ) * pow( t, -2.0 ) *
        ( b - A - t * Z_t ) / A;
      double p4 = G_bAt * ( -Z_t - t * Z_td )/A;

      out = p1 + p2 + p3 + p4;
    }

  }

  // Check for incorrect values
  if ( out < 0.0 ) out = 0.0;

  // If log-density is desired
  if ( ln == 1 ) out = log( out );

  return( out );
}

// Lookup - 21
//' @rdname rlba_1acc
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dlba_1acc( Rcpp::NumericVector t,
                               Rcpp::NumericVector A,
                               Rcpp::NumericVector b,
                               Rcpp::NumericVector alpha,
                               Rcpp::NumericVector beta,
                               int ver = 0,
                               int ln = 0 ) {

  int N_t = t.size(); // Number of observations
  int N_A = A.size(); // Number of parameters
  int N_b = b.size();
  int N_alpha = alpha.size();
  int N_beta = beta.size();

  // Increment variables for loop
  int t_inc = 0;
  int A_inc = 0;
  int b_inc = 0;
  int alpha_inc = 0;
  int beta_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create( N_t, N_A, N_b,
                                            N_alpha, N_beta ) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector t_v(N);
  Rcpp::NumericVector A_v(N);
  Rcpp::NumericVector b_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector beta_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    t_v(nv) = t(t_inc);
    A_v(nv) = A(A_inc);
    b_v(nv) = b(b_inc);
    alpha_v(nv) = alpha(alpha_inc);
    beta_v(nv) = beta(beta_inc);

    t_inc = t_inc + 1;
    A_inc = A_inc + 1;
    b_inc = b_inc + 1;
    alpha_inc = alpha_inc + 1;
    beta_inc = beta_inc + 1;
    if (N_t==t_inc) t_inc = 0;
    if (N_A==A_inc) A_inc = 0;
    if (N_b==b_inc) b_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_beta==beta_inc) beta_inc = 0;
  }

  // Calculate distribution function
  for (int n = 0; n < N; n++) {
    out(n) = dlba_1acc_scl( t_v(n), A_v(n), b_v(n), alpha_v(n),
        beta_v(n), ver, ln );
  }

  return( out );
}

// Lookup - 22
// Scalar version of quantile function for one LBA accumulator
// distribution using linear interpolation.

double qlba_1acc_scl( double p, double A, double b,
                      double alpha, double beta, int ver,
                      double mxT, int em_stop,
                      double err ) {

  // Initialize output
  double cur_t = 0.0;

  if ( p > 0.0 ) {

    // Define an initial set of values
    std::vector<double> iT(5);
    for (int i = 1; i < 5; i++) {
      iT[i] = iT[i-1] + mxT/5;
    }

    // Determine the associated CDF values
    std::vector<double> iPrb(5);
    for (int i = 0; i < 5; i++) {
      iPrb[i] = plba_1acc_scl( iT[i], A, b, alpha, beta, ver );
    }

    // Determine the initial window that the point falls between
    int p0 = minMax( p, iPrb, 0 );
    int p1 = minMax( p, iPrb, 1 );

    double lbp = iPrb[p0]; double lbt = iT[p0];
    double ubp = iPrb[p1]; double ubt = iT[p1];

    double prev_t = ubt; double prev_prb = ubp;
    cur_t = linInterp( p, lbp, ubp, lbt, ubt );
    double prb = plba_1acc_scl( cur_t, A, b, alpha, beta, ver );

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
      prb = plba_1acc_scl( cur_t, A, b, alpha, beta, ver );

      cnt = cnt + 1;
      epsilon = std::abs( p - prb );

    }

  }

  return( cur_t );
}

// Lookup - 23
//' @rdname rlba_1acc
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector qlba_1acc(Rcpp::NumericVector p,
                              Rcpp::NumericVector A,
                              Rcpp::NumericVector b,
                              Rcpp::NumericVector alpha,
                              Rcpp::NumericVector beta,
                              int ver = 0,
                              double mxT = 10.0, int em_stop = 20,
                              double err = .0001 ) {

  int N_p = p.size(); // Number of observations
  int N_A = A.size(); // Number of parameters
  int N_b = b.size();
  int N_alpha = alpha.size();
  int N_beta = beta.size();

  // Increment variables for loop
  int p_inc = 0;
  int A_inc = 0;
  int b_inc = 0;
  int alpha_inc = 0;
  int beta_inc = 0;

  // Determine the longest input vector
  int N = max( Rcpp::NumericVector::create( N_p, N_A, N_b,
                                            N_alpha, N_beta ) );

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for the parameters
  Rcpp::NumericVector p_v(N);
  Rcpp::NumericVector A_v(N);
  Rcpp::NumericVector b_v(N);
  Rcpp::NumericVector alpha_v(N);
  Rcpp::NumericVector beta_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    p_v(nv) = p(p_inc);
    A_v(nv) = A(A_inc);
    b_v(nv) = b(b_inc);
    alpha_v(nv) = alpha(alpha_inc);
    beta_v(nv) = beta(beta_inc);

    p_inc = p_inc + 1;
    A_inc = A_inc + 1;
    b_inc = b_inc + 1;
    alpha_inc = alpha_inc + 1;
    beta_inc = beta_inc + 1;
    if (N_p==p_inc) p_inc = 0;
    if (N_A==A_inc) A_inc = 0;
    if (N_b==b_inc) b_inc = 0;
    if (N_alpha==alpha_inc) alpha_inc = 0;
    if (N_beta==beta_inc) beta_inc = 0;
  }

  // Calculate the distribution function
  for (int n = 0; n < N; n++) {
    out(n) = qlba_1acc_scl( p_v(n), A_v(n), b_v(n),
        alpha_v(n), beta_v(n),
        ver, mxT, em_stop, err );
  }

  return( out );
}

// Lookup - 24
//' The LBA Model
//'
//' Random generation, density, distribution, and quantile functions
//' for a two accumulator version of the LBA model (Brown & Heathcote,
//' 2008; Terry et al., 2015).
//'
//' @param N the number of observations to simulate.
//' @param rt a vector of response times (rt > 0).
//' @param ch a vector of choices (ch = {0,1}).
//' @param A1 the upper boundary for the uniformly distributed start
//'   points (A1 >= 0) for choices == 1.
//' @param b1 the threshold (b1 >= A1) for choices == 1.
//' @param alpha1 a location/shape parameter for the distribution
//'   of drift rates for choices == 1.
//' @param beta1 a scale parameter for the distribution of drift
//'   rates for choices == 1.
//' @param A0 the upper boundary for the uniformly distributed start
//'   points (A0 >= 0) for choices == 0.
//' @param b0 the threshold (b0 >= A0) for choices == 0.
//' @param alpha0 a location/shape parameter for the distribution
//'   of drift rates for choices == 0.
//' @param beta0 a scale parameter for the distribution of drift
//'   rates for choices == 0.
//' @param ln indicates whether the log-likelihood should be returned,
//'   where 1 = True, 0 = False (the default).
//' @param ver indicates which variant of drift rate distribution
//'   for the LBA should be used.
//' @param rl if 1, the residual latency impacts the decision rule.
//' @param mxRT the maximum RT response time value that the algorithm is applied to.
//' @param em_step the maximum number of iterations for the linear interpolation.
//' @param err the desired degree of precision for the linear interpolation.
//' @param joint If 1, indicates that the probabilities are based on the joint
//'   distribution function.
//' @param parYes if set to 1, the code is run in parallel.
//'
//' @section Notes:
//' For unequal vector lengths, values are recycled. For random draws,
//' inadmissible values return NA. For a standard LBA with normally
//' distributed drifts, the distribution is truncated so that drift
//' rates will only be positive.
//'
//' @section References:
//' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of
//'   choice response time: Linear ballistic accumulation. Cognitive
//'   psychology, 57(3), 153-178.
//' Terry, A., Marley, A. A. J., Barnwal, E.-J., Wagenmakers, Heathcote,
//'   A., & Brown, S. D. (2015). Generalising the drift rate distribution
//'   for linear ballistic accumulators. Journal of Mathematical
//'   Psychology, 68, 49 - 58.
//' Singman, H., Brown, S., Gretton, M., Heathcote, A., Voss, A.,
//'   Voss, J., & Terry, A. (2016). rtdists: Response Time
//'   Distributions.  R package version 0.6-6.
//'
//' @examples
//' Forthcoming
//'
//' @export
// [[Rcpp::export]]

Rcpp::NumericMatrix rlba ( int N,
                           Rcpp::NumericVector A1,
                           Rcpp::NumericVector b1,
                           Rcpp::NumericVector alpha1,
                           Rcpp::NumericVector beta1,
                           Rcpp::NumericVector tau1,
                           Rcpp::NumericVector A0,
                           Rcpp::NumericVector b0,
                           Rcpp::NumericVector alpha0,
                           Rcpp::NumericVector beta0,
                           Rcpp::NumericVector tau0,
                           int rl = 0,
                           int ver = 0 ) {

  int N_A1 = A1.size(); // Number of parameters
  int N_b1 = b1.size();
  int N_alpha1 = alpha1.size();
  int N_beta1 = beta1.size();
  int N_tau1 = tau1.size();
  int N_A0 = A0.size();
  int N_b0 = b0.size();
  int N_alpha0 = alpha0.size();
  int N_beta0 = beta0.size();
  int N_tau0 = tau0.size();

  // Increment variables for loop
  int A1_inc = 0.0;
  int b1_inc = 0.0;
  int alpha1_inc = 0.0;
  int beta1_inc = 0.0;
  int tau1_inc = 0.0;
  int A0_inc = 0.0;
  int b0_inc = 0.0;
  int alpha0_inc = 0.0;
  int beta0_inc = 0.0;
  int tau0_inc = 0.0;

  // Set output vector
  Rcpp::NumericMatrix out(N,2);

  // Create vectors for parameters
  Rcpp::NumericVector A1_v(N);
  Rcpp::NumericVector b1_v(N);
  Rcpp::NumericVector alpha1_v(N);
  Rcpp::NumericVector beta1_v(N);
  Rcpp::NumericVector tau1_v(N);
  Rcpp::NumericVector A0_v(N);
  Rcpp::NumericVector b0_v(N);
  Rcpp::NumericVector alpha0_v(N);
  Rcpp::NumericVector beta0_v(N);
  Rcpp::NumericVector tau0_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    A1_v(nv) = A1(A1_inc);
    b1_v(nv) = b1(b1_inc);
    alpha1_v(nv) = alpha1(alpha1_inc);
    beta1_v(nv) = beta1(beta1_inc);
    tau1_v(nv) = tau1(tau1_inc);
    A0_v(nv) = A0(A0_inc);
    b0_v(nv) = b0(b0_inc);
    alpha0_v(nv) = alpha0(alpha0_inc);
    beta0_v(nv) = beta0(beta0_inc);
    tau0_v(nv) = tau0(tau0_inc);

    A1_inc = A1_inc + 1;
    b1_inc = b1_inc + 1;
    alpha1_inc = alpha1_inc + 1;
    beta1_inc = beta1_inc + 1;
    tau1_inc = tau1_inc + 1;
    A0_inc = A0_inc + 1;
    b0_inc = b0_inc + 1;
    alpha0_inc = alpha0_inc + 1;
    beta0_inc = beta0_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_A1 == A1_inc) A1_inc = 0;
    if (N_b1 == b1_inc) b1_inc = 0;
    if (N_alpha1 == alpha1_inc) alpha1_inc = 0;
    if (N_beta1 == beta1_inc) beta1_inc = 0;
    if (N_tau1 == tau1_inc) tau1_inc = 0;
    if (N_A0 == A0_inc) A0_inc = 0;
    if (N_b0 == b0_inc) b0_inc = 0;
    if (N_alpha0 == alpha0_inc) alpha0_inc = 0;
    if (N_beta0 == beta0_inc) beta0_inc = 0;
    if (N_tau0 == tau0_inc) tau0_inc = 0;
  }

  // Generate draws
  for (int n = 0; n < N; n++) {

    // Simulate values from each accumulator
    double t1 = rlba_1acc_scl( A1_v(n),
                               b1_v(n),
                               alpha1_v(n),
                               beta1_v(n),
                               ver );
    double t0 = rlba_1acc_scl( A0_v(n),
                               b0_v(n),
                               alpha0_v(n),
                               beta0_v(n),
                               ver );

    if ( Rcpp::NumericVector::is_na(t1) |
         Rcpp::NumericVector::is_na(t0) ) {
      out(n,0) = NA_REAL;
      out(n,1) = NA_REAL;
    } else {
      // If the residual latency impacts the decision rule
      if ( rl == 1 ) {
        t1 = t1 + tau1_v(n);
        t0 = t0 + tau0_v(n);
      }

      double rt = -1.0;
      double ch = -1.0;

      // Determine the winning accumulator
      if (t1 < t0) {
        rt = t1;
        if ( rl == 0 ) rt = rt + tau1_v(n);
        ch = 1;
      }
      if (t1 > t0) {
        rt = t0;
        if ( rl == 0 ) rt = rt + tau0_v(n);
        ch = 0;
      }

      // Save the results
      out(n,0) = rt; out(n,1) = ch;
    }

  }

  // Set column names
  colnames(out) = Rcpp::CharacterVector::create("rt", "ch");

  return ( out );
}

// Lookup - 25
// A scalar version for the dlba function

double dlba_scl( double rt, double ch,
                 double A1, double b1,
                 double alpha1, double beta1, double tau1,
                 double A0, double b0,
                 double alpha0, double beta0, double tau0,
                 double rl, int ln, int ver ) {

  double dt; double pt;
  if ( rl == 1.0 ) {
    // If residual latency impacts decisions
    dt = rt - ch*tau1 - (1-ch)*tau0;
    pt = rt - ch*tau0 - (1-ch)*tau1;
  } else {
    // If residual latency is separate
    dt = rt - ch*tau1 - (1-ch)*tau0;
    pt = rt - ch*tau1 - (1-ch)*tau0;
  }

  double dA = ch*(A1) + (1-ch)*A0;
  double pA = ch*(A0) + (1-ch)*A1;
  double db = ch*(b1) + (1-ch)*b0;
  double pb = ch*(b0) + (1-ch)*b1;
  double dalpha = ch*(alpha1) + (1-ch)*alpha0;
  double palpha = ch*(alpha0) + (1-ch)*alpha1;
  double dbeta = ch*(beta1) + (1-ch)*beta0;
  double pbeta = ch*(beta0) + (1-ch)*beta1;

  double cdf = 1.0;
  double out = log(0.0);

  if ( (A1 > 0.0 ) & (A1 <= b1) &
       (A0 > 0.0 ) & (A0 <= b0) &
       (dt > 0.0) ) {
    if (pt > 0.0) cdf = 1.0 - plba_1acc_scl(pt,pA,pb,palpha,pbeta,ver);
    out = dlba_1acc_scl(dt,dA,db,dalpha,dbeta,ver,1) + log(cdf);
  }
  if (ln==0) out = exp( out );

  return( out );
}

// Lookup - 26
//' @rdname rlba
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector dlba ( Rcpp::NumericVector rt,
                           Rcpp::NumericVector ch,
                           Rcpp::NumericVector A1,
                           Rcpp::NumericVector b1,
                           Rcpp::NumericVector alpha1,
                           Rcpp::NumericVector beta1,
                           Rcpp::NumericVector tau1,
                           Rcpp::NumericVector A0,
                           Rcpp::NumericVector b0,
                           Rcpp::NumericVector alpha0,
                           Rcpp::NumericVector beta0,
                           Rcpp::NumericVector tau0,
                           double rl = 0.0, int ln = 0,
                           int ver = 0 ) {

  int N_rt = rt.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_A1 = A1.size(); // Number of parameters
  int N_b1 = b1.size();
  int N_alpha1 = alpha1.size();
  int N_beta1 = beta1.size();
  int N_tau1 = tau1.size();
  int N_A0 = A0.size();
  int N_b0 = b0.size();
  int N_alpha0 = alpha0.size();
  int N_beta0 = beta0.size();
  int N_tau0 = tau0.size();

  // Increment variables for loop
  int rt_inc = 0.0;
  int ch_inc = 0.0;
  int A1_inc = 0.0;
  int b1_inc = 0.0;
  int alpha1_inc = 0.0;
  int beta1_inc = 0.0;
  int tau1_inc = 0.0;
  int A0_inc = 0.0;
  int b0_inc = 0.0;
  int alpha0_inc = 0.0;
  int beta0_inc = 0.0;
  int tau0_inc = 0.0;

  // Determine the longest input vector
  int N = max(Rcpp::NumericVector::create(N_rt, N_ch, N_A1, N_b1,
                                          N_alpha1, N_beta1, N_tau1,
                                          N_A0, N_b0, N_alpha0,
                                          N_beta0, N_tau0 ));

  // Set output vector
  Rcpp::NumericVector out(N);

  // Create vectors for parameters
  Rcpp::NumericVector rt_v(N);
  Rcpp::NumericVector ch_v(N);
  Rcpp::NumericVector A1_v(N);
  Rcpp::NumericVector b1_v(N);
  Rcpp::NumericVector alpha1_v(N);
  Rcpp::NumericVector beta1_v(N);
  Rcpp::NumericVector tau1_v(N);
  Rcpp::NumericVector A0_v(N);
  Rcpp::NumericVector b0_v(N);
  Rcpp::NumericVector alpha0_v(N);
  Rcpp::NumericVector beta0_v(N);
  Rcpp::NumericVector tau0_v(N);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {
    rt_v(nv) = rt(rt_inc);
    ch_v(nv) = ch(ch_inc);
    A1_v(nv) = A1(A1_inc);
    b1_v(nv) = b1(b1_inc);
    alpha1_v(nv) = alpha1(alpha1_inc);
    beta1_v(nv) = beta1(beta1_inc);
    tau1_v(nv) = tau1(tau1_inc);
    A0_v(nv) = A0(A0_inc);
    b0_v(nv) = b0(b0_inc);
    alpha0_v(nv) = alpha0(alpha0_inc);
    beta0_v(nv) = beta0(beta0_inc);
    tau0_v(nv) = tau0(tau0_inc);

    rt_inc = rt_inc + 1;
    ch_inc = ch_inc + 1;
    A1_inc = A1_inc + 1;
    b1_inc = b1_inc + 1;
    alpha1_inc = alpha1_inc + 1;
    beta1_inc = beta1_inc + 1;
    tau1_inc = tau1_inc + 1;
    A0_inc = A0_inc + 1;
    b0_inc = b0_inc + 1;
    alpha0_inc = alpha0_inc + 1;
    beta0_inc = beta0_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_rt==rt_inc) rt_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_A1 == A1_inc) A1_inc = 0;
    if (N_b1 == b1_inc) b1_inc = 0;
    if (N_alpha1 == alpha1_inc) alpha1_inc = 0;
    if (N_beta1 == beta1_inc) beta1_inc = 0;
    if (N_tau1 == tau1_inc) tau1_inc = 0;
    if (N_A0 == A0_inc) A0_inc = 0;
    if (N_b0 == b0_inc) b0_inc = 0;
    if (N_alpha0 == alpha0_inc) alpha0_inc = 0;
    if (N_beta0 == beta0_inc) beta0_inc = 0;
    if (N_tau0 == tau0_inc) tau0_inc = 0;
  }

  // Calculate likelihood
  for (int n = 0; n < N; n++) {
    out(n) = dlba_scl(rt_v(n), ch_v(n), A1_v(n), b1_v(n),
        alpha1_v(n), beta1_v(n), tau1_v(n), A0_v(n), b0_v(n),
        alpha0_v(n), beta0_v(n), tau0_v(n), rl, ln, ver);
  }

  return ( out );
}

// Lookup - 27
// A variant of the scalar version of the density function to
// integrate over.

double int_dlba_scl( double x, void * params) {

  // Extract parameters
  std::vector<double> par = *(std::vector<double> *) params;

  // Initialize output
  double out = 0.0;

  // Calculate the density for the Wald race model
  int ver = par[12];
  out = dlba_scl( x, par[0],
                  par[1], par[2], par[3], par[4], par[5],
                  par[6], par[7], par[8], par[9], par[10],
                  par[11], 0, ver );

  return out;
}

// Lookup - 28
// A function that numerically integrates the density function in
// order to determine the distribution function.

double plba_scl( std::vector<double> par, double a,double b) {

  // Turn off GSL error handler
  gsl_set_error_handler_off ();

  double result = 0.0;

  if ( a >= 0.0 ) {

    // Allocate memory
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double error;

    gsl_function F;
    F.function = &int_dlba_scl;
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

// Lookup - 29
// RcppParallel worker function

struct plbaWorker : public Worker
{
  // Input matrix
  const RMatrix<double> input;

  // Destination matrix
  RVector<double> output;

  // initialize with source and destination
  plbaWorker(const Rcpp::NumericMatrix input,
                  Rcpp::NumericVector output)
    : input(input), output(output) {}

  // function call operator working for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {

    for(std::size_t j = begin; j < end; j++) {

      double cur_rt = input(j,0); // Extract RT

      std::vector<double> par(13); // Extract parameters
      for (int i = 1; i < 14; i++) { par[i-1] = input(j,i); }

      output[j] = plba_scl( par, 0.0, cur_rt );
    }
  }
};

// Lookup - 30
//' @rdname rlba
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector plba ( Rcpp::NumericVector rt,
                           Rcpp::NumericVector ch,
                           Rcpp::NumericVector A1,
                           Rcpp::NumericVector b1,
                           Rcpp::NumericVector alpha1,
                           Rcpp::NumericVector beta1,
                           Rcpp::NumericVector tau1,
                           Rcpp::NumericVector A0,
                           Rcpp::NumericVector b0,
                           Rcpp::NumericVector alpha0,
                           Rcpp::NumericVector beta0,
                           Rcpp::NumericVector tau0,
                           double rl = 0.0, double ver = 0.0,
                           int parYes = 1 ) {

  int N_rt = rt.size(); // Number of response times
  int N_ch = ch.size(); // Number of choices
  int N_A1 = A1.size(); // Number of parameters
  int N_b1 = b1.size();
  int N_alpha1 = alpha1.size();
  int N_beta1 = beta1.size();
  int N_tau1 = tau1.size();
  int N_A0 = A0.size();
  int N_b0 = b0.size();
  int N_alpha0 = alpha0.size();
  int N_beta0 = beta0.size();
  int N_tau0 = tau0.size();

  // Increment variables for loop
  int rt_inc = 0.0;
  int ch_inc = 0.0;
  int A1_inc = 0.0;
  int b1_inc = 0.0;
  int alpha1_inc = 0.0;
  int beta1_inc = 0.0;
  int tau1_inc = 0.0;
  int A0_inc = 0.0;
  int b0_inc = 0.0;
  int alpha0_inc = 0.0;
  int beta0_inc = 0.0;
  int tau0_inc = 0.0;

  // Determine the longest input vector
  int N = max(Rcpp::NumericVector::create(N_rt, N_ch, N_A1, N_b1,
                                          N_alpha1, N_beta1, N_tau1,
                                          N_A0, N_b0, N_alpha0,
                                          N_beta0, N_tau0 ));
  // Structure of matrix
  // prm(,0) = rt; prm(,1) = ch;
  // prm(,2) = A1; prm(,3) = b1;
  // prm(,4) = alpha1; prm(,5) = beta1;
  // prm(,6) = tau1;
  // prm(,7) = A0; prm(,8) = b0;
  // prm(,9) = alpha0; prm(,10) = beta0;
  // prm(,11) = tau0;
  // prm(,12) = rl; prm(,13) = ver;

  // Set output vector
  Rcpp::NumericVector output(N);

  // Set input matrix for parameters
  Rcpp::NumericMatrix input(N,14);

  // Loop through observations
  for (int nv = 0; nv < N; nv++) {

    input(nv,0) = rt(rt_inc);
    input(nv,1) = ch(ch_inc);
    input(nv,2) = A1(A1_inc);
    input(nv,3) = b1(b1_inc);
    input(nv,4) = alpha1(alpha1_inc);
    input(nv,5) = beta1(beta1_inc);
    input(nv,6) = tau1(tau1_inc);
    input(nv,7) = A0(A0_inc);
    input(nv,8) = b0(b0_inc);
    input(nv,9) = alpha0(alpha0_inc);
    input(nv,10) = beta0(beta0_inc);
    input(nv,11) = tau0(tau0_inc);
    input(nv,12) = rl;
    input(nv,13) = ver;

    rt_inc = rt_inc + 1;
    ch_inc = ch_inc + 1;
    A1_inc = A1_inc + 1;
    b1_inc = b1_inc + 1;
    alpha1_inc = alpha1_inc + 1;
    beta1_inc = beta1_inc + 1;
    tau1_inc = tau1_inc + 1;
    A0_inc = A0_inc + 1;
    b0_inc = b0_inc + 1;
    alpha0_inc = alpha0_inc + 1;
    beta0_inc = beta0_inc + 1;
    tau0_inc = tau0_inc + 1;
    if (N_rt==rt_inc) rt_inc = 0;
    if (N_ch==ch_inc) ch_inc = 0;
    if (N_A1 == A1_inc) A1_inc = 0;
    if (N_b1 == b1_inc) b1_inc = 0;
    if (N_alpha1 == alpha1_inc) alpha1_inc = 0;
    if (N_beta1 == beta1_inc) beta1_inc = 0;
    if (N_tau1 == tau1_inc) tau1_inc = 0;
    if (N_A0 == A0_inc) A0_inc = 0;
    if (N_b0 == b0_inc) b0_inc = 0;
    if (N_alpha0 == alpha0_inc) alpha0_inc = 0;
    if (N_beta0 == beta0_inc) beta0_inc = 0;
    if (N_tau0 == tau0_inc) tau0_inc = 0;
  }

  // Calculate likelihood
  if (parYes == 0) {

    for (int j = 0; j < N; j++) {

      double cur_rt = input(j,0);
      std::vector<double> par(13);
      for (int i = 1; i < 14; i++) { par[i-1] = input(j,i); }

      output(j) = plba_scl( par, 0.0, cur_rt );
    }
  } else {

    // Function call operator that works for the specified
    // range (begin/end)
    plbaWorker mt(input, output);

    // Call parallelFor to do the work
    parallelFor(0, N, mt);
  }

  return ( output );
}
