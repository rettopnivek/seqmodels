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
rtdist packages (Forthcoming)
Singmann, H., Brown, S., Gretton, M., Heathcote, A., Voss, A., Voss, J.,
  & Terry, A. (2016). rtdists: Response Time Distributions. R package
version 0.6-6.

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

### TO DO ###
- Add density functions
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
