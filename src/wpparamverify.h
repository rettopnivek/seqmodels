// Include guard to protect against multiple definitions error
#ifndef __WPPARAMVERIFY__
#define __WPPARAMVERIFY__

#include <Rcpp.h>

/*
Purpose:
Defines a function to verify valid parameter inputs for
the wiener process that can be included in multiple files.

Index:
Lookup - 01:  wp_param_verify
*/

// Lookup - 01

inline bool wp_param_verify( std::vector<double> prm,
                             std::vector<int> index ) {
  /*
  Purpose:
  Function to verify whether inputs are admissable for the
  two-boundary wiener process.
  Arguments:
  prm   - A vector of parameters
  index - The type of parameter being inputted, where...
  0 = skip
  1 = response time
  2 = choice
  3 = alpha
  4 = theta
  5 = xi
  6 = tau
  7 = sigma
  8 = probability
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
  // Variable for whether a response time is present
  double rt_yes = 0.0;

  for ( int i = 0; i < sz; i++ ) {

    if ( index[i] > 0 ) n_pass += 1;

    // Check whether inputted response time is a real number
    // and greater than 0
    if ( index[i] == 1 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) {
        chk += 1;
        rt_yes += prm[i];
      }
    }

    // Check whether inputted choice is a real number and
    // whether it is 0 or 1
    if ( index[i] == 2 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( ( prm[i] == 0.0 ) || ( prm[i] == 1.0 ) ) ) chk += 1;
    }

    // Check whether boundary separation is a real number and
    // greater than 0
    if ( index[i] == 3 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether starting point is a real number and whether
    // it is bounded between 0 and 1
    if ( index[i] == 4 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) && ( prm[i] < 1.0 ) ) chk += 1;
    }

    // Check whether drift rate is a real number
    if ( index[i] == 5 ) {
      if ( !Rcpp::NumericVector::is_na( prm[i] ) ) chk += 1;
    }

    // Check whether residual latency is a real number and
    // whether it is greater than or equal to 0 and whether
    // it is less than the inputted response time
    if ( index[i] == 6 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] >= 0.0 ) ) {
        // Check whether residual latency is less than the
        // given response time
        if ( rt_yes > 0.0 ) {
          if ( rt_yes > prm[i] ) chk += 1;
        } else {
          chk += 1;
        }
      }
    }

    // Check whether coefficient of drift is a real number and
    // whether it is greater than 0
    if ( index[i] == 7 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( prm[i] > 0.0 ) ) chk += 1;
    }

    // Check whether probability is a real number and whether it
    // is appropriately bounded
    if ( index[i] == 8 ) {
      if ( ( !Rcpp::NumericVector::is_na( prm[i] ) ) &&
           ( (prm[i] >= 0.0) | (prm[i] <= 1.0) ) ) chk += 1;
    }

  }

  // Test how many checks were passed
  out = chk == n_pass;

  return( out );
}

#endif // End include guard
