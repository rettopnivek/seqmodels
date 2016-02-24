#include <Rcpp.h>
#include <vector>

/*
Purpose:
Miscellaneous Rcpp code for use in other functions.

Index:
Lookup - 01: minMax
Lookup - 02: linInterp

*/

// Lookup - 01
// A function that takes a vector and a comparison value, and determines the
// position of the first value in the vector that is either smaller or larger
// than the comparison value.

int minMax( double comp, std::vector<double> x, int Direction ) {

  // Determine which values are larger or smaller depending on the variable
  // direction
  int sz = x.size();
  std::vector<int> chk( sz );
  for (int i = 0; i < sz; i++) {
    chk[i] = 0;
    if ( (Direction == 1) && (x[i] > comp) ) {
      chk[i] = 1;
    }
    if ( (Direction == 0) && (x[i] < comp) ) {
      chk[i] = 1;
    }
  }

  // Find the first value that is smaller/larger
  int out;
  int stp = 0;
  if ( Direction == 1 ) {
    out = 0;
    int inc = sz - 1;
    while ( (stp==0) & (inc >= 0) ) {
      if (chk[inc] == 0) { out = inc+1; stp = 1; }
      inc = inc - 1;
    }
    if ( out == sz ) out = out - 1;
  } else {
    out = sz - 1;
    int inc = 0;
    while ( (stp==0) & (inc < sz) ) {
      if (chk[inc] == 0) { out = inc-1; stp = 1; }
      inc = inc + 1;
    }
    if ( out == -1 ) out = 0;
  }

  return( out );
}


// Lookup - 02
// A function that carries out a linear interpolation of an x-value given a set
// of cartesion coordinates ( x and y values ) and a new y-value

double linInterp( double yN, double y0, double y1,
                  double x0, double x1 ) {

  double b1 = ( y1 - y0 ) / ( x1 - x0 ); // Slope
  double b0 = y1 - b1*x1; // Intercept
  double num = yN - b0;
  double xN = ( num )/b1;

  return( xN );
}
