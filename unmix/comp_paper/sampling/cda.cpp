/** C++ implementation of hit-and-run sampling. */

#include <Rcpp.h>
// #include <RcppArmadillo.h>  // this includes Rcpp.h

#include <iostream>

using namespace std;
using namespace Rcpp;

/** Picks a random amount to go in some direction.
  Args:
    x - current position
    d - a direction vector (needn't be normalized)
  Returns: a random scalar r, such that all elements
    of x + ra >= 0. (Note that r can be negative.) */
double jumpDist(NumericVector x, NumericVector dir) {
  double lo = -1, hi = 1;

  // loop through the elements
  for(int i=0; i<x.size(); i++) {

    double d = dir[i];
    // if this entry of the direction is nonzero...
    if (1) {
      double r = x[i] / dir[i];
      if (r > 0 && -r >= lo)
        lo = -r;
      if (r < 0 && -r <= hi)
        hi = -r;
    }
  }

  // finally, pick a number between a and b
  return Rf_runif(lo, hi);
}

/** Does hit-and-run sampling.
    Args:
    A, b - matrix and vector giving the constraint A x = b .
    x0 - initial solution (which should satisfy A x0 = b)
    Z - matrix whose rows span the nullspace of A
    numSamples - number of samples to return
    thinning - interval at which to save samples
    (used to reduce autocorrelation. FIXME currently does nothing)
    Returns: numSamples samples. */
// [[Rcpp::export]]
NumericMatrix cdaCpp(NumericMatrix A,
		     NumericVector b,
		     NumericVector x0,
		     NumericMatrix Z,
		     int numSamples,
		     int thinning) {
  RNGScope scope;

  int n = A.ncol();   // number of cells

  NumericVector x(n), // current location
    dir(n),           // direction of jump
    g(Z.nrow());       // accumulated jumps in nullspace
  x = x0;
  NumericMatrix samples(numSamples, n);

  for(int iter = 0; iter < numSamples; iter++) {

    // which random direction in which to jump
    int z = random() % Z.nrow();
    NumericVector dir = Z(z, _);

    // amount to go in that direction
    double a = jumpDist(x0, dir);

    // update current location (in both coordinate systems)
    g[z] += a;
    x += a * dir;

    // FIXME: renormalize, to avoid rounding error?


    samples(iter, _) = x;
  }

  return wrap(samples);
}

