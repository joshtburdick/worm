#pragma once

#include <vector>

using namespace std;

/** Functions for doing (and undoing) a wavelet transform. Pretty inefficient. */
template <class T> struct WaveletTransform {

  /** The data to transform. This is transformed in-place. */
  vector<T> & x;

  /** Buffers for storing anything prefixed by 0 or 1, respectively. */
  vector<T> x0;
  vector<T> x1;
  
  /** Constructor. */
  WaveletTransform(vector<T> & _x) : x(_x) {
  }

  /** Does the transform. */
  void transform(int numBits) {
    // starts at 2, since we leave the first bit unchanged
    for(int b = numBits - 2; b >= 0; b--)
      transformBit(b);
  }

private:
  /** Transforms one bit. */
  void transformBit(int bit) {
    x0.clear();
    x1.clear();

    // useful bitwise masks
    T a = 1L << bit;
    T b = a - 1L;

    // split LSBs into two groups
    for(uint64_t i = 0; i < x.size(); i++) {
      if (a & x[i])
        x1.push_back(b & x[i]);
      else
        x0.push_back(b & x[i]);
    }

    // copy the LSBs back in that order
    uint64_t n0 = x0.size();
    for(uint64_t i = 0; i < n0; i++)
      x[i] = (b ^ x[i]) | x0[i];
    for(uint64_t i = 0; i < x1.size(); i++)
      x[i+n0] = (b ^ x[i+n0]) | x1[i];
  }
};

