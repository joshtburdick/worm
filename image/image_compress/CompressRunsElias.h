#pragma once

#include <stdint.h>
#include <map>

#include "EliasDeltaCoder.h"

using namespace std;

/** Compresses runs of a particular bit, using Elias delta encoding. */
template <class T>

  /** What we're counting bits of. */
  vector<T> x;

  

  /** Where to write numbers to. */
  EliasDeltaCoder c;

  /** Which bit we will count. */
  int bitIndex;




  void countRuns(vector<T> &x, int bitIndex) {

  // Where to write numbers to.
  EliasDeltaCoder c;

  // The bit we're currently counting. We initialize this to whatever
  // the first bit is (since the coder can't encode a 0.) The first
  // bit will have to be stored separately.
  bool bit = (x[i] & (0x1 << bitIndex)) != 0;

  // the count of that bit
  int runLength = 0;

  // loop through the entries of x
  for(uint32_t i=0; i<x.size(); i++) {
    bool bit1 = (x[i] & (0x1 << bitIndex)) != 0;
    if (bit == bit1)
      runLength++;
    else {
      c.writeNumber(runLength);
      runLength = 1;
      bit = bit1;
    }
  }
  c.writeNumber(runLength);

  return r;
}


