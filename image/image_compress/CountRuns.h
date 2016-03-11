#pragma once

#include <stdint.h>
#include <map>

#include "RunCounter.h"

/** XXX this is probably available somewhere. */
#define bitTest(a, i) ((a & (0x1 << (i))) != 0)

using namespace std;

/** Counts runs of one bit (original version.) */
template <class T>
map<uint32_t, uint32_t> countRuns1(vector<T> &x, int bitIndex) {
  map<uint32_t, uint32_t> r;

  bool bit = 0;
  int runLength = 0;
  for(uint32_t i=0; i<x.size(); i++) {
    bool bit1 = (x[i] & (0x1 << bitIndex)) != 0;
    if (bit == bit1)
      runLength++;
    else {
      r[runLength]++;
      runLength = 1;
      bit = bit1;
    }
  }
  r[runLength]++;

  return r;
}

/** Counts runs of one bit (new version.) */
template <class T>
map<uint32_t, uint32_t> countRuns(vector<T> &x, int bitIndex) {

  RunCounter rc;

  for(uint32_t i=0; i<x.size(); i++)
    rc.countBit( bitTest(x[i], bitIndex) );

  return rc.getCounts();
}

/** Counts runs of one bit (wavelet-transformed.) 
  XXX I don't think this works. */
template <class T>
map<uint32_t, uint32_t> countRunsWavelet(vector<T> &x, int bitIndex, int maxBits) {

  // if this is the leftmost bit, just count the runs as before
  if (bitIndex == (maxBits-1))
    return countRuns(x, bitIndex);

  RunCounter rc;

  // first, count bits in which the previous bit is 0  
  for(uint32_t i=0; i<x.size(); i++)
    if ( ! bitTest(x[i], bitIndex+1) )
      rc.countBit( bitTest(x[i], bitIndex) );

  // ... then when the previous bit is 1
  for(uint32_t i=0; i<x.size(); i++)
    if (bitTest(x[i], bitIndex+1))
      rc.countBit( bitTest(x[i], bitIndex) );

  return rc.getCounts();
}




