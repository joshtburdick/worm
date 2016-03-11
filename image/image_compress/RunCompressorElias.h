#pragma once

#include <stdint.h>
#include <map>

#include "EliasDeltaCoder.h"

// possibly using this
#include "EliasGammaCoder.h"


using namespace std;

/** Compresses runs of a particular bit, using Elias gamma encoding. */
template <class T>
struct RunCompressorElias {

  /** Where to write numbers to. */
  vector<EliasDeltaCoder> c;

  RunCompressorElias() {
    c.resize(8);
  }

  /** Counts runs of all bits. After calling this,
    c[bitIndex] should contain the compressed bit
    representations of each bit slice. */
  void countRuns(vector<uint8_t> & x) {
    for(int i=0; i<8; i++)
      countRuns(x, i);
  }

  /** Counts runs of a particular bit. */
  void countRuns(vector<uint8_t> & x, int bitIndex) {

    // The bit we're currently counting. We pretend that
    // the bitstring starts with a 0, since the Elias
    // code can't code a 0.
    bool bit = 0;

    // the count of that bit
    int runLength = 1;

    // loop through the entries of x
    for(uint32_t i=0; i<x.size(); i++) {
      bool bit1 = (x[i] & (0x1 << bitIndex)) != 0;
      if (bit == bit1)
        runLength++;
      else {
        c[bitIndex].writeNumber(runLength);
        runLength = 1;
        bit = bit1;
      }
    }

    // write out the last chunk
    c[bitIndex].writeNumber(runLength);
  }

};

