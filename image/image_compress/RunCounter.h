#pragma once

/** Utility to count the lengths of runs of the same bit. */

#include <stdint.h>
#include <map>

using namespace std;

class RunCounter {

  /* The bit we just saw. */
  bool previousBit;

  /* How many times we've seen it. */
  int runLength;  

  /* How many runs of a given length we've seen. */
  map<uint32_t, uint32_t> lengthHistogram;

public:
  /* Constructor. */
  RunCounter() {
    previousBit = false;
    runLength = 0;
  }

  /** Counts one bit. */
  void countBit(bool bit) {
    if (bit == previousBit)
      runLength++;
    else {
      lengthHistogram[runLength]++;
      runLength = 1;
      previousBit = bit;
    }
  }

  /** Gets the bit counts (this includes counting the last run of numbers.) */
  map<uint32_t, uint32_t> getCounts() {
    lengthHistogram[runLength]++;
    return lengthHistogram;
  }

};

