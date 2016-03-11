#pragma once

#include <stdint.h>

#include <boost/dynamic_bitset.hpp>

using namespace boost;

/** Elias delta decoding, based on
  http://en.wikipedia.org/wiki/Elias_delta_coding */
class EliasDeltaDecoder {

  /** Input is read from here. */
  dynamic_bitset<uint64_t> bits;

  /** Which bit we're at (this could be replaced with some
   * sort of iterator, but for now, not bothering.) */
  uint64_t b;

public:
  EliasDeltaDecoder(dynamic_bitset<uint64_t> & _bits) {
    bits = _bits;
    b = 0;
  }

  /** Gets one number from the input stream.
    Returns: a number from the end of the stream. Note that
    results are undefined if you read past the end of the stream. */
  uint64_t readNumber() {
    int num = 1, len = 1, lengthOfLen = 0;

    // read the length-of-the-length
    // potentially dangerous with malformed files.
    while (!bits[b++])
      lengthOfLen++;

    // read the length
    for (int i = 0; i < lengthOfLen; i++) {
      len <<= 1;
      if (bits[b++])
        len |= 1;
    }

    // read the actual bits
    for (int i = 0; i < len-1; i++) {
      num <<= 1;
      if (bits[b++])
        num |= 1;
     }

     return num;
  }
};

