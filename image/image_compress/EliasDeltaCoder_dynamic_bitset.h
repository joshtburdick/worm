#pragma once

#include <stdint.h>

#include <boost/dynamic_bitset.hpp>

using namespace boost;

/** Elias delta coding, based on
  http://en.wikipedia.org/wiki/Elias_delta_coding
  If I understood the C++ better, this could "write"
  to a stream; not bothering with this for now. */
struct EliasDeltaCoder {

  /** Output is written here. */
  dynamic_bitset<uint64_t> bits;

  EliasDeltaCoder() {
  }

  /** Postpends a number onto the end. */
  void writeNumber(int num) {
    int len = 0;
    int lengthOfLen = 0;

    // calculate 1+floor(log2(num))
    for (int temp = num; temp > 0; temp >>= 1)
      len++;

    // calculate floor(log2(len))
    for (int temp = len; temp > 1; temp >>= 1)
      lengthOfLen++;

    for (int i = lengthOfLen; i > 0; i--)
      bits.push_back(0);
    for (int i = lengthOfLen; i >= 0; i--)
      bits.push_back((len >> i) & 1);
    for (int i = len-2; i >= 0; i--)
      bits.push_back((num >> i) & 1);
  }
};


