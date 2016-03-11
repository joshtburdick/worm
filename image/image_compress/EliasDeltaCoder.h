#pragma once

#include <stdint.h>

#include <vector>
// #include <boost/dynamic_bitset.hpp>

/** Elias delta coding, based on
  http://en.wikipedia.org/wiki/Elias_delta_coding
  If I understood the C++ better, this could "write"
  to a stream; not bothering with this for now. */
struct EliasDeltaCoder {

  /** Output is written here. */
  vector<uint64_t> bits;

  EliasDeltaCoder() {
    a = 0L;
    b = 0;
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
      writeBit(0);
    for (int i = lengthOfLen; i >= 0; i--)
      writeBit((len >> i) & 1);
    for (int i = len-2; i >= 0; i--)
      writeBit((num >> i) & 1);
  }

  /** Gets the vector of words. */
  vector<uint64_t> & getBits() {
    bits.push_back(a);
    return bits;
  }

private:

  /** Current word being written. */
  uint64_t a;

  /** Index of bit being written in this word. */
  int b;

  /** Tacks on one bit. */
  void writeBit(bool x) {
    if (x) a |= 1L << b;
    b++;
    if (b == 64) {
      bits.push_back(a);
      a = 0L;
      b = 0;
    }      
  }

};

