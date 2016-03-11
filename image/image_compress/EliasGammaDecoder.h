#pragma once

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include <vector>

using namespace std;

/** Elias gamma decoding, based on
  http://en.wikipedia.org/wiki/Elias_gamma_coding */
struct EliasGammaDecoder {

  /** Vector of words we're iterating through.
    (Slower than using pointers, but arguably safer.) */
  vector<uint64_t> & bits;

  /** Index of the current bit. */
  int i;

  /** The current word being read. As we read in bits,
    we clear them in this. */
  uint64_t b;

  /** Number of bits left to read in the current word. */
  int numBitsLeft;

  EliasGammaDecoder(vector<uint64_t> & _bits) : bits(_bits) {
    i = 0;
    b = bits[0];
    numBitsLeft = 64;
  }

  /** Reads in one number from the stream.
    Note that the behaviour of this is undefined if
    you read past the end of the stream. */
  uint64_t readNumber() {

    // count number of leading 0's
    int numBits = readZeros();
// cout << "numBits = " << numBits << endl;

    // get that many bits
    return readBits(numBits);
  }

private:

  /** Counts the number of leading zeros. */
  uint64_t readZeros() {
    int numZeros = 0;

    // first, if there are no bits set in the current word
    if (b == 0L) {
      numZeros = numBitsLeft;  // count these zeros
      b = bits[ ++i ];         // read the next word
      numBitsLeft = 64;        // and reset pointer to next bit
    }

    // at this point, there should be at least one bit set
    int a = msb(b);
    assert( a > 0 );
// cout << "numBitsLeft = " << numBitsLeft << ", a = " << a << endl;
    // add these zeros to the count
    numZeros += numBitsLeft - a;

    // and advance the index of the next bit
    numBitsLeft = a;

    return numZeros;
  }

  /** Reads some number of bits (most-significant bits
    first; starting from the "leftmost" bit in b which
    hasn't been read already.
    Args:
      n - number of bits to read (1 <= n <= 64)
    Side effects: advances current pointer n bits
    Returns: n bits */
  uint64_t readBits(int n) {
// cout << "in readBits: numBitsLeft = " << numBitsLeft << ", n = " << n << endl;   
    uint64_t a = 0L;

    // will we need all the bits in b ?
    if (n >= numBitsLeft) {

      // extract those bits (from LSBs of b)
      a = b << (n - numBitsLeft);

      // advance counters
      n -= numBitsLeft;
      numBitsLeft = 64;
      
      // read the next word
      b = bits[ ++i ];
    }

    // if we need to read in any more bits
    if (n > 0) {

      // at this point, we should only be reading bits
      // from this word
      a |= b >> (numBitsLeft - n);

      numBitsLeft -= n;

      // mask out bits we just read from b
      b &= (1L << numBitsLeft) - 1L;
    }

    return a;
  }

  /** Hopefully there's a faster version of this... */
  int msb(uint64_t n) {
    for(int i=63; i>=0; i--) {
      if (n & (1L << i))
        return i+1;
    }

    return 0;
  }

};

