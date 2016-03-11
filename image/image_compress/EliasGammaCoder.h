#pragma once

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include <vector>

using namespace std;

/** Elias gamma coding, based on
  http://en.wikipedia.org/wiki/Elias_gamma_coding */
struct EliasGammaCoder {

  /** Output is written here. */
  vector<uint64_t> bits;

  /** The current word being written. */
  uint64_t b;

  /** Number of bits available in the current word. */
  int numBitsLeft;

  EliasGammaCoder() {
    clear();
  }

  /** Clears this. */
  void clear() {
    bits.resize(0);
    b = 0L;
    numBitsLeft = 64;
  }

  /** Appends a number onto the end. */
  void writeNumber(uint64_t a) {
    assert( a > 0 );
// cout << "writing " << a << endl;
    // compute the number of bits required
    int numBits = msb( a );
// cout << " " << numBits << " bits" << endl;

    // encode that by prepending zeros
    writeZeros( numBits );

    // write out the actual number
    writeBits(a, numBits);
// cout << " bits left = " << numBitsLeft << endl;
  }

  /** Must be called after you finish writing bits,
    to "flush the cache" of the last word. */
  void writeLastWord() {
    if (numBitsLeft < 64) {
      bits.push_back(b);
    }
  }

private:

  /** Postpends some number of 0's. */
  void writeZeros(int n) {

    // if there's no room for the 0's...
    if (n > numBitsLeft) {

      // tack on a new empty word
      bits.push_back(b);
      b = 0L;
      n -= numBitsLeft;
      numBitsLeft = 64;
    }

    // at this point, the remaining 0's should
    // fit in this word    
    numBitsLeft -= n;
  }

  /** Appends some number of bits (the LSBs of x.)
    Args:
      a - number with bits to be appended
      n - number of bits to append (<= 64)
    Side effects: appends that many of the LSBs of x. */
  void writeBits(uint64_t a, int n) {
    // if the bits won't all fit in this word
    if (n > numBitsLeft) {
// cout << "writeBits overflow: n = " << n << " numBitsLeft = " << numBitsLeft << endl;
      // write the bits that will fit
//      int numToWrite = n - numBitsLeft;
      b |= (a >> (n - numBitsLeft));

      // clear the bits that we just wrote
      a &= ((1L << n) - 1L);

      // update count of bits left to write
      n -= numBitsLeft;

      // start a new empty word
      bits.push_back(b);
      b = 0L;
      numBitsLeft = 64;
    }

    // the remaining bits should now fit in the last word
    b |= a << (numBitsLeft - n);
    numBitsLeft -= n;
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

