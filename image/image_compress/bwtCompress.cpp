
#include <stdlib.h>

#include <istream>
#include <fstream>
#include <iostream>
#include <vector>

// #include <lfs.h>
#include <sais.h>
#include "CountRuns.h"

// these are deprecated
#include "EliasDeltaCoder.h"
#include "EliasDeltaDecoder.h"

#include "BWTEliasCompressor.h"

#include "EliasGammaDecoder.h"

using namespace std;

/** Quick test of the Elias coder/decoder. */
void eliasCodecTest(int n, int maxNumber) {
  vector<uint64_t> x;
  EliasGammaCoder c;

  cout << "adding numbers:" << endl;
  for(int i=0; i<n; i++) {
    uint64_t a = (uint64_t) random() % maxNumber + 1;
    x.push_back(a);
    c.writeNumber(a);
    cout << a << " ";
  }
  cout << endl;
  c.writeLastWord();

cout << "c.bits.size = " << c.bits.size() << endl;

  cout << "decompressing:" << endl;
  EliasGammaDecoder dc(c.bits);
  for(int i=0; i<n; i++) {
    uint64_t a = dc.readNumber();
    cout << a << " ";
  }
  cout << endl;

}

/** Reads a file into a vector. Modified from
  http://www.cplusplus.com/reference/iostream/istream/read/
  Possibly not used. */
void readStringIntoBuffer (string inFileName, vector<uint8_t> &buffer) {
  ifstream is;
  is.open(inFileName.c_str(), ios::binary);

  // get length of file:
  is.seekg (0, ios::end);
  int length = is.tellg();
  is.seekg (0, ios::beg);

  // resize buffer: XXX this may provide too much room
  buffer.resize( length );

  // read data as a block:
  is.read((char *) &(buffer.front()), length);
  is.close();
}

/** Compresses a file. */
void compress(istream &is, ostream &os) {
  BWTEliasCompressor comp(is, os);
  comp.writeCompressed();
}

int main(int argc, char ** argv) {

  // for now, just testing
//  eliasCodecTest(100,10000);

  compress(cin, cout);

  return 0;
}

