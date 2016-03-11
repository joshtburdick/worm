/** Counts runs of bits at a given position in an array. */

#include <iostream>
#include <fstream>
#include <vector>

// #include <lfs.h>
#include <sais.h>
#include "CountRuns.h"
#include "WaveletTransform.h"

/** Reads a file into a vector. Modified from
http://www.cplusplus.com/reference/iostream/istream/read/ */
void readStringIntoBuffer (string inFileName, vector<uint16_t> &buffer) {
  ifstream is;
  is.open(inFileName.c_str(), ios::binary);

  // get length of file:
  is.seekg (0, ios::end);
  int length = is.tellg();
  is.seekg (0, ios::beg);

  // resize buffer: XXX this may provide too much room
  buffer.resize( length / 2 + 1 );

  // read data as a block:
  is.read((char *) &(buffer.front()), length);
  is.close();
}

int main(int argc, char ** argv) {

  if (argc != 2) {
    cerr << "usage: countRuns FILE" << endl;
    return 1;
  }

  // read file into a vector
  string fileName(argv[1]);
  vector<uint16_t> buffer;
  cout << "reading from " << fileName << endl;
  readStringIntoBuffer(fileName, buffer);
  
  // BWT-transform it
  sais_u16_bwt(&(buffer.front()), &(buffer.front()),
    NULL, buffer.size(), 65536);

  // print counts for each bit
  cout << "bitIndex\tlength\tcount" << endl;
  for(int bitIndex = 0; bitIndex < 32; bitIndex++) {

    // get counts
    map<uint32_t, uint32_t> count =
      countRunsWavelet<uint16_t>(buffer, bitIndex, 16);

    // print them
    for(map<uint32_t, uint32_t>::iterator p =
        count.begin(); p != count.end(); p++) {
      cout << bitIndex << "\t" <<
        p->first << "\t" << p->second << endl;
    }
  }

  return 0;
}

