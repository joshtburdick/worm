#pragma once

#include <stdint.h>
#include <map>

// #include "RunDecompressorElias.h"


using namespace std;

/** Compresses a file using the BWT and (bitwise)
    Elias gamma coding. */
struct BWTEliasCompressor {

  BWTEliasCompressor(istream &_is, ostream &_os) : is(_is), os(_os) {
  }

  /** Writes out a compressed representation. */
  void writeCompressed() {

    readBlock();

    // BWT-transform it
    pidx = sais_u8_bwt(&(buf.front()), &(buf.front()),
      NULL, buf.size(), 256);

    cout << pidx << endl;

    // count runs of each bit
    rc.countRuns(buf);

    // write out compressed representation
    writeBlockCompressed();


  }

private:

  /** where to read from */
  istream &is;

  /** where to write to */
  ostream &os;

  /** the buffer */
  vector<uint8_t> buf;

  /** the "start index" */
  uint64_t pidx;

  /** compresses runs */
  RunCompressorElias<uint8_t> rc;

  /** Reads a block into buf. */
  void readBlock() {

    // (XXX ignoring all sorts of efficiency)
    // (this seems to actually be a hot spot)
    while (!is.eof())
      buf.push_back(cin.get());
  }

  /** Writes a block compressed. */
  void writeBlockCompressed() {

    // first, write a magic number
    os << "#BWTElias compressed file\n";

    // then, the length of the file (in bytes)
    writeWord(buf.size());

    // then, the pidx (index of the start of the string)
    writeWord(pidx);

    // then, the individual bit planes
    for(int i=0; i<8; i++) {
      vector<uint64_t> & b = rc.c[i].bits;
      writeWord(b.size());
      for(size_t j=0; j<b.size(); j++)
        writeWord(b[j]);
    }
  }

  /** Utility to write a uint64_t to the output stream. */
  void writeWord(uint64_t a) {
    // FIXME check byte order
    os.write((const char *) &a, sizeof(uint64_t));
  }

};

