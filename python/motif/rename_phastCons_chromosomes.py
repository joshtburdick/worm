#!/usr/bin/python3
# Renames phastCons bigWig chromosome names to what WS220 uses.

import os
import re
import subprocess

dataDir = "/murrlab/seq/igv/conservation/"

inputFile = dataDir + "ce10.phastCons7way.bw"
outputFile = "/var/tmp/WS220.phastCons7way.bw"

tmpFile1 = "/var/tmp/conservation1.wig"
tmpFile2 = "/var/tmp/conservation2.wig"
chromSizes = "/murrlab/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.sizes"

# then, compute conservation at those locations
subprocess.call(["bigWigToWig", inputFile, tmpFile1])

# rename chromosomes
inF = open(tmpFile1, "r")
outF = open(tmpFile2, "w")
for s in inF:
  # for speed, only alter lines that are at least this long
  if (len(s) >= 7):
    s = re.sub(" chrom=chrM", " chrom=MtDNA", s)
    s = re.sub(" chrom=chr", " chrom=", s)
  outF.write(s)
inF.close()
outF.close()

# convert back to bigWig
subprocess.call(["wigToBigWig", tmpFile2, chromSizes, outputFile])

# clean up
os.remove(tmpFile1)
os.remove(tmpFile2)

