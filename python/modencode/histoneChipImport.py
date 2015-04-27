#!/usr/bin/python
# 
# Converts histone ChIP DNA-seq .wig.gz files to .bw, including parsing
# filenames.

import os
import re
import subprocess

# FIXME: should these be command-line args?
wigGzDir = "/media/jburdick/disk2/data/ftp/data.modencode.org/C.elegans/Histone-Modification/ChIP-seq/coverage-graph_wiggle/"
bwDir = "/media/jburdick/disk2/histone_chip_seq_bw_new/"

sizesFile = "/murrlab/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.sizes"

# Uncompresses a wig.gz file, and alters it slightly, so that hopefully
# wigToBigWig reads it OK. This includes:
# - renaming chromosomes from, e.g., chrI => I
# - removing overlaps, which otherwise cause wigToBigWig to fail
#   (arguably the thing to do).
# Args:
#   wigGzFile - the wig.gz file to read
#   outputBedGraph - where to write output to
# Side effects: writes a bedGraph file 
def removeOverlaps(wigGzFile, outputBedGraph):
  # files to read and write
  wig = subprocess.check_output(["gunzip", "-c", wigGzFile])
  f = open(outputBedGraph, "w")

  # most recent position that we've seen, and written out
  lastPos = -1000

  for s in wig.splitlines():
    if (len(s) >= 15):
      # rename chromosomes, and reset the "last window seen" counter
      if ("chrom" in s):
#        s = re.sub(" chrom=chrM", " chrom=MtDNA", s)
        s = re.sub(" chrom=chr", " chrom=", s)
        lastPos = -1000
      f.write(s + "\n")
    else:
      a = s.split("\t")
      if len(a) == 2:
        pos = int(a[0])
        # FIXME distance shouldn't be hardwired, ideally
        if (pos - lastPos >= 10):
          f.write(a[0] + "\t" + a[1] + "\n")
          lastPos = pos
  f.close()

# Converts one .wig.gz file.
def wigGzConvert(wigGzFile, bwFile):
  fTmp = bwDir + '/tmp.bedGraph'
#  os.system("gunzip -c \"" + wigGzFile + "\" | sed -e 's/chrom=chr/chrom=/g' > " + fTmp)
  removeOverlaps(wigGzFile, fTmp)
  r = subprocess.call(["wigToBigWig", "-clip", fTmp, sizesFile, bwFile])
  if (r==0):
    os.system("rm " + fTmp)
  else:
    print("wigToBigWig failed for " + wigGzFile)

for f in os.listdir(wigGzDir):
  print(f)
#  if re.search(r'.combined.', f):
#    fNew = re.sub(r'.wig.gz$', '.bw', f)
#    print("converting " + f + " to " + fNew)
  m = re.match(r'^([^:]+):Developmental-Stage=([^#]+)#.*:ChIP-seq:Rep-(.):(ChIP|input):Cele_WS220:modENCODE_\d+:(.*).wig.gz$', f)
  if m:
    tf = m.group(1)
    stage = m.group(2)
    replicate = m.group(3)
    chipOrInput = m.group(4)
    suffix = m.group(5)
    if True:
      print(tf + " " + stage + " " + replicate + " " + chipOrInput + " " + suffix)
      os.system("mkdir -p " + bwDir + "/" + stage)
      outputFile = bwDir + "/" + stage + "/" + tf + "_rep" + replicate + "_" + chipOrInput + "_" + suffix + ".bw"
      wigGzConvert(wigGzDir + "/" + f, outputFile)

