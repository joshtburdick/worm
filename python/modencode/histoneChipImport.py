#!/usr/bin/python
# 
# Converts histone ChIP DNA-seq .wig.gz files to .bw, including parsing
# filenames.

import os
import re

# FIXME: should these be command-line args?
wigGzDir = "/media/jburdick/disk2/data/ftp/data.modencode.org/C.elegans/Histone-Modification/ChIP-seq/coverage-graph_wiggle/"
bwDir = "/media/jburdick/disk2/histone_chip_seq_bw/"

sizesFile = "/murrlab/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.sizes"

# Converts one .wig.gz file.
def wigGzConvert(wigGzFile, bwFile):
  fTmp = re.sub('.bw$', '.bedGraph', bwFile)
  os.system("gunzip -c " + wigGzFile + " | sed -e 's/chrom=chr/chrom=/g' > " + fTmp)
  os.system("wigToBigWig -clip " + fTmp + " " + sizesFile + " " + bwFile)
  os.system("rm " + fTmp)

for f in os.listdir(wigGzDir):
#  if re.search(r'.combined.', f):
#    fNew = re.sub(r'.wig.gz$', '.bw', f)
#    print("converting " + f + " to " + fNew)
  m = re.match(r'^([^:]+):Developmental-Stage=([^#]+)#.*:ChIP-seq:.*:(ChIP|input):.*aligned_linear_10bp.*', f)
  if m:
    tf = m.group(1)
    stage = m.group(2)
    chipOrInput = m.group(3)
    print(tf + " " + stage + " " + chipOrInput)
    os.system("mkdir -p " + bwDir + "/" + stage)
    outputFile = bwDir + "/" + stage + "/" + tf + "." + chipOrInput + ".bw"
    wigGzConvert(wigGzDir + "/" + f, outputFile)

