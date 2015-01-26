#!/usr/bin/python
# Renames ChIP peak calls from modENCODE, and translates them
# to BED format.

import gzip
import math
import os
import re
import subprocess

# ChIP-TF peaks
# inputDir = "/home/jburdick/tmp/chip_TF_peak/"
inputDir = "/home/jburdick/data/ftp/data.modencode.org/C.elegans/Transcriptional-Factor/ChIP-seq/computed-peaks_gff3/"

# where to write GFF3 files
# outputDir = "/home/jburdick/tmp/chip_TF_peak_bed/"
outputDir = "/murrlab/seq/igv/chip.TF/peak/"

try:
  os.makedirs(outputDir)
except OSError as exc:
  pass

# Parses a ChIP experiment's name.
def parseExperimentName(name):
  a = name.split(":")
  gene = a[0]
  replicate = a[3].replace("Rep-", "")
  f = a[1].split("#")
  stage = f[0].split("=")[1]
  return {'stage':stage, 'gene':gene, 'replicate':replicate}


# Converts a peak file (in .gff3.gz format) to a BED file.
# This includes converting the score to an integer.
def convertPeakGff3ToBed(inputGff3GzFile, outputBedFile):
  inFile = gzip.GzipFile(inputGff3GzFile, "r")
  outFile = open("tmp1.bed", "w")

  i = 0
  for line in inFile:
    # skip comments
    if (line.startswith("#")):
      continue
    s = line.split("\t")
    a = int(s[3])
    b = int(s[4])

    # rename chromosome slightly
    chr = "chr" + s[0]
    if s[0] == "MtDNA":
      chr = "chrM"

    # convert q-value to an integer
    s1 = float(s[5])
    if s1 > 1e-100:
      # note multiplication by 10, since these are often small
      score = - math.trunc( 10 * math.log10(float(s[5])))
    else:
      score = 1000      # max. allowed by BED standard

    # we arbitrarily put these on the "+" strand
    # for now, numbering these (although it arguably doesn't
    # look as nice in IGV)
    outFile.write(chr + "\t" + str(a) + "\t" + str(b) +
      "\t" + str(i) + "\t" + str(score) + "\t" + "+" + "\n")
    i = i + 1
  inFile.close()
  outFile.close()
  
  # and sort this
  # XXX shouldn't use shell=True
  subprocess.call("bedtools sort -i tmp1.bed > " + outputBedFile, shell=True)


for f in os.listdir(inputDir):
  if re.match(".*combined.*", f) and re.match(".*\.gff3\.gz", f):
#    print f
    n = parseExperimentName(f)
    # output filename is like the original filename,
    # but slightly simplified
    outputFile = (n['gene'] + "_" + n['stage']   # + "_rep" + n['replicate'] 
      + ".bed")
    print outputFile
    convertPeakGff3ToBed(inputDir + "/" + f,
      outputDir + "/" + outputFile)




