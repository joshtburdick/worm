#!/usr/bin/python3
# For each occurence of a given motif, finds the nearest
# gene it's upstream of (and the distance to that gene),
# but doesn't compute motif conservation. Designed for counting
# motifs in, e.g., C.briggsae.

# Requires bedtools and Jim Kent's software to be in $PATH.

import os
import re
import subprocess
import sys

motifBamPath = sys.argv[1]
outputDir = sys.argv[2]

# BED file containing the locations of the upstream regions
# (presumably one per gene, not filtered by conservation)
upstreamBed = sys.argv[3]

# Gives each line in a BED file a unique name.
# Also:
# - XXX omits any motifs with location <= 0
def nameMotifs(inputBedFile, outputBedFile):
  inFile = open(inputBedFile, "r")
  outFile = open(outputBedFile, "w")
  i = 0
  for line in inFile:
    s = line.split("\t")
    a = int(s[1])
    b = int(s[2])
    if int(s[1]) > 0:
      outFile.write(s[0] + "\t" + str(a) + "\t" + str(b) +
        "\t" + str(i) + "\t" + s[4] + "\t" + s[5])
  #"\t0\t+\t" +
  #      str(a) + "\t" + str(b) + "\t0\t1\t1\t0\n")
      i = i + 1
  inFile.close()
  outFile.close()

# Computes motif distances and conservation for one motif.
# Args:
#   motif - name of motif to compute results for
# Side effects: writes a file in outputDir, with columns FIXME
def computeMotifDistAndConservation(name):

  # first, unpack motif .bam to a bigBed file
  subprocess.call(["bedtools", "bamtobed",
    "-i", motifBamPath + "/" + name + ".bam"],
    stdout=open("tmpMotif1.bed", "w"))

  # give each line in this a unique name (and change chromosome names)
  nameMotifs("tmpMotif1.bed", "tmpMotif2.bed")

  # compute overlaps of those with upstream intergenic regions
  subprocess.call(["bedSort",
    "tmpMotif2.bed", "tmpMotif2.bed"])
  subprocess.call(["bedtools",
    "intersect",
    "-wa", "-wb",
    "-a", upstreamBed,
    "-b", "tmpMotif2.bed"],
    stdout=open(outputDir + "/" + name + "_upstreamMotifCons.tsv", "w"))

  # compress output
  subprocess.call(["gzip", "-f",
    outputDir + "/" + name + "_upstreamMotifCons.tsv"])

  # clean up
#  subprocess.call(["rm", "tmpMotif1.bed", "tmpMotif2.bed"])

subprocess.call(["mkdir", "-p", outputDir])

if True:
  for f in os.listdir(motifBamPath):
    if re.match(".*.bam$", f):
      m = f.replace(".bam", "")
      print(m)
      computeMotifDistAndConservation(m)

