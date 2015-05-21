#!/usr/bin/python
# For each occurence of a given motif, finds the nearest
# gene it's upstream of (and the distance to that gene),
# and the conservation at that motif.

# Requires bedtools and Jim Kent's software to be in $PATH.

import os
import re
import subprocess
import sys

motifBamPath = sys.argv[1]
outputDir = sys.argv[2]

# bigWig file of conservation
# XXX probably should rename chromosomes in this
conservationBigWig = "/var/tmp/data/WS220.phastCons7way.bw"

# BED file containing the locations of the upstream regions
# (presumably one per gene, not filtered by conservation)
upstreamBed = "/home/jburdick/gcb/git/tf/motif/motifCount/regions/upstream_5kb_WS220.bed";

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
        "\t" + str(i) + "\t" + s[4] + "\t" + s[5] + "\n")
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

  # then, compute conservation at those locations
  subprocess.call(["bigWigAverageOverBed",
    conservationBigWig, "tmpMotif2.bed",
    "tmpMotifCons1.tsv",  # the TSV file isn't actually generated, AFAIK
    "-bedOut=tmpMotifCons1.bed"])

  # compute overlaps of those with upstream intergenic regions
  subprocess.call(["bedSort",
    "tmpMotifCons1.bed", "tmpMotifCons1.bed"])
  subprocess.call(["bedtools",
    "intersect",
    "-wa", "-wb",
    "-a", upstreamBed,
    "-b", "tmpMotifCons1.bed"],
    stdout=open(outputDir + "/" + name + "_upstreamMotifCons.tsv", "w"))

  # compress output
  subprocess.call(["gzip", "-f",
    outputDir + "/" + name + "_upstreamMotifCons.tsv"])

  # clean up
  subprocess.call(["rm", "tmpMotif1.bed", "tmpMotif2.bed",
    "tmpMotifCons1.tsv", "tmpMotifCons1.bed"])

subprocess.call(["mkdir", "-p", outputDir])

if True:
  for f in os.listdir(motifBamPath):
    if re.match(".*.bam$", f):
      m = f.replace(".bam", "")
      print(m)
      computeMotifDistAndConservation(m)

