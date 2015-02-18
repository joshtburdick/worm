#!/usr/bin/python
# For each occurence of a given motif, finds the nearest
# gene it's upstream of (and the distance to that gene),
# and the conservation at that motif.

import os
import re
import subprocess
import sys

# where to write results
# outputDir = "distAndConservation/5kb/"
# outputDir = "distAndConservation/5kb_de_novo/"
# outputDir = "distAndConservation/5kb_hughes/"
# outputDir = "distAndConservation/5kb_TCF_LEF/"
# outputDir = "/home/jburdick/tmp/distAndConservation/5kb_hughes_20141202/"
outputDir = "/home/jburdick/tmp/distAndConservation/jolma2013_shuffled/"

# where .bam files of motif locations are
# motifBamPath = "/murrlab/seq/igv/motif/known/"
# motifBamPath = "/home/jburdick/tmp/fimo_denovo/"
# motifBamPath = "/home/jburdick/tmp/fimo_hughes/"
# motifBamPath = "/murrlab/seq/igv/motif/TCF_LEF/"
# motifBamPath = "/home/jburdick/tmp/fimo_hughes_20141202/"
motifBamPath = "/home/jburdick/tmp/fimo/jolma2013_shuffled/"

# bigWig file of conservation
# XXX probably should rename chromosomes in this
conservationBigWig = "/murrlab/seq/igv/conservation/ce10.phastCons7way.bw"

# BED file containing the locations of the upstream regions
# (presumably one per gene, not filtered by conservation)
upstreamBed = "/home/jburdick/gcb/git/tf/motif/motifCount/regions/upstream_5kb_WS220.bed";

# paths to various tools (??? just put these in PATH?)
bedtoolsPath = "/home/jburdick/bin/x86_64/"
kentToolsPath = "/home/jburdick/bin/x86_64/"

# Gives each line in a BED file a unique name.
# Also:
# - renames chromosomes
# - XXX omits any motifs with location <= 0
def nameMotifs(inputBedFile, outputBedFile):
  inFile = open(inputBedFile, "r")
  outFile = open(outputBedFile, "w")
  i = 0
  for line in inFile:
    s = line.split("\t")
    a = int(s[1])
    b = int(s[2])
    # also rename chromosome slightly
    chr = "chr" + s[0]
    if s[0] == "MtDNA":
      chr = "chrM"
    if int(s[1]) > 0:
      outFile.write(chr + "\t" + str(a) + "\t" + str(b) +
        "\t" + str(i) + "\t" + s[4] + "\t" + s[5] + "\n")
  #"\t0\t+\t" +
  #      str(a) + "\t" + str(b) + "\t0\t1\t1\t0\n")
      i = i + 1
  inFile.close()
  outFile.close()

# Undoes the previous chromosome renaming.
def renameChromosomes(inputBedFile, outputBedFile):
  inFile = open(inputBedFile, "r")
  outFile = open(outputBedFile, "w")
  for line in inFile:
    s = line.split("\t")
    s[0] = s[0][3:100]
    if (s[0] == "M"):
      s[0] = "MtDNA"
    outFile.write("\t".join(s))
  inFile.close()
  outFile.close()

# Computes motif distances and conservation for one motif.
# Args:
#   motif - name of motif to compute results for
# Side effects: writes a file in outputDir, with columns FIXME
def computeMotifDistAndConservation(name):

  # first, unpack motif .bam to a bigBed file
  subprocess.call([bedtoolsPath + "/bedtools", "bamtobed",
    "-i", motifBamPath + "/" + name + ".bam"],
    stdout=open("tmpMotif1.bed", "w"))

  # give each line in this a unique name (and change chromosome names)
  nameMotifs("tmpMotif1.bed", "tmpMotif2.bed")

  # then, compute conservation at those locations
  subprocess.call([kentToolsPath + "/bigWigAverageOverBed",
    conservationBigWig, "tmpMotif2.bed",
    "tmpMotifCons1.tsv",  # the TSV file isn't actually generated, AFAIK
    "-bedOut=tmpMotifCons1.bed"])

  # XXX change chromosome names back
  renameChromosomes("tmpMotifCons1.bed", "tmpMotifCons2.bed")

  # compute overlaps of those with upstream intergenic regions
  subprocess.call([bedtoolsPath + "/bedSort",
    "tmpMotifCons2.bed", "tmpMotifCons2.bed"])
  subprocess.call([bedtoolsPath + "/bedtools",
    "intersect",
    "-wa", "-wb",
    "-a", upstreamBed,
    "-b", "tmpMotifCons2.bed"],
    stdout=open(outputDir + "/" + name + "_upstreamMotifCons.tsv", "w"))

  # compress output
  subprocess.call(["gzip", "-f",
    outputDir + "/" + name + "_upstreamMotifCons.tsv"])

  # clean up
  subprocess.call(["rm", "tmpMotif1.bed", "tmpMotif2.bed",
    "tmpMotifCons1.tsv", "tmpMotifCons1.bed", "tmpMotifCons2.bed"])

subprocess.call(["mkdir", "-p", outputDir])

# computeMotifDistAndConservation("FOXB1_DBD_3")

# get list of motifs to include
# motifs = subprocess.check_output("cut -f1 /home/jburdick/gcb/git/cluster/hierarchical/hier.200.clusters/uniqueKnownMotifEnrichment_5kb_0.5cons.tsv | tail --lines=+2 | sort | uniq",
#   shell=True).split("\n")

# nope, running it on everything
#motifs = [ f.replace(".bam", "")

# computeMotifDistAndConservation("bp.5kb.cons.hier.200clusters_186_5")

if True:
  for f in os.listdir(motifBamPath):
    if re.match(".*.bam$", f):
      m = f.replace(".bam", "")
      print(m)
      computeMotifDistAndConservation(m)

