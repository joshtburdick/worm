#!/usr/bin/python
# For each occurence of a given motif, finds the nearest
# gene it's upstream of (and the distance to that gene),
# and the conservation at that motif.

import subprocess
import sys

# where to write results
outputDir = "distAndConservation/"

# where bigWig files of motif locations are
motifBigWigPath = "/murrlab/seq/igv/motif/meme/"

# bigWig file of conservation
# XXX probably should rename chromosomes in this
conservationBigWig = "/murrlab/seq/igv/conservation/ce10.phastCons7way.bw"

# BED file containing the locations of the upstream regions
# (presumably one per gene, not filtered by conservation)
upstreamBed = "/home/jburdick/gcb/git/tf/motif/motifCount/upstreamRegionsWS220_5kb_nogenes.bed"

# paths to various tools (??? just put these in PATH?)
bedtoolsPath = "/home/jburdick/bin/x86_64/"
kentToolsPath = "/home/jburdick/bin/x86_64/"

# Gives each line in a BED file a unique name.
def nameMotifs(inputBedFile, outputBedFile):
  inFile = open(inputBedFile, "r")
  outFile = open(outputBedFile, "w")
  i = 0
  for line in inFile:
    s = line.split("\t")
    # widen intervals slightly
    a = int(s[1]) - 1
    b = int(s[2])
    # also rename chromosome slightly
    chr = "chr" + s[0]
    if s[0] == "MtDNA":
      chr = "chrM"
    outFile.write(chr + "\t" + str(a) + "\t" + str(b) +
      "\t" + str(i) + "\n")
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

  # first, unpack motif bigWig to a bigBed file
  subprocess.call([kentToolsPath + "/bigWigToBedGraph",
    motifBigWigPath + "/" + name + ".bw",
    "tmpMotif1.bed"])

  # give each line in this a unique name
  nameMotifs("tmpMotif1.bed", "tmpMotif2.bed")

  # then, compute conservation at those locations
  subprocess.call([kentToolsPath + "/bigWigAverageOverBed",
    conservationBigWig, "tmpMotif2.bed",
    "tmpMotifCons1.tsv",
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

# computeMotifDistAndConservation("FOXB1_full")

# get list of motifs to include
motifs = subprocess.check_output("cut -f1 /home/jburdick/gcb/git/cluster/hierarchical/hier.200.clusters/uniqueKnownMotifEnrichment_5kb_0.5cons.tsv | tail --lines=+2 | sort | uniq",
  shell=True).split("\n")

for m in motifs:
  if (m != ""):
    print("motif: " + m)
    computeMotifDistAndConservation(m)

