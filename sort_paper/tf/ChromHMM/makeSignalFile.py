#!/usr/bin/python3
# Converts bigWig files of histone marks to the 'signal file'
# format expected by ChromHMM.

import math
import os
import re
import subprocess

import pdb

# Number of bp to use for the window.
windowSize = 50

# Size of each chromosome.
chromSizesFile = '/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.sizes'

# Directory containing the bigWig files to include.
# bwDir = "/media/jburdick/disk2/histone_chip_seq_bw_new/Early-Embryos/"
bwDir = '/murrlab/seq/igv/histone.ChIP.seq/'

# Directory in which to write output
outputDir = "signalFiles/"

os.system("mkdir -p " + outputDir)

# Gets sizes of chromosomes from a text file, so as to know
# how many windows there are in each chromosome.
chromSizes = {}
for line in open(chromSizesFile):
  a = line.replace('\n', '').split('\t')
  chromSizes[ a[0] ] = int(a[1])

# Gets data for one bigWig file.
# Args:
#   bigWigFile - the file from which to get counts
#   chrom - which chromosome to read
# Returns: an array of counts from that file
def getCounts(bigWigFile, chrom):
  numBins = math.floor( chromSizes[chrom] / windowSize )
  b = windowSize * numBins
  # XXX universal_newlines means this will return a string, not a buffer
  s = subprocess.check_output(["bigWigSummary", bigWigFile,
    chrom, "0", str(b), str(numBins)], universal_newlines=True)
  s = s.replace('n/a', '0')
  s = s.replace('\n', '')
  return([ math.floor(float(a)) for a in s.split('\t') ])

# Writes counts in a set of bigWig files.
# Args:
#   f - file object to write counts to
#   bwFiles - list of filenames
#   chrom - the chromosome to write
# Side effects: writes to that file.
def writeCounts(f, bwFiles, chrom):

  # read in all the counts
  x = [getCounts(b, chrom) for b in bwFiles] 
  # write them out  
  m = len(x[0])
  n = len(x)
  for i in range(0, m):
    for j in range(0, n-1):
      f.write( str(x[j][i]) + '\t' )
    f.write( str(x[n-1][i]) + '\n' )

# Writes out one 'signal file' of counts.
# Args:
#   experimentNames - (hopefully readable) names of the experiments
#   chipInputPairs - list of pairs of files
#   outBase - base name of files to write out
#   cell - name of the cell
# Side effects: writes the file.
def writeSignalFiles(experimentNames, chipInputPairs, outBase, cell):
  for chrom in ['I', 'II', 'III', 'IV', 'V', 'X']:

    # first, the signal
    f = open(outputDir + '/' + outBase + "_" + chrom + '_signal', 'w')
    f.write(cell + '\t' + chrom + '\n')
    f.write('\t'.join(experimentNames) + '\n')
    writeCounts(f, [bwDir + x["ChIP"] for x in chipInputPairs], chrom)
    f.close()

    # then, the control
    f = open(outputDir + '/' + outBase + '_' + chrom + '_controlsignal', 'w')
    f.write(cell + '\t' + chrom + '\n')
    f.write('\t'.join(experimentNames) + '\n')
    writeCounts(f, [bwDir + x["input"] for x in chipInputPairs], chrom)
    f.close()

# Gets paired ChIP and input experiment names from a
# list of filenames.
# Args:
#   files - a list of filenames
# Returns: a list of records with fields "ChIP" and "input".
def getChIPAndInputNames(files):

  # get the chip files, and corresponding input file
  chipFiles = [a for a in files if ' ChIP ' in a]
  filePairs = [{'ChIP': a,
    'input': a.replace(' ChIP ', ' input ')}
    for a in chipFiles]
  filePairs = [f for f in filePairs
    if f["ChIP"] in files
    and f["input"] in files]
  return(filePairs)

# Gets experiment names.
# Args:
#   a - a list of files
# Returns: 
def getExperimentNames(a):
  chipExperiments = [ f.replace('_rep1_ChIP.bw', '') for f in a if '_rep1_ChIP.bw' in f ]
  inputExperiments = [ f.replace('_rep1_input.bw', '') for f in a if '_rep1_input.bw' in f ]
  r = [ e for e in chipExperiments if e in inputExperiments ]
  r.sort()
  return(r)

# XXX hack to include two directories' worth of files
a = (["Early-Embryos/" + x for x in os.listdir(bwDir + "/Early-Embryos")] +
  ["Embryos-mixed-stages/" + x for x in os.listdir(bwDir + "/Embryos-mixed-stages")])
chipInputPairs = getChIPAndInputNames(a)
experimentNames = [a["ChIP"].replace(".bw", "").replace("Early-Embryos/", "EE ").replace("Embryos-mixed-stages/", "Emixed ") for a in chipInputPairs]

writeSignalFiles(experimentNames, chipInputPairs, "Embryo_50bp", "Embryo")


