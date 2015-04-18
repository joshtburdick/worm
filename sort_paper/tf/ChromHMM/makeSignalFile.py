#!/usr/bin/python3
# Converts bigWig files of histone marks to the 'signal file'
# format expected by ChromHMM.

import math
import subprocess

# Number of bp to use for the window.
windowSize = 1000000

# Gets sizes of chromosomes from a text file, so as to know
# how many windows there are in each chromosome.
chromSizesFile = '/home/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.fa.sizes'
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
  print(s)
  s = s.replace('n/a', '0')
  s = s.replace('\n', '')
  return(s.split('\t'))

# quick test
f = "/media/jburdick/disk2/histone_chip_seq_bw_old/Early-Embryos/H3_rep1_ChIP.bw"
a = getCounts(f, "I")
print(a)

