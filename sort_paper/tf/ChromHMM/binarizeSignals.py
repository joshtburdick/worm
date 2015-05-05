#!/usr/bin/python3
# Binarizes the signals from ChromHMM.

import os
import subprocess

jar = '/home/jburdick/lib/java/ChromHMM/ChromHMM.jar'

# number of bp in a window
windowSize = 50

# sizes of chromosomes, for wigToBigWig
sizesFile = "/murrlab/jburdick/data/seq/Caenorhabditis_elegans.WS220.64.dna.toplevel.sizes"

# Writes out one bigWig file.
# Args:
#   signalFileBase - base name (including directory) of signal file
#   signalColumn - index of the column to write
#   outputDir - where to write output
# Side effects: writes a bigWig file in outputDir
def writeWiggle1(signalFileBase, signalColumn, outputDir):
  fTmp = "tmp.wig"

  # get name of signal column
  with open(signalFileBase + "_I_binary.txt") as signalFile:
    s = signalFile.readline()
    s = signalFile.readline().strip("\n")
    signal = s.split("\t")[signalColumn]

  # write out the wiggle file
  with open("tmp.wig", "w") as wigFile:
    wigFile.write('track type=wiggle_0 name="' + signal + '" \n')

    # loop through the chromosomes
    for chr in ["I", "II", "III", "IV", "V", "X"]:

      # get signal data. XXX we read in the whole file, which is inefficient,
      # but easier to think about
      with open(signalFileBase + "_" + chr + "_binary.txt") as signalFile:
        s = [x.strip('\n') for x in signalFile.readlines()]
      s = s[2:]

      wigFile.write("fixedStep chrom=" + chr + " start=1" +
        " step=" + str(windowSize) + " span=" + str(windowSize) + "\n")

      for a in s:
        wigFile.write(a.split("\t")[signalColumn] + "\n")
  
  # convert .wig to .bw
  bwFile = outputDir + "/" + signal + ".bw"
  r = subprocess.call(["wigToBigWig", "-clip", fTmp, sizesFile, bwFile])
  if (r==0):
    os.system("rm " + fTmp)
  else:
    print("wigToBigWig failed for " + bwFile)

# binarize the data
# subprocess.call(['java', '-jar', jar, 'BinarizeSignal', '-c', 'signalFiles',
#  'signalFiles', 'binarized']) 

os.makedirs("binarizedBigWig", exist_ok=True)

# as a check, write out wiggle files of the binarized signals
writeWiggle1("binarized/EarlyEmbryo", 0, "binarizedBigWig")


