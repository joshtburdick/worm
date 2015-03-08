# Fiji script to merge stacks of FISH, one stack per embryo.
#
# Loosely based on code at
# http://www.ini.uzh.ch/~acardona/fiji-tutorial/#RGB-stack-to-hyperstack

import os
import sys
from ij import IJ, ImagePlus, ImageStack, CompositeImage
from ij.io import FileSaver
from ij.process import LUT

# directory containing images
imgPath = "/media/jburdick/disk2/data/FISH/2015_02_24 ElegansEmbryos/"

# directory in which to write output
outputDir = "/media/jburdick/disk2/data/FISH/merged/"

# names of the channels to add, in order of colors given in
# the CompositeImage source code (since I haven't been able
# to set LUTs "by hand")
# XXX omitting "trans", as it has a different number of planes
# channelNames = ["cy", "gfp", "dapi", "trans", "alexa", "tmr"]
channelNames = ["cy", "gfp", "dapi", "alexa", "tmr"]

# number of stacks to include
numStacks = 10   # XXX debugging

# IJ.debugMode = True
try:
  os.makedirs(outputDir)
except OSError:
  print "failure in making directory\n"

# loop through the stacks
for s in xrange(1, numStacks+1):

  # the IJ objects; hopefully we'll be able to close
  # these, and free up memory
  ijs = {}

  # the new stack, and the number of slices
  # (both set after we see an image)
  newStack = None
  numSlices = 0

  # pad the "stack" string with 0's
  if s < 10:
    s1 = "00" + str(s)
  else:
    s1 = "0" + str(s)

  # read in all the channels for this stack
  imgStack = {}
  for ch in channelNames:
    f = imgPath + "/" + ch + s1 + ".tif"
    print("reading " + f)

    ij = IJ.openImage(f)
    ijs[ch] = ij

    print(ij)
    imgStack[ch] = ij.getImageStack()
    numSlices = ij.getNSlices()

    # possibly create the stack to hold the new images
    # (now that we know its size)
    if (newStack is None):
      newStack = ImageStack(ij.width, ij.height)

  # add the slices
  for slice in xrange(1, numSlices+1):  
    for ch in channelNames:

      # Extract the channel as FloatProcessor
      ip = imgStack[ch].getProcessor(slice).toFloat(0, None)

      # add to the new stack
      newStack.addSlice(ch + " " + str(slice), ip)

  # create the merged stack
  imp = ImagePlus("merged " + s1, newStack)

  numChannels = len(channelNames)
  numSlicesTotal = newStack.getSize()
# imp.setDimensions(len(channelNames), numSlicesTotal, 1)
#    numSlicesTotal / len(channelNames))
  imp.setDimensions(len(channelNames), numSlices, 1)

  comp = CompositeImage(imp, CompositeImage.COLOR)
#  comp.setCalibration(img1.getCalibration().copy())
  IJ.save(comp, outputDir + s1 + ".tif")

  # free up images we saved
  for ch in channelNames:
    ijs[ch].close()
  imp.close()
  comp.close()

