# Fiji script to merge stacks.
#
# Loosely based on code at
# http://www.ini.uzh.ch/~acardona/fiji-tutorial/#RGB-stack-to-hyperstack

import sys
from ij import IJ, ImagePlus, ImageStack, CompositeImage
from ij.io import FileSaver
from ij.process import LUT

from jarray import array

from java.awt import Color

# from ij import ImageStack

if 0:
  if len(sys.argv) < 2:
    print "need more args\n";
    exit(1)

  print "first arg = " + sys.argv[1] + "\n"

imgPath = "/home/jburdick/gcb/data/image/FISH/20150208_Elicia_hlh-6_GFPprobe/"
outputPath = "mergedImages";

# read in images
img1 = IJ.openImage(imgPath + "/" + "alexa005.tif")
img2 = IJ.openImage(imgPath + "/" + "dapi005.tif")


ip1 = ImagePlus("ip1", img1.getImageStack())
ip2 = ImagePlus("ip2", img2.getImageStack())

# XXX temp hack
ip1.setCalibration(img1.getCalibration().copy())  

# create stack containing the images
stack1 = ip1.getImageStack()
stack2 = ip2.getImageStack()


stackNew = ImageStack(img1.width, img1.height)

for i in xrange(1, img1.getNSlices()+1):  

  # Extract the channels as FloatProcessor  
  i1 = stack1.getProcessor(i).toFloat(0, None)
  i2 = stack2.getProcessor(i).toFloat(0, None)

  # Add both to the new stack
  stackNew.addSlice("ch1 " + str(i), i1)
  stackNew.addSlice("ch2 " + str(i), i2)

imp2 = ImagePlus("merged", stackNew)
imp2.setCalibration(img1.getCalibration().copy())  

print(imp2.toString())
nChannels = 2              # two color channels  
nSlices = stackNew.getSize() # the number of slices of the original stack  
nFrames = 1                # only one time point   
imp2.setDimensions(nChannels, img1.getNSlices(), nFrames)
print(imp2.toString())

comp = CompositeImage(imp2, CompositeImage.COMPOSITE)  
# comp.setCalibration(img1.getCalibration().copy())

from jarray import array
lut1 = comp.createLutFromColor(Color.red)
lut2 = comp.createLutFromColor(Color.blue)


luts = array([lut1, lut2], LUT)

print(luts)

print(comp.getLuts())

# print comp.createLutFromColor(Color.blue)

comp.setLuts(luts)
comp.resetDisplayRanges()

# IJ.save(comp, "/home/jburdick/merged.tif")
fs = FileSaver(comp)
fs.saveAsTiffStack("/home/jburdick/merged.tif")

comp.show() 

