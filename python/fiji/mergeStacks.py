# Fiji script to merge stacks.
#
# Loosely based on code at
# http://www.ini.uzh.ch/~acardona/fiji-tutorial/#RGB-stack-to-hyperstack

import sys
from ij import IJ, ImagePlus, ImageStack, CompositeImage
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
print "img1 = " + img1.toString() + "\n"



# create stack containing the images
stack1 = img1.getImageStack()
stack2 = img2.getImageStack()

stackNew = ImageStack(img1.width, img1.height)

for i in xrange(1, 4):     # img1.getNSlices()+1):  
  # Extract the channels as FloatProcessor  
  i1 = stack1.getProcessor(i).toFloat(0, None)
  i2 = stack2.getProcessor(i).toFloat(1, None)

  # Add both to the new stack
  stackNew.addSlice("ch1", i1)
  stackNew.addSlice("ch2", i2)

imp2 = ImagePlus("merged", stackNew)
imp2.setCalibration(img1.getCalibration().copy())  
# IJ.run(imp2, "Green", "")
print(imp2.toString())

nChannels = 2              # two color channels  
nSlices = stackNew.getSize() # the number of slices of the original stack  
nFrames = 2                # only one time point   
imp2.setDimensions(nChannels, nSlices, nFrames)  
comp = CompositeImage(imp2, CompositeImage.COMPOSITE)  
comp.show()  


