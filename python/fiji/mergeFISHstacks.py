# Fiji script to merge stacks.
#
# Loosely based on code at
# http://www.ini.uzh.ch/~acardona/fiji-tutorial/#RGB-stack-to-hyperstack

import sys
from ij import IJ, ImagePlus, ImageStack, CompositeImage
from ij.io import FileSaver
from ij.process import LUT

imgPath = "/home/jburdick/gcb/data/image/FISH/20150208_Elicia_hlh-6_GFPprobe/"

# IJ.debugMode = True

# names of the channels to add, in order of colors given in
# the CompositeImage source code (since I haven't been able
# to set LUTs "by hand")
channelNames = ["cy", "gfp", "dapi", "trans", "alexa", "tmr"]

# read in images
img1 = IJ.openImage(imgPath + "/" + "alexa005.tif")
img2 = IJ.openImage(imgPath + "/" + "dapi005.tif")


ip1 = ImagePlus("ip1", img1.getImageStack())
ip2 = ImagePlus("ip2", img2.getImageStack())

# XXX temp hack
# ip1.setCalibration(img1.getCalibration().copy())  

# create stack containing the images
stack1 = ip1.getImageStack()
stack2 = ip2.getImageStack()

stackNew = ImageStack(img1.width, img1.height)

for i in xrange(1, 4):  # img1.getNSlices()+1):  
  for(ch in
  # Extract the channels as FloatProcessor  
  i1 = stack1.getProcessor(i).toFloat(0, None)
  i2 = stack2.getProcessor(i).toFloat(0, None)

  # Add both to the new stack
  stackNew.addSlice("ch1 " + str(i), i1)
  stackNew.addSlice("ch2 " + str(i), i2)

imp2 = ImagePlus("merged", stackNew)
# imp2.setCalibration(img1.getCalibration().copy())  

print(imp2.toString())
nChannels = 6
nSlices = stackNew.getSize() # the number of slices of the original stack  
nFrames = 10             # number of time points
imp2.setDimensions(nChannels, stackNew.getSize(), nFrames)
print(imp2.toString())
print(imp2.getStatistics())

comp = CompositeImage(imp2, CompositeImage.COLOR)  
comp.setCalibration(img1.getCalibration().copy())

# comp.show()

IJ.save(comp, "/home/jburdick/merged.tif")

# XXX various attempts at saving the result

# fs = FileSaver(comp)
# fs.saveAsTiffStack("/home/jburdick/merged.tif")

# IJ.runPlugIn("loci.plugins.LociExporter",
#   "imageid=" + str(img1.getID()) + " " +
#   "outfile=/home/jburdick/merged.tif " +
#   "splitz=false splitc=false splitt=false saveroi=false ")  # +
#  "compression-type=Uncompressed")



