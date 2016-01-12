package edu.murrlab.image.translate;

import java.util.*;

import java.awt.image.*;
import java.io.*;
import java.util.zip.*;
import java.util.regex.*;
// import javax.media.jai.*;
// import com.sun.media.jai.codec.*;

/* javadoc for these is at 
http://hudson.openmicroscopy.org.uk/job/BIOFORMATS-trunk/javadoc/? */
import loci.common.*;
import loci.common.services.*;
import loci.formats.*;
import loci.formats.in.*;
import loci.formats.meta.IMetadata;
import loci.formats.ome.*;
import loci.formats.out.*;
import loci.formats.services.OMEXMLService;

import ome.xml.model.primitives.*;
// FIXME if I rewrite this to use a more recent version of
// loci_tools.jar, then I'll need to use this
// import ome.units.quantity.Length;

/** Imports a .lif file (hopefully) into a multiplane TIFF.
 * Hopefully this will be easier to deal with than a large number
 * of small files (and AceTree can deal with this.)
 * 
 * Heavily based on the GetPhysicalMetadata example.
 * @author jburdick
 *
 */
public class LeicaImporterMultiplane {

	/** Directory to store output directories (by default, this is
	 * the same directory that contains the .lif file.) */
	String outputPath = null;
	
	/** The name of the series. */
	String seriesName = null;
	
	/** Stores metadata. */
	OMEXMLMetadata meta;
	
	/** The raw XML data. (Probably it's better to use the OMXML file
	 * when possible, but this should at least theoretically contain
	 * all the annotation in the .lif file.) */
	private String rawAnnotation;
	
	/** Reads image info. */
	ImageReader reader;
	
	/** Number of series. */
	int seriesCount;

	/** Directories in which to put the two channels. */
	String[] channelName = {"green", "red"};
	
	/** The "job ID", set by getImageNames(). */
	// String jobID = null;   probably not used
	
	/** The XY position strings, set by getImageNames(). */
	Vector<String> xyPosition = new Vector<String>();
	
    private ServiceFactory factory;

    private OMEXMLService service;

    /** File of images we're currently writing out. */
    private ImageWriter writer;
    
    /** Buffer for images (which we avoid allocating more than once.) */
    private byte[] img;
    
    /** Constructor. */
	public LeicaImporterMultiplane(String lifFile)
			throws DependencyException, FormatException, IOException, ServiceException {
		
		// set some path names
		File f = new File(lifFile);
		this.seriesName = f.getName().replaceAll(".lif$", "");
		this.outputPath = f.getAbsoluteFile().getParent();
				
	    // read in all the XML metadata
	    LeicaLIFXMLReader lr = new LeicaLIFXMLReader(new File(lifFile));
	    this.rawAnnotation = lr.getXML();
	    
	    factory = new ServiceFactory();
	    service = factory.getInstance(OMEXMLService.class);
	    
	    // create OME-XML metadata store	    
	    meta = service.createOMEXMLMetadata();

	    // create format reader
	    reader = new ImageReader();
		reader.setMetadataStore(meta);	    

	    // initialize file
	    System.out.println("Initializing " + lifFile);
	    reader.setId(lifFile);

	    seriesCount = reader.getSeriesCount();
	}

	/** Reads through the images, and sets jobID and xyPosition. */
	public void getImageNames() {
		
		// count of each sequence name
		HashMap<String, Integer> count = new HashMap<String, Integer>();
		
		// loop through the series
		for(int s = 0; s < seriesCount; s++) {
			
			// get image name
			reader.setSeries(s);
			Hashtable<String, Object> h = reader.getSeriesMetadata();
			String imageName = (String) h.get("Image name");
			
			// record the XY position (used for separating out
			// images from different embryos)
			String xy = getXY(imageName);
			if (xy != null && !xyPosition.contains(xy))
				xyPosition.add( xy );
		}
	}
	
	/** Writes out all of the series. */
	public void writeAllSeries() throws IOException, FormatException, ServiceException {
	
		int embryoNumber = 1;
		
		// loop through the XY position strings
		for (String xy : xyPosition) {
			writeEmbryo(xy, embryoNumber);
			embryoNumber++;
		}
	}
	
	/** Test of whether a given image name is one that we should
	 * write out.
	 * @param imageName  name of the image
	 * @param xy the "coordinates" of the image */
	private boolean isSeriesToWrite(String imageName, String xy) {
		return !imageName.contains("Sequence/Image") &&
						!imageName.contains("Sequence/AF") &&
						!imageName.contains("AF Job") &&
						!imageName.contains("DriftAF") &&
						xy.equals( getXY(imageName) );
	}
	
	/** Gets the size of one embryo's images (number of slices, and
	 * timepoints.) This is needed because we need to specify this
	 * when we open the output image stack.
	 * XXX this is sort of annoyingly slow.
	 * 
	 * @param xy  coordinates of the stack to get
	 * 
	 * @return four-element array, containing z and t dimensions
	 */
	private int[] getEmbryoStackDimensions(String xy) {
		// the X, Y, Z, and T dimensions
		int[] r = {1, 1, 1, 0};
			
		// loop through the series
		System.out.println("seriesCount = " + seriesCount);
		for(int s = 0; s < seriesCount; s++) {
			
			// get name of the image
			reader.setSeries(s);
			Hashtable<String, Object> h = reader.getSeriesMetadata();
			String imageName = (String) h.get("Image name");
			
			// if this is part of the appropriate sequence, include it
			// in the dimensions
			if (isSeriesToWrite(imageName, xy)) {
				System.out.println("including " + imageName);
				r[0] = reader.getSizeX();
				r[1] = reader.getSizeY();
				if (reader.getSizeZ() > r[2])
					r[2] = reader.getSizeZ();
				r[3]++;		
			}	
		}		
		
		return r;
	}
	
	/** Writes one embryo's worth of movies (along with the metadata.)
	 * @param xy  the "XY position string" to include
	 * @param embryoNumber  which embryo this is
	 * @throws FormatException
	 * @throws IOException
	 */
	public void writeEmbryo(String xy, int embryoNumber)
			throws IOException, FormatException, ServiceException {
		String outputPath1 = outputPath + "/" + seriesName + "_L" + embryoNumber;
		System.out.println("writing " + outputPath1 + "   XY = " + xy);
		
		// get dimensions
		int[] sizeXYZT = getEmbryoStackDimensions(xy);
		System.out.println("size: X=" + sizeXYZT[0] + " Y=" + sizeXYZT[1] +
				" Z=" + sizeXYZT[2] + " T=" + sizeXYZT[3]);
		
		// allocate the image buffer (this assumes one byte per pixel)
		img = new byte[ sizeXYZT[0] * sizeXYZT[1] ];
		
		// create directory, and open merged TIF file
		new File(outputPath1 + "/dats").mkdirs();
		
		String multiplaneTiffName = outputPath1 + "merged.tif";
		openImageWriterFile(multiplaneTiffName, sizeXYZT[2], sizeXYZT[3]);
	
		// we copy the metadata into each embryo directory, even though it
		// seems to be the same for each series
		writeMetadata(outputPath1 + "/dats/");
		
		// file to write a list of images in this series
		PrintWriter imageList = new PrintWriter(outputPath1 + "/dats/imageList.tsv");
		imageList.println("time\timage");
		
		// counter for current time point (now 0-based)
		int t = 0;
		
		// loop through the series
		for(int s = 0; s < seriesCount; s++) {
			
			// get name of the image
			reader.setSeries(s);
			Hashtable<String, Object> h = reader.getSeriesMetadata();
			String imageName = (String) h.get("Image name");
			
			// if this is part of the appropriate sequence, write it
			// ("names of things to skip" is similar to what's in
			// LeicaTifRename.pl)
			if (isSeriesToWrite(imageName, xy)) {
				System.out.println("writing image " + imageName + "   XY = " + getXY(imageName));
				imageList.println(t + "\t" + imageName);
				writeSeries(outputPath1, embryoNumber, s, t);
				t++;
			}	
		}
		
		imageList.close();
		writer.close();
	}

	/** Writes one series (which is actually only one timepoint.)
	 * @param outputPath  directory in which to write the file
	 * @param embryoNumber  which embryo this is (so as to add
	 *   "_L1", "_L2", "_L3", etc. to the name)
	 * @param series  the index of the series in the file
	 * @param t  the (0-based) timepoint number with which to label this. */
	private void writeSeries(String outputPath, int embryoNumber,
			int series, int t) throws IOException, FormatException, ServiceException {
		reader.setSeries(series);

		int sizeZ = reader.getSizeZ();
		int maxZ = sizeZ - 1;    // this is used for reversing the z-stacks
		
		// total number of planes
		int planesPerTimeSlice = 2 * sizeZ;
		
		System.out.println("planesPerTimeSlice = " + planesPerTimeSlice);
		
		// because the order in which to write these may be different,
		// we store the images temporarily
		// probably not used: byte[][] imgs = new byte[ 2 * numPlanes ][];
		// byte[] img = new byte[ reader.getSizeX() * reader.getSizeY() ];
		
		// loop through the images in the series
		for(int i = 0; i < reader.getImageCount(); i++) {
			
			// get channel and plane info
			int channel = meta.getPlaneTheC(series, i).getValue();
			
			// note that this reverses the z-stacks
			int plane = maxZ - meta.getPlaneTheZ(series, i).getValue();
			
			// this may not be used
			String imageName = channelName[ channel ] + "_" +
					seriesName + "_L" + embryoNumber +
					"-t" + String.format("%03d", (t+1)) +
					"-p" + String.format("%02d", plane + 1);
			
			// 0-based index of the image to write out.
			int imgIndex = planesPerTimeSlice * t + 2 * plane + channel;
			
			// store this image
//			imgs[imgIndex] = new byte[ numBytesPerImage ];
//			reader.openBytes(i, imgs[imgIndex]);
			reader.openBytes(i, img);
			writer.saveBytes(imgIndex % 100, img);   // XXX testing again
		}
		
		// actually write out the planes
		/*
		for(int p=0; p<imgs.length; p++) {
			int imageIndexInStack = (2*numPlanes) * t + p;
			writer.saveBytes(imageIndexInStack, imgs[p]);
			System.out.println("saved plane " + imageIndexInStack);
		}
		*/
	}
	
	/** Opens a file using the ImageWriter; this should be called before
	 * starting to write out images.
	 * @param filename the name of the output file
	 * @param numZ number of z stacks
	 * @param numT number of timepoints */
	private void openImageWriterFile(String filename, int numZ, int numT)
			throws IOException, FormatException, ServiceException {
	
		// this assumes the reader is already open
		int sizeX = reader.getSizeX();
		int sizeY = reader.getSizeY();
		
		// XXX this is somewhat complicated; it's based on 
		// loci-tools' components/bio-formats/utils/MinimumWriter.java
		IMetadata meta = service.createOMEXMLMetadata();
	    MetadataTools.populateMetadata(meta, 0, (String) null, false, "XYZCT",
	    	      FormatTools.getPixelTypeString(FormatTools.UINT8),
	    	      sizeX, sizeY, numZ, 2, numT, 1);

	    // FIXME loci_tools v.5.1.x uses a new API here, which I don't
	    // yet understand (workaround for now is to use v.5.0.x)
	  	meta.setPixelsPhysicalSizeX(new PositiveFloat(1.0), 0);
	    meta.setPixelsPhysicalSizeY(new PositiveFloat(1.0), 0);
	    
	    // storing in little-endian order (as that's what the Perl script
	    // does; although it may not matter)
	    meta.setPixelsBinDataBigEndian(Boolean.FALSE, 0, 0);
	    
   		writer = new ImageWriter();
   		writer.setMetadataRetrieve(meta);
   		writer.setCompression(TiffWriter.COMPRESSION_LZW);
		writer.setId(filename);
		
		System.out.println("can do stacks = " + writer.canDoStacks());
		writer.setWriteSequentially(false);
	}
	
	/** Writes out the metadata for the current series.
	 * ??? gzip this? */
	private void writeMetadata(String outputPath) throws IOException {
		
		// print the OMXML annotation (this is an open standard, but
		// leaves out some of the annotation)
		PrintStream ps = new PrintStream(new GZIPOutputStream(
				new FileOutputStream(outputPath + "/omxml.xml.gz")));
		ps.println(meta.dumpXML());
		ps.close();
		
		// print the raw annotation
		ps = new PrintStream(new GZIPOutputStream(
				new FileOutputStream(outputPath + "/rawAnnotation.xml.gz")));
		ps.println(rawAnnotation);
		ps.close();
	}
	
	/** Utility to convert from an image name to "JobID"
	 * (name of a subgroup of images.) */
	public String getJobID(String s) {
		// strip off everything after the last '/'
		return Pattern.compile("/.*$").matcher(s).replaceAll("/");
	}

	/** Parses out the "xy location" from an image name.
	 * @param imageName  e.g., "Sequence/A 0_S 0_U 0_V 0_X 0_Y 0_Job 1_002"
	 * @return in this case, "X 0_Y 0" */
	private String getXY(String imageName) {
		// ??? arguably this pattern could be written a bit differently, to specifically
		// match either, e.g., "U 1" or "U11"
		Pattern p = Pattern.compile(".*/A ?\\d+_S ?\\d+_U ?\\d+_V ?\\d+_(X ?\\d_Y ?\\d)_.*");
		Matcher m = p.matcher(imageName);
		if (!m.matches())
			System.err.println("failed to parse " + imageName);
		return (m.matches() ? m.group(1) : null);
	}
		
	public static void main(String[] args) throws Exception {
		System.out.println("started at " + new Date().toString());
		// XXX print errors to console, or a log file?
		// org.apache.log4j.BasicConfigurator.configure();
		// XXX testing
		if (Math.cos(0) < 0) {
			if (args.length != 1) {
				System.err.println("imports Leica images");
				System.err.println("Usage: leica_import.pl file.lif");
				System.err.println("writes series in the same directory");
				System.exit(1);
			}
			String lifFile = args[0];
		}
		String lifFile = "/var/tmp/20151027_tlp-1_enh4K.lif";

		System.out.println("converting " + lifFile);
    	LeicaImporterMultiplane li = new LeicaImporterMultiplane(lifFile);
    	li.getImageNames();
      	li.writeAllSeries();
      	System.out.println("finished at " + new Date().toString());
	}
}
