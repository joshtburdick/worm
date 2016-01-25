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

/** Imports a .lif file (hopefully) into multiplane TIFFs, with
 * one stack per timepoint.
 * Hopefully this will be easier to deal with than a large number
 * of small files (and AceTree can deal with this.)
 *
 * @author jburdick
 */
public class LeicaImporterLSM {

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

	/** Directories in which to put the two channels (deprecated). */
	// String[] channelName = {"green", "red"};
	
	/** The "job ID", set by getImageNames(). */
	// String jobID = null;   probably not used
	
	/** The XY position strings, set by getImageNames(). */
	Vector<String> xyPosition = new Vector<String>();
	
	/** The dimensions of the embryo (X, Y, Z, and T, respectively.) */
	private int[] sizeXYZT = new int[4];
	
    private ServiceFactory factory;

    private OMEXMLService service;

    /** File of images we're currently writing out. */
    private ImageWriter writer;

    /** Buffer for images (which we avoid allocating more than once.) */
    private byte[][] imgs;
    
    /** Constructor. */
	public LeicaImporterLSM(String lifFile)
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
			/*
			Hashtable<String, Object> h = reader.getSeriesMetadata();
			String imageName = (String) h.get("Image name");
			*/
			// possibly faster way of getting the image name
			String imageName = (String) reader.getSeriesMetadataValue("Image name");
			
			// record the XY position (used for separating out
			// images from different embryos)
			String xy = getXY(imageName);
			if (xy != null && !xyPosition.contains(xy))
				xyPosition.add( xy );
		}
	}
	
	/** Writes out all of the series. */
	public void writeAllSeries() throws IOException, FormatException, ServiceException {
	
		// since numbering of these starts at 1
		// (e.g. names end with "_L1", "_L2", "_L3", etc.)
		int embryoNumber = 1;
		
		// loop through the XY position strings
		for (String xy : xyPosition) {
			writeEmbryo(xy, embryoNumber);
			embryoNumber++;
		}
	}
	
	/** Test of whether the current image (as set by
	 * Reader.setSeries()) is one that we should write out,
	 * based on the image name and dimensions (note that these
	 * should be set before calling this.)
	 * @param xy the "coordinates" of the image 
	 * @return */
	private boolean isSeriesToWrite(String xy) {
		String imageName = (String) reader.getSeriesMetadataValue("Image name");
		
		// first, check on the image name
		boolean r = (!imageName.contains("Sequence/Image") &&
						!imageName.contains("Sequence/AF") &&
						!imageName.contains("AF Job") &&
						!imageName.contains("DriftAF") &&
						xy.equals( getXY(imageName)));
		
		if (!r)
			return false;

		// if we're at the first timepoint, then don't check dimensions
		if (sizeXYZT[3] == 0)
			return true;
		
		// if we're past the first timepoint, check if the dimensions
		// match those of the first timepoint (and both channels are present)
		else return (reader.getSizeX() == sizeXYZT[0] &&
				reader.getSizeY() == sizeXYZT[1] &&
				reader.getSizeZ() == sizeXYZT[2] &&
				reader.getSizeC() == 2);
	}
	
	/** Gets the size of one embryo's images (number of slices, and
	 * timepoints.) This is needed because we need to specify this
	 * when we open the output image stack.
	 * XXX this is sort of annoyingly slow. In theory, we could read
	 * through the entire .lif file once, and get the dimensions for
	 * each series at the same time.
	 * 
	 * @param xy  coordinates of the stack
	 * 
	 * @return four-element array, containing z and t dimensions
	 */
	private void setEmbryoStackDimensions(String xy) {

		// loop through the series
		System.out.println("seriesCount = " + seriesCount);
		for(int s = 0; s < seriesCount; s++) {
			
			// get name of the image
			reader.setSeries(s);
			String imageName = (String) reader.getSeriesMetadataValue("Image name");
			
			// if this is the first image in a sequence, get dimensions
			// FIXME in theory, we could just get these from the first image, and
			// assume they're all the same throughout a given stack (which
			// seems pretty reasonable)
			if (isSeriesToWrite(xy)) {
				System.out.println("including " + imageName);
				// if this is the first image included, get XYZ dimensions
				if (sizeXYZT[3] == 0) {
					sizeXYZT[0] = reader.getSizeX();
					sizeXYZT[1] = reader.getSizeY();
					sizeXYZT[2] = reader.getSizeZ();
				}
				
				// increment time
				sizeXYZT[3]++;	
			}	
		}		
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
		setEmbryoStackDimensions(xy);
		System.out.println("size: X=" + sizeXYZT[0] + " Y=" + sizeXYZT[1] +
				" Z=" + sizeXYZT[2] + " T=" + sizeXYZT[3]);
		
		// allocate the image buffer (this assumes one byte per pixel)
		imgs = new byte[ 2 * sizeXYZT[2] ][];
		for(int z=0; z < imgs.length; z++)
			imgs[z] = new byte[ sizeXYZT[0] * sizeXYZT[1] ];
		
		// create directory, and create a writer
		new File(outputPath1 + "/dats").mkdirs();
	
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
			String imageName = (String) reader.getSeriesMetadataValue("Image name");
			
			// if this is part of the appropriate sequence, write it
			// ("names of things to skip" is similar to what's in
			// LeicaTifRename.pl)
			if (isSeriesToWrite(xy)) {
//				System.out.println("writing image " + imageName + "   XY = " + getXY(imageName));
				System.out.println(new Date().toString() + " writing image " + imageName);
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

		// seek to that series, and get some dimensions
		reader.setSeries(series);
		int sizeZ = reader.getSizeZ();
		int maxZ = sizeZ - 1;    // this is used for reversing the z-stacks
//		System.out.println("planesPerTimeSlice = " + planesPerTimeSlice);
		
		// loop through the images in the series
		// because the order in which to write these may be different,
		// we store the images temporarily
		for(int i = 0; i < reader.getImageCount(); i++) {
			
			// get channel and plane info (currently swapped,
			// for compatibility with current StarryNite)
			int channel = 1 - meta.getPlaneTheC(series, i).getValue();
			
			// note that this reverses the z-stacks
			int plane = maxZ - meta.getPlaneTheZ(series, i).getValue();
			
			// 0-based index (in this timepoint) of the image to write out
			int imageIndex = reader.getIndex(plane, channel, 0);
			
			// this seems to write things in the right order, but using the
			// getIndex() instead
			// int imageIndex = 2 * plane + channel;		
			
			// read this image
			System.out.println("reading plane with imageIndex = " + imageIndex);
			reader.openBytes(i, imgs[ imageIndex ]);
		}
		
		// actually write out the planes
		String oneTimepointStackFilename =
				outputPath + File.separator + seriesName + "_t" + t + ".ome.tif";
	
		openImageWriterFile(oneTimepointStackFilename, sizeZ);
		System.out.println("opened " + oneTimepointStackFilename);
		for(int p=0; p<imgs.length; p++) {
			writer.saveBytes(p, imgs[p]);
			System.out.println("saved plane " + p);
		}
		writer.close();
	}
	
	/** Opens a file using the ImageWriter; this should be called before
	 * starting to write out images.
	 * @param filename the filename
	 * @param numZ number of z stacks */
	protected void openImageWriterFile(String filename, int numZ)
			throws IOException, FormatException, ServiceException {
	
		// this assumes the reader is already open
		int sizeX = reader.getSizeX();
		int sizeY = reader.getSizeY();
		
		// XXX this is somewhat complicated; it's based on 
		// loci-tools' components/bio-formats/utils/MinimumWriter.java
		IMetadata meta = service.createOMEXMLMetadata();
	    MetadataTools.populateMetadata(meta, 0, (String) null, false,
	    		"XYCZT",   // was "XYZCT"
	    	      FormatTools.getPixelTypeString(FormatTools.UINT8),
	    	      sizeX, sizeY, numZ, 2, 1, 1);

	    // FIXME loci_tools v.5.1.x uses a new API here, which I don't
	    // yet understand (workaround for now is to use v.5.0.x)
	    /*
	  	meta.setPixelsPhysicalSizeX(new PositiveFloat(1.0), 0);
	    meta.setPixelsPhysicalSizeY(new PositiveFloat(1.0), 0);
	    */
	    
	    // storing in little-endian order (as that's what the Perl script
	    // does; although it may not matter)
	    meta.setPixelsBinDataBigEndian(Boolean.FALSE, 0, 0);
	    
   		writer = new ImageWriter();
   		writer.setMetadataRetrieve(meta);
   		
   		// not sure when this needs to be called, but it seems to have to be
   		// before the file is opened
   		writer.setWriteSequentially(true);

//   		writer.setCompression(TiffWriter.COMPRESSION_LZW);
//   		writer.setCompression(TiffWriter.COMPRESSION_UNCOMPRESSED);
   		
   		writer.setId(filename);
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
				System.err.println("Imports Leica images, one stack per timepoint");
				System.err.println("Usage: leica_import_LSMtime.pl file.lif");
				System.err.println("Writes series in the same directory");
				System.exit(1);
			}
			String lifFile = args[0];
		}
		String lifFile = "/var/tmp/20151027_tlp-1_enh4K.lif";

		System.out.println("converting " + lifFile + " to one-stack-per-timepoint format");
    	LeicaImporterLSM li = new LeicaImporterLSM(lifFile);

    	li.getImageNames();
      	li.writeAllSeries();
      	System.out.println("finished at " + new Date().toString());
	}
}
