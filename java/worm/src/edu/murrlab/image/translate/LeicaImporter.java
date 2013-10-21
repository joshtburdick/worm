package edu.murrlab.image.translate;

import java.util.*;

import java.awt.image.*;
import java.io.*;
import java.util.zip.*;
import java.util.regex.*;
import javax.media.jai.*;
import com.sun.media.jai.codec.*;

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

/** Imports a .lif file (hopefully) into format that
 * AceTree and StarryNite can deal with.
 * 
 * Heavily based on the GetPhysicalMetadata example.
 * @author jburdick
 *
 */
public class LeicaImporter {

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
	String[] channelDir = {"tif", "tifR"};
	
	/** The "job ID", set by getImageNames(). */
	// String jobID = null;   probably not used
	
	/** The XY position strings, set by getImageNames(). */
	Set<String> xyPosition = new TreeSet<String>();
	
    private ServiceFactory factory;

    private OMEXMLService service;

    /** Buffer for images (which we'll reuse.) */
    private byte[] img;
    
    /** Constructor. */
	public LeicaImporter(String lifFile)
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
			
			// increment count of "images with this job ID"
			/*
			String j = getJobID(imageName);
			Integer n = count.get(j);
			if (n != null)
				count.put(j, n.intValue() + 1);
			else
				count.put(j, 1);
			*/
			
			// record the XY position (used for separating out
			// images from different embryos)
			String xy = getXY(imageName);
			if (xy != null)
				xyPosition.add( xy );
		}
			
		// find largest sequence name
		/*
		int max = 0;
		for(Integer c : count.values())
			if (c > max)
				max = c;
		for(String k : count.keySet())
			if (count.get(k) == max) {
				this.jobID = k;
				return;
			}
		*/
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
		
		// create directories
		new File(outputPath1 + "/dats").mkdirs();
		new File(outputPath1 + "/tif").mkdirs();
		new File(outputPath1 + "/tifR").mkdirs();
		
		// we copy the metadata into each embryo directory, even though it
		// seems to be the same for each series
		writeMetadata(outputPath1 + "/dats/");
		
		// file to write a list of images in this series
		PrintWriter imageList = new PrintWriter(outputPath1 + "/dats/imageList.tsv");
		imageList.println("time\timage");
		
		// counter for current time point
		int t = 1;
		
		// loop through the series
		for(int s = 0; s < seriesCount; s++) {
			
			// get name of the image
			reader.setSeries(s);
			Hashtable<String, Object> h = reader.getSeriesMetadata();
			String imageName = (String) h.get("Image name");
			
			// if this is part of the sequence, write it
			if (!imageName.contains("Sequence/Image") &&
					!imageName.contains("Sequence/AF") &&
					!imageName.contains("AF Job") &&
					!imageName.contains("DriftAF") &&
//					imageName.startsWith(jobID) &&
					xy.equals( getXY(imageName) )) {
				System.out.println("writing image " + imageName + "   XY = " + getXY(imageName));
				imageList.println(t + "\t" + imageName);
				
				writeSeries(outputPath1, embryoNumber, s, t);
				
				t++;
			}	
		}
		
		imageList.close();
	}

	/** Writes one series (which is actually only one timepoint.)
	 * @param outputPath  directory in which to write the file
	 * @param embryoNumber  which embryo this is (so as to add
	 *   "_L1", "_L2", "_L3", etc. to the name)
	 * @param series  the index of the series in the file
	 * @param t  the timepoint number with which to label this. */
	private void writeSeries(String outputPath, int embryoNumber,
			int series, int t) throws IOException, FormatException, ServiceException {
		reader.setSeries(series);

		for(int i = 0; i < reader.getImageCount(); i++) {
			
			// get channel and plane info
			int channel = meta.getPlaneTheC(series, i).getValue();
			int plane = meta.getPlaneTheZ(series, i).getValue();
			
			// original, presumably slower version of this
			//			byte[] img = reader.openBytes(i);

			// all of this is to avoid re-allocating the byte buffer
			// XXX note that this assumes one byte per pixel
			int numBytes = reader.getSizeX() * reader.getSizeY();
			if (img == null || img.length != numBytes)
				img = new byte[ numBytes ];
			reader.openBytes(i, img);
			
			/*
			System.out.println("series " + series + " image " + i + ": " +
					" C = " + meta.getPlaneTheC(series, i) +
					" T = " + meta.getPlaneTheT(series, i) +
					" Z = " + meta.getPlaneTheZ(series, i));
			*/
			
			String dir = channelDir[ channel ];
			String f = outputPath + "/" + dir + "/" +
					seriesName + "_L" + embryoNumber +
					"-t" + String.format("%03d", t) +
					"-p" + String.format("%02d", plane + 1) + ".tif";
			writeCompressedTiff(img, f);
		}
	}
	
	/** Writes a compressed TIFF file. */
	private void writeCompressedTiff(byte[] data, String filename)
			throws IOException, FormatException, ServiceException {
		
		int sizeX = reader.getSizeX();
		int sizeY = reader.getSizeY();

		/* One way of saving TIFs, just using JDK. However, it
		 * is slow, and doesn't support LZW compression, so
		 * it's currently not used.
		BufferedImage b = new BufferedImage(sizeX, sizeY, BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster wr = b.getRaster();
		int[] data1 = new int[ data.length ];
		for(int i=0; i<data.length; i++)
			data1[i] = data[i];
		wr.setPixels(0,  0,  sizeX,  sizeY, data1);

		// specifying compression
		TIFFEncodeParam param = new TIFFEncodeParam();
		// FIXME: COMPRESSION_LZW isn't supported
		param.setCompression(TIFFEncodeParam.COMPRESSION_DEFLATE);
		
		// this is what actually stores the file
	    RenderedOp op = JAI.create("filestore", b,
                filename, "TIFF", param);
        */
		
		// XXX this is somewhat complicated; it's based on 
		// loci-tools' components/bio-formats/utils/MinimumWriter.java
		IMetadata meta = service.createOMEXMLMetadata();
	    MetadataTools.populateMetadata(meta, 0, null, false, "XYZCT",
	    	      FormatTools.getPixelTypeString(FormatTools.UINT8), sizeX, sizeY, 1, 1, 1, 1);
	    meta.setPixelsPhysicalSizeX(new PositiveFloat(1.0), 0);
	    meta.setPixelsPhysicalSizeY(new PositiveFloat(1.0), 0);

   		ImageWriter writer = new ImageWriter();
   		writer.setMetadataRetrieve(meta);
   		writer.setCompression(TiffWriter.COMPRESSION_LZW);
		writer.setId(filename);
	    writer.saveBytes(0, data);
	    writer.close();
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
		Pattern p = Pattern.compile(".*/A \\d+_S \\d+_U \\d+_V \\d+_(X \\d_Y \\d)_.*");
		Matcher m = p.matcher(imageName);
		return (m.matches() ? m.group(1) : null);
	}

	/** Test utility to dump image information. */
	private void dumpInfo() {
		for(int series = 0; series < 10; series++) {
			reader.setSeries(series);
			System.out.println("series = " + series);
			
			// XXX I don't actually know the legit dimensions of
			// this, so just looping through them all			
			for(int imageIndex = 0; imageIndex < 10; imageIndex++)
				for(int planeIndex = 0; planeIndex < 10; planeIndex++) {
					// ignore out-of-bounds
					try {
						System.out.println(imageIndex + "," + planeIndex + ": " +
							" C=" + meta.getPlaneTheC(imageIndex,planeIndex).getValue() +
							" T=" + meta.getPlaneTheT(imageIndex,planeIndex).getValue() +
							" Z=" + meta.getPlaneTheZ(imageIndex,planeIndex).getValue() +
							"  pos = " + meta.getPlanePositionX(imageIndex, planeIndex) + "," +
									meta.getPlanePositionY(imageIndex, planeIndex) + "," +
									meta.getPlanePositionZ(imageIndex, planeIndex));
					}
					catch (Exception e) {
					}
				}
		}
	}
	
	
	public static void main(String[] args) throws Exception {
		// XXX print errors to console, or a log file?
		// org.apache.log4j.BasicConfigurator.configure();
		if (args.length != 1) {
			System.err.println("imports Leica images");
			System.err.println("Usage: leica_import.pl file.lif");
			System.err.println("writes series in the same directory");
			System.exit(1);
		}
		String lifFile = args[0];
//		String lifFile = "/media/disk2/jburdick/image/testLIF_20130822_JIM221.lif";
		
		System.out.println("converting " + lifFile);
    	LeicaImporter li = new LeicaImporter(lifFile);
    	// li.dumpInfo();
    	li.getImageNames();
      	li.writeAllSeries();
	}
}
