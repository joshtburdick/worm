package worm.image;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.util.regex.*;
import java.awt.geom.*;
import java.awt.image.*;

import javax.imageio.*;
import javax.media.jai.*;
import javax.media.j3d.*;
import com.sun.j3d.utils.image.*;

/** A set of images. This basically just reads the images in. */
public class ImageSet {

	/** Base location of images. This includes the pathname
	 * (if any), and the series name. */
	public String baseName;
	
	/** The images, keyed first by time, then z. */
	public TreeMap<Integer, TreeMap<Integer, BufferedImage>> img
	= new TreeMap<Integer, TreeMap<Integer, BufferedImage>>();
	
	/** The textures, keyed first by time, then z.
	 * XXX this is essentially a cache, which may not be used.
	 * XXX also, this is somewhat deprecated, as we're also
	 * storing the images. */
	public TreeMap<Integer, TreeMap<Integer, Texture>> tex
		= new TreeMap<Integer, TreeMap<Integer, Texture>>();

	/** Voxel size, in micrometers (for x, y, and z.) 
	 * FIXME this shouldn't be hardcoded */
	public double resolution[] = { 0.087, 0.087, 0.504 };
	
	/** Constructor, which doesn't read in any images. */
	public ImageSet(String baseName) {
		this.baseName = baseName;
	}
	
	/** Reads an image.
	 * XXX It would be good to know precisely which .jars need to be
	 * where for this to work.) */
	private BufferedImage readImage(File f) throws IOException {
		PlanarImage planarImage = JAI.create("fileload", f.getAbsolutePath());
	    BufferedImage img = planarImage.getAsBufferedImage();
		return img;
	}
	
	/** Constructs a (possibly mipmapped) Texture2D of a given image. 
	 * Deprecated. */
	private Texture readTexture2D(File f) throws IOException {		
		BufferedImage i = readImage(f);
		System.out.println(i);
		
		// for now, not including mipmaps (though they may be useful later)
		TextureLoader tl = new TextureLoader(i,
			TextureLoader.BY_REFERENCE,  //  | TextureLoader.GENERATE_MIPMAP,
			null);
		return tl.getTexture();
	}
	
	/** Constructs this from a directory. 
	 * Deprecated (reading all the images at once probably isn't a good plan.) */
	public ImageSet(File imageDir) throws IOException {
	
		// loop through contents of directory
		if (imageDir.isDirectory()) {
			for(File f : imageDir.listFiles()) {
				System.out.println(f);
				BufferedImage i = readImage(f);
				Texture t = readTexture2D(f);

				int[] tz = parseFilename(f.getName());
				if (tz != null) {
					
					if (!img.containsKey(tz[0]))
						img.put(tz[0], new TreeMap<Integer, BufferedImage>());
					img.get(tz[0]).put(tz[1], i);			
							
					if (!tex.containsKey(tz[0]))
						tex.put(tz[0], new TreeMap<Integer, Texture>());
					tex.get(tz[0]).put(tz[1], t);
				}
				else {
					System.err.println("failed to parse filename " + f.getName());
				}
			}
		}
	}

	/** Pattern used for extracting numbers from filenames. */
    Pattern filenameParse =
    	Pattern.compile("-t([0-9][0-9][0-9])-p([0-9][0-9])\\.(tif|png)$");

	/** Gets which image this is from the filename.
	 * This may need to be generalized at some point, and may not be used.
	 * @param filename  an image filename
	 * @return the time and plane. 
	 * Possibly deprecated. */
	private int[] parseFilename(String filename) {
		int[] r = new int[2];
		Matcher m = filenameParse.matcher(filename);
		if (m.find()) {
			r[0] = Integer.parseInt(m.group(1));
			r[1] = Integer.parseInt(m.group(2));
			return r;
		}
		else
			return null;
	}
	
	/** Constructs the filename for a given time and z. */
	private String imageFileName(int t, int z) {
		return baseName + String.format("-t%03d-p%02d.tif", t, z);
	}
	
	/** Gets all images at a particular time point.
	 * (Later this may include caching.)
	 * @param t  the timepoint of images to fetch
	 * @return  a hash, keyed by plane number, of the images. */
	public SortedMap<Integer, BufferedImage> getImagesAtTime(int t) {
		TreeMap<Integer, BufferedImage> images = new TreeMap<Integer, BufferedImage>();
System.err.println("reading images from " + baseName + ", time " + t);		
		// this reads in images until it fails to read in a file,
		// and then assumes there aren't any more images
		for(int z=1; z<=99; z++) {
			File imageFile = new File(imageFileName(t, z));
//			System.out.println("trying to read " + imageFile);
			if (!imageFile.exists())
				return images;
			try {
				BufferedImage img = readImage(imageFile);
				images.put(z, img);
			} catch (IOException e) {
				return images;
			}
		}
		
		return images;
	}
}
