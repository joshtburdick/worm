package worm.image;

import java.util.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.regex.*;

import javax.imageio.*;
import javax.media.jai.*;
import javax.media.j3d.*;
import com.sun.j3d.utils.image.*;

/** A set of images. */
public class ImageSet {

	/** The images, keyed first by time, then z. */
	TreeMap<Integer, TreeMap<Integer, Texture>> img
		= new TreeMap<Integer, TreeMap<Integer, Texture>>();

	/** Reads an image.
	 * 	(Despite previous comment, this now seems to be reading compressed
	 * TIFFs; need to check this though. Also, it would be good to know
	 * precisely which .jars need to be where for this to work.) */
	private BufferedImage readImage(File f) throws IOException {
		PlanarImage planarImage = JAI.create("fileload", f.getAbsolutePath());
	    BufferedImage img = planarImage.getAsBufferedImage();
		return img;
	}
	
	/** Constructs a mipmapped Texture2D of a given image. */
	private Texture readTexture2D(File f) throws IOException {		
		BufferedImage img = readImage(f);
		System.out.println(img);
		TextureLoader tl = new TextureLoader(img,
			TextureLoader.GENERATE_MIPMAP | TextureLoader.BY_REFERENCE,
			null);
		return tl.getTexture();
	}
	
	/** Constructs this from a directory. */
	public ImageSet(File imageDir) throws IOException {
	
		// loop through contents of directory
		if (imageDir.isDirectory()) {
			for(File f : imageDir.listFiles()) {
				System.out.println(f);
				Texture t = readTexture2D(f);

				int[] tz = parseFilename(f.getName());
				if (tz != null) {
					if (!img.containsKey(tz[0]))
						img.put(tz[0], new TreeMap<Integer, Texture>());
					img.get(tz[0]).put(tz[1], t);
				}
				else {
					System.err.println("failed to parse file " + f.getName());
				}
			}
		}
	}

	/** Pattern used for extracting numbers from filenames. */
    Pattern filenameParse =
    	Pattern.compile("-t([0-9][0-9][0-9])-p([0-9][0-9])\\.(tif|png)$");

	/** We get which image this is from the filename.
	 * This may need to be generalized at some point. */
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
}
