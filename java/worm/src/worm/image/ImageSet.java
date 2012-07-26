package worm.image;

import java.awt.image.BufferedImage;
import java.io.*;

import javax.imageio.*;
import javax.media.jai.*;
import javax.media.j3d.*;
import com.sun.j3d.utils.image.*;

/** A set of images. */
public class ImageSet {

	/** Reads an image.
	 * FIXME: currently this doesn't read compressed TIFFs :( */
	private BufferedImage readImage(File f) throws IOException {
		PlanarImage planarImage = JAI.create("fileload", f.getAbsolutePath());
	    BufferedImage img = planarImage.getAsBufferedImage();
		System.out.println("img = " + img);
		return img;
	}
	
	/** Constructs a mipmapped Texture2D of a given image. */
	private Texture readTexture2D(File f) throws IOException {		
		BufferedImage img = readImage(f);
		
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
				System.out.println(t);
				
				
			}
		}
		

	}

	
	
	
}
