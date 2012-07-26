package worm.image;

import java.io.*;

public class ImageViewer {

	/** Views a set of images.
	 * 
	 * @param args Currently, just the pathname of the files to read.
	 */
	public static void main(String[] args) {

		try {
			File imagePath = new File(args[0]);
			ImageSet is = new ImageSet(imagePath);
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

}
