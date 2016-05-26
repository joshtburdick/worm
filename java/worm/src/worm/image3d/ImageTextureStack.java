package worm.image3d;

// possibly see https://www.java.net//node/671267

import java.applet.Applet;
import java.awt.*;
import java.awt.geom.*;
import java.awt.image.*;

import javax.media.j3d.*;
import javax.vecmath.*;
import java.util.*;

import com.sun.j3d.utils.image.TextureLoader;

import worm.image.*;

/** Java3D object, displaying images as a stack of textures.
 * This maintains its own stack of RGBA image planes, on the
 * assumption that the viewed image will differ somewhat from
 * the image stacks.
 * May eventually provide a buffer, for combining different signals. */
public class ImageTextureStack {
	
	/** The image buffer which we write to, and which will be copied
	 * to the actual Java3D images. */
	private BufferedImage[] img;
	
	/** Group node containing all of the planes. */
	private Group objRoot = new Group();
	
	/** The planes. */
	private Shape3D[][] shapes;
	
	/** The actual ImageComponents which are displayed (copying to these
	 * actually updates the images.) The first index is which direction this
	 * set of planes is in. Currently only includes viewing from one direction;
	 * eventually, may add the others (totalling six) later. 
	 * XXX deprecated; this doesn't seem to be working */
//	private ImageComponent2D ic[][];
	
	/** Creates a Shape3D with given corner endpoints. */
	private Shape3D createShape3D(Point3f p0, Point3f p1, Point3f p2, Point3f p3,
			int width, int height) {
		
		// create a plane
		QuadArray qa = new QuadArray(4,
				QuadArray.COORDINATES | GeometryArray.TEXTURE_COORDINATE_2);
		qa.setCoordinate(0, p0);
		qa.setCoordinate(1, p1);
		qa.setCoordinate(2, p2);
		qa.setCoordinate(3, p3);
		
		// set texture coordinates
		qa.setTextureCoordinate(0, 0, new TexCoord2f(0.0f,0.0f));
		qa.setTextureCoordinate(0, 1, new TexCoord2f(1.0f,0.0f)); 
	    qa.setTextureCoordinate(0, 2, new TexCoord2f(1.0f,1.0f));
	    qa.setTextureCoordinate(0, 3, new TexCoord2f(0.0f,1.0f));

	    BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
//	    System.out.println("img = " + img);
		
		// for testing
	    /*
		Graphics2D g = (Graphics2D) img.getGraphics();
		g.setColor(Color.darkGray);
		g.fillRect(0,  0,  img.getWidth(), img.getHeight());
		*/
	    Texture t = new TextureLoader(img, TextureLoader.ALLOW_NON_POWER_OF_TWO).getTexture();

	    // make this writeable
/* XXX not sure this works */
	    ImageComponent ic = t.getImage(0);
	    ic.setCapability(ImageComponent.ALLOW_IMAGE_WRITE);
	    
 	   	Appearance a = new Appearance();
 //	   	a.setCapability(Appearance.ALLOW_TEXTURE_WRITE);
		TextureAttributes ta = new TextureAttributes();
		a.setTextureAttributes(ta);
		a.setTexture(t);
		
		// make this semi-transparent?		
 		TransparencyAttributes tra = new TransparencyAttributes();
 		tra.setTransparencyMode (tra.NICEST);
 		tra.setTransparency (1f);  // should be 1?
 		
 		// alternative version, from
 		// http://www.java3d.org/appearance.html
 		// ??? I don't quite understand the 0.5f as second arg.
 		TransparencyAttributes t_attr =
 		   new TransparencyAttributes(TransparencyAttributes.BLENDED,0.5f,
 			TransparencyAttributes.BLEND_SRC_ALPHA,
 			TransparencyAttributes.BLEND_ONE);
 		
 		a.setTransparencyAttributes (t_attr);
		return new Shape3D(qa, a);
	}
	
	/** Constructs an ImageTextureStack of given dimensions,
	 * initially blank.
	 * 
	 * @param numVoxels  number of voxels in each dimension
	 * @param resolution  size of each voxel (in micrometers)
	 */
	public ImageTextureStack(int[] numVoxels, float[] resolution) {
		shapes = new Shape3D[6][];
		shapes[0] = new Shape3D[ numVoxels[2] ];
		img = new BufferedImage[ numVoxels[2] ];
		
		// find next higher power of two for image sizes (since
		// textures seem to require that)
		// XXX this hopefully won't be needed; we'll either assume non-power-of-2 textures
		// are working, or let TextureLoader rescale stuff
		int[] textureSize = new int[3];
		float[] s = new float[3];
		for(int i=0; i<3; i++) {
			textureSize[i] = 1 << ((int) Math.ceil( Math.log(numVoxels[i]) / Math.log(2) ));
			s[i] = textureSize[i] * resolution[i];
		}

		// first, create the images in the buffer
		// XXX shouldn't be as big
		for(int i=0; i<numVoxels[2]; i++) {
			img[i] = new BufferedImage(numVoxels[0], numVoxels[1],
    			BufferedImage.TYPE_INT_ARGB);

		}
		
		// create the actual textures (and ImageComponents)
//		ic = new ImageComponent2D[6][];
//		ic[0] = new ImageComponent2D[numVoxels[2]];
		
		for(int i=0; i<numVoxels[2]; i++) {

			// currently coordinates are in micrometers
			float x1 = numVoxels[0] * resolution[0];
			float y1 = numVoxels[1] * resolution[1];
			float z1 = i * resolution[2];

			shapes[0][i] = createShape3D(new Point3f(0,0,z1),
					new Point3f(x1,0,z1),
					new Point3f(x1,y1,z1),
					new Point3f(0,y1,z1),
					numVoxels[0], numVoxels[1]);
			
//			ic[0][i] = (ImageComponent2D) shapes[i].getAppearance().getTexture().getImage(0);
			objRoot.addChild(shapes[0][i]);		
		}
	}
		
	/** Constructs a transparent version of an image.
	 * Possibly not used. */
	private BufferedImage makeTransparent(BufferedImage b) {
		/*
		 ColorModel cm = b.getColorModel();
		 boolean isAlphaPremultiplied = cm.isAlphaPremultiplied();
		 WritableRaster raster = b.copyData(null);
		 return new BufferedImage(cm, raster, isAlphaPremultiplied, null);
		 */
		BufferedImage b1 = new BufferedImage(b.getWidth(), b.getHeight(), BufferedImage.TYPE_INT_ARGB);
		Graphics2D g = (Graphics2D) b1.getGraphics();
		/*
		AlphaComposite composite = AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 1f);
		g.setComposite(composite);		
		g.drawImage(b, new AffineTransform(1f,0f,0f,1f, 0, 0), null);
		*/
		float[] color = {1,0,0,1};   // FIXME
		addColorizing(b.getRaster(), b1.getRaster(), color);
		
		return b1;
	}
	
	/** Adds a grayscale image to a Raster, as some color. */
	public void addColorizing(java.awt.image.Raster from,
			java.awt.image.WritableRaster to,
			float[] color) {
		// ??? why did this have two channels here?
//		float[][] m = {{3f,0},{0,0},{0,0},{1f,0}};
		
		// this superficially seems to be working better
//		float[][] m = {{3f},{0},{0},{1f}};
		float[][] m = {{color[0]}, {color[1]}, {color[2]}, {color[3]}};
		
		BandCombineOp op = new BandCombineOp(m, null);
		op.filter(from, to);
	}
		
	/** Clear buffer. */
	public void clearBuffer() {
		/*
		for(int i=0; i<img.length; i++) {
			Graphics2D g = (Graphics2D) img[i].getGraphics();
			g.setColor(Color.black);
			g.fillRect(0, 0, img[i].getWidth(), img[i].getHeight());
		}
		*/
	}

	/** Adds a greyscale image to the buffer. */
	public void addToBuffer(SortedMap<Integer, BufferedImage> a, float[] color) {
		for(int i=0; i<img.length; i++) {
			/*
			System.out.println(i);
			System.out.print(a.get(i) + " ");
			System.out.println(img[i]);
			*/
			float[][] m = {{color[0]}, {color[1]}, {color[2]}, {color[3]}};
			
			if (a.get(i) != null && img[i] != null) {
				addColorizing(a.get(i).getRaster(), img[i].getRaster(), color);

				// another approach
				// Graphics2D g = (Graphics2D) img[i].getGraphics();
				// g.drawImage(a.get(i), new ColorConvertOp(null), 0, 0);
				// g.drawImage(a.get(i), 0, 0, Color.black, null);
			}
			
			// XXX testing	
/*		
			Graphics2D g = (Graphics2D) img[i].getGraphics();
			g.setColor(Color.darkGray);
			g.fillRect(0, 0, img[i].getWidth(), img[i].getHeight());
*/		
		}
	}
	
	/** Updates the actual Java3D ImageComponents with the current buffer. */
	public void showBuffer() {
	
		// for now, only including the first direction
		Shape3D[] s = shapes[0];
		for(int i=0; i<s.length; i++) {
			ImageComponent2D ic = (ImageComponent2D) s[i].getAppearance().getTexture().getImage(0);
//			ic1[i].setSubImage(img[i], img[i].getWidth(), img[i].getHeight(), 0, 0, 0, 0);
//			System.out.println("image size = " + img[i].getWidth() + " x " + img[i].getHeight());
//			System.out.println("ic size = " + ic.getWidth() + " x " + ic.getHeight());
			// ic1[i].setSubImage(img[i], 100, 100, 0, 0, 0, 0);
			ic.set(img[i]);
		}
		
	}
	
	public void showBufferSlow() {
		/*
		// for now, only including the first direction
		Shape3D[] s = shapes[0];
		for(int i=0; i<s.length; i++) {
			Appearance a = s[i].getAppearance();
			Texture t = new TextureLoader(img[i], TextureLoader.GENERATE_MIPMAP).getTexture();
			a.setTexture(t);
		}
		*/

	}
	
	/** Seriously deprecated. */
	/*
	public void showBufferOld() {
		// create buffer to rescale to the size of the texture map
		Texture t = shapes[1].getAppearance().getTexture();
		BufferedImage scaledBuf = new BufferedImage(t.getWidth(), t.getHeight(),
				BufferedImage.TYPE_4BYTE_ABGR);
		
		for(int i=0; i<img.length; i++) {
			System.out.println(i);
			if (img[i] != null & shapes[i] != null) {
//				addToBuffer(a.get(i).getRaster(), img[i].getRaster(), color);
				
				// scale the buffer; Make image always std_height tall
				AffineTransform tx = new AffineTransform();
				tx.scale((float) scaledBuf.getWidth() / (float) img[i].getWidth(),
						(float) scaledBuf.getHeight() / (float) img[i].getHeight());
				AffineTransformOp op = new AffineTransformOp(tx, AffineTransformOp.TYPE_BILINEAR);
				op.filter(img[i], scaledBuf);
				
				// then copy scaled image to the actual texture
				ImageComponent2D ic = (ImageComponent2D) shapes[i].getAppearance().getTexture().getImage(0);
				ic.set(scaledBuf);
			}
			
		}
	}
	*/
	/** Accesses the group. */
	public Node get() {
		return objRoot;
	}
	
	/** XXX test of drawing onto the images after it's already been created. */
	/*
	public void testDrawing() {
	
		for(int i=1; i<30; i++) {
			Texture t = shapes[i].getAppearance().getTexture();
			ImageComponent2D ic = (ImageComponent2D) t.getImage(0);
			Graphics2D g = (Graphics2D) ic.getImage().getGraphics();
			g.setColor(Color.CYAN);
			g.drawRect(1, 1, 100, 100);
		}
		System.out.println("drew on textures");
	}
	*/
}
