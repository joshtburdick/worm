package worm.image3d;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.applet.*;
import java.util.*;
import javax.swing.*;
import javax.media.j3d.*;

import java.net.*;
import javax.vecmath.*;

import com.sun.j3d.utils.applet.MainFrame; 
import com.sun.j3d.utils.geometry.ColorCube;
import com.sun.j3d.utils.universe.*;
import com.sun.j3d.utils.behaviors.mouse.*;

import com.sun.j3d.utils.image.TextureLoader;

import worm.image.*;
import worm.image3d.*;

/** A text label which rotates to face the viewer. */
public class Label {

	/** The location of this. */
	private Vector3f location = new Vector3f(0,0,0);
	
	/** The string to display. */
	private String text = "";
	
	/** The BranchGroup containing this. */
	private BranchGroup bg = new BranchGroup();
	
	/** The group node containing the plane. */
	private TransformGroup tg;
	
	/** The Shape node. */
	private OrientedShape3D shape;
	
	/** Image on which to draw the string. */
	private BufferedImage img;
	
	/** Buffer where this will actually be drawn. */
	private ImageComponent2D ic;
	
	/** Constructor.
	 * 
	 * @param size  size of the long edge of the label.
	 */
	public Label(float size) {
		tg = new TransformGroup();
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
	    
	    // the image to draw on
		img = new BufferedImage(256, 32, BufferedImage.TYPE_INT_ARGB);
	    TextureLoader tl = new TextureLoader(img);
	    Texture t = tl.getTexture();

		// create a plane
		QuadArray qa = new QuadArray(4,
				QuadArray.COORDINATES | GeometryArray.TEXTURE_COORDINATE_2);
		float h = (float) img.getHeight() / (float) img.getWidth();
		qa.setCoordinate(0, new Point3f(0,0,0));
		qa.setCoordinate(1, new Point3f(size,0,0));
		qa.setCoordinate(2, new Point3f(size,size*h,0));
		qa.setCoordinate(3, new Point3f(0,size*h,0));
		
		// set texture coordinates
		qa.setTextureCoordinate(0, 0, new TexCoord2f(0.0f,0.0f));
		qa.setTextureCoordinate(0, 1, new TexCoord2f(1.0f,0.0f)); 
	    qa.setTextureCoordinate(0, 2, new TexCoord2f(1.0f,1.0f));
	    qa.setTextureCoordinate(0, 3, new TexCoord2f(0.0f,1.0f));

	    // trying to make this writeable; may not be necessary
	    ic = (ImageComponent2D) t.getImage(0);
//	    ic.setCapability(ImageComponent.ALLOW_IMAGE_WRITE);
	    
 	   	Appearance a = new Appearance();
		TextureAttributes ta = new TextureAttributes();
		a.setTextureAttributes(ta);
		a.setTexture(t);
		
		// make this semi-transparent?
		TransparencyAttributes tra = new TransparencyAttributes();
 		tra.setTransparencyMode (tra.BLENDED);
 		tra.setTransparency (0.5f);
 		a.setTransparencyAttributes (tra);
 		
 		// create the shape
		shape = new OrientedShape3D(qa, a, OrientedShape3D.ROTATE_ABOUT_POINT, new Point3f(0,0,0));
		tg.addChild(shape);
	}
	
	/** Modifies the location and text of this label. */
	public void setLocationAndText(Vector3f location, String text) {
		this.location = location;
		this.text = text;
		
		Transform3D tr = new Transform3D();
		tg.getTransform(tr);
		tr.setTranslation(location);
		tg.setTransform(tr);
			
		// draw the text
		Graphics2D g = (Graphics2D) img.getGraphics();
		g.setColor(new Color(0,0,0,0));
		// g.setColor(new Color(0.3f,0.4f,0.5f,0.5f));  // for debugging
		g.fillRect(0, 0, img.getWidth(), img.getHeight());
		g.setColor(new Color(1,1,1,0.99f));
		g.setFont(new Font("SansSerif", Font.PLAIN, 24));
		g.drawString(text, 0, 24);
		
//		ic.set(img);
	}
	
	public Node getNode() {
		return tg;
	}
}
