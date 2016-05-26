package worm.image3d;

import java.awt.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import com.sun.j3d.utils.geometry.*;

/** Convenience wrapper for a Sphere representing a cell,
 * whose location, color, and visibility can be set. */
public class CellSphere {

	/** The transform. */
	Transform3D transform = new Transform3D();
	
	/** The group node. */
	TransformGroup tg;
	
	/** The appearance. */
	Appearance appearance;
	
	/** The material. */
	Material material;
	
	/** The sphere. */
	Sphere sphere;
	
	/** Constructor. */
	public CellSphere() {
		tg = new TransformGroup();
	    tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
	    tg.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		tg.setTransform(transform);
		
		material = new Material();
		material.setShininess(1f);
		material.setAmbientColor(0,0,0);
		material.setDiffuseColor(0,0,0);
		material.setEmissiveColor(0f,0.0f,0.0f);
		material.setSpecularColor(1,0,0);
		
		appearance = new Appearance();
		appearance.setMaterial(material);
		
		TransparencyAttributes ta =
			new TransparencyAttributes(TransparencyAttributes.BLENDED,
					0.5f,
					TransparencyAttributes.BLEND_SRC_ALPHA,
					TransparencyAttributes.BLEND_ONE);
		appearance.setTransparencyAttributes(ta);
		
		// make this invisible, initially
//		appearance.getRenderingAttributes().setVisible(false);

		sphere = new Sphere(1f,
				Sphere.GENERATE_NORMALS,
				appearance);
		tg.addChild(sphere);
	}
	
	public Node getNode() {
		return tg;		
	}
	
	/** Sets location. */
	public void setLocation(Vector3f x, Vector3d size) {
		transform.setTranslation(x);
		transform.setScale(size);
		tg.setTransform(transform);
	}
	
	/** Sets location (and sets size to a sphere.) */
	public void setLocation(Vector3f x, double size) {
		transform.setTranslation(x);
		Vector3d s = new Vector3d(size, size, size);
		transform.setScale(s);
		tg.setTransform(transform);
	}
	
	/** Sets color. */
	public void setColor(Color3f c) {
		appearance.getMaterial().setSpecularColor(c);
	}
	
	/** Sets visibility. */
	public void setVisible(boolean visible) {
		appearance.getRenderingAttributes().setVisible(visible);
	}
	
}
