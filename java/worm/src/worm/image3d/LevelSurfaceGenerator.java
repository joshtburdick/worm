package worm.image3d;

import java.awt.image.*;
import java.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import com.sun.j3d.utils.geometry.*;

/** Generates a level surface.
 * Currently, this only implements Marching Tetrahedrons, and isn't
 * optimized. */
public class LevelSurfaceGenerator {

	/** The images in which to find a level surface. */
	TreeMap<Integer, BufferedImage> img;
	
	/** The level to find. */
	double level;
	
	public LevelSurfaceGenerator(TreeMap<Integer, BufferedImage> img, double level) {
		this.img = img;
		this.level = level;
	}
	
	/** Interpolates a point. Possibly not used.
	 * 
	 * @param x  the value of x
	 * @param f  the value of f(x), which should be > level
	 * @param f1  the value of f(x+1), which should be <= level
	 * @return  interpolated value between x and x+1
	 */
	private double interpolate(int x, double f, double f1) {
		return x + f / (f - f1);
	}
	
	/** Constructs a mesh with some resolution.
	 * 
	 * @param img  the stack of images to construct a level set from
	 * @param level  the level of the surface
	 * @return  a geometry representing that level set
	 */
	public Geometry constructMesh(TreeMap<Double, BufferedImage> img) {
	
		// compute mesh, by looping through pairs of adjacent images
		Vector<Slice> slices = new Vector<Slice>();
		int totalVertices = 0;
		for(int z = 1; z<3; z++) {    // XXX just doing a few slices
			System.out.println("doing slice " + z);
			Slice s = new Slice(img.get(z), img.get(z+1), z);
			totalVertices += s.getVertices().size();
		}
		
		// place in a triangle array. Note that we're relying on various
		// Java3D utilities to reduce the number of vertices (which may
		// be slow.)
		TriangleArray g =
				new TriangleArray(totalVertices, GeometryArray.COORDINATES);
		int i = 0;
		for(Slice s : slices) {
			Vector<Point3f> v = s.getVertices();
			g.setCoordinates(i, (Point3f[]) v.toArray());
			i += v.size();
		}
		
		// add normals
		GeometryInfo gi = new GeometryInfo(g);
		NormalGenerator ng = new NormalGenerator();
		ng.generateNormals(gi);
		
		return gi.getGeometryArray();
	}
	
	/** Utility class to create the mesh for two planes.
	 * This is intended to avoid excessive dereferencing of the
	 * table of images (and potentially simplify parallelizing.) */
	class Slice {
	
		/** The Rasters of the two images we're currently dealing with. */
		java.awt.image.Raster img1, img2;
		
		/** Current x, y, and z coordinates. z is the coordinate of
		 * img1; img2 is at z+1. */
		int x, y, z;
		
		/** The level at which to slice. */
		float level;
		
		/** Values at the corners of the current cube (cached for speed.)
		 * Numbering is as in Paul Bourke's code, linked from
		 * http://paulbourke.net/geometry/polygonise/ */
		float[] f = new float[8];
		
		/** Vertices of triangles constructed in this slice. This should
		 * end up having a length divisible by three. */
		Vector<Point3f> vertices = new Vector<Point3f>();
		
		/** Constructor. */
		public Slice(BufferedImage img1, BufferedImage img2, int z) {
			this.img1 = img1.getRaster();
			this.img2 = img2.getRaster();
			this.z = z;
		}

		/** Scans through the images, adding parts where the surface is. */
		public void scanImages() {
			for(x=0; x<img1.getWidth()-1; x++)
				for(y=0; y<img1.getHeight()-1; y++) {
					checkTetrahedron(0,1,2,3);   // FIXME
					
				}
		}
		
		/** This should only be called after geometry is computed. */
		public Vector<Point3f> getVertices() {
			return vertices;
		}
		
		/** Checks one tetrahedron. */
		private void checkTetrahedron(int v0, int v1, int v2, int v3) {
		
			// summarize which vertices are less than the level
			int triIndex = 0;
			if (f[v0] < level) triIndex |= 1;
			if (f[v1] < level) triIndex |= 2;
			if (f[v2] < level) triIndex |= 4;
			if (f[v3] < level) triIndex |= 8;

			// loop through the cases
			switch (triIndex) {
			case 0x00:
			case 0x0f:
				// in these cases, all vertices are inside or outside,
				// so we don't add anything to the surface
				break;
			case 0x01:
				addTriangle(v1, v2, v3);
				

				break;
			}
			
		}
		
		/** Fills in the points at the corners. */
		private void getFunctionAtCorners() {
			// XXX this is arguably a bit inefficient
			f[0] = img1.getSampleFloat(x, y, 0);
			f[1] = img1.getSampleFloat(x+1, y, 0);
			f[2] = img1.getSampleFloat(x+1, y+1, 0);
			f[3] = img1.getSampleFloat(x, y+1, 0);
			
			f[4] = img2.getSampleFloat(x, y, 0);
			f[5] = img2.getSampleFloat(x+1, y, 0);
			f[6] = img2.getSampleFloat(x+1, y+1, 0);
			f[7] = img2.getSampleFloat(x, y+1, 0);
		}
		
		/** Adds a triangle, with corners at midpoints between the given
		 * vertices.
		 */
		private void addTriangle(int v1, int v2, int v3) {
			Point3f p1 = getPoint(v1), p2 = getPoint(v2), p3 = getPoint(v3);
			
			// add each midpoint
			Point3f w = new Point3f();
			w.interpolate(p1, p2, 1);
			vertices.add(w);
			w.interpolate(p2, p3, 1);
			vertices.add(w);
			w.interpolate(p1, p3, 1);
			vertices.add(w);
		}
		
		/** Gets a single point, based on its index in the above numbering. */
		private Point3f getPoint(int i) {
			return new Point3f( (i==1 || i==2 || i==5 || i==6) ? x+1 : x,
					(i == 2 || i==3 || i==6 || i==7) ? y+1 : y,
					i >= 4 ? z+1 : z);
		}

	}
	
}
