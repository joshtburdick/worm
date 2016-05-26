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

import worm.image.*;
import worm.tree.*;

public class ImagePlaneViewer extends Applet
implements ActionListener {

	/** Base pathname of files (including images and cell locations.) */
	public String basePathname = "";
	
	/** Where images are. */
	// public String imagePathname = "/home/jburdick/gcb/data/image/tifs/green_small/";
	
	/** Manages reading images from files. */
	public ImageSet redImages, greenImages;
	
	/** The cell locations. */
    HashMap<String, Cell> cellTree;
    
  	/** Root of the scene graph. */
	BranchGroup objRoot;
	
	/** The Java3D object which actually displays the images. */
	public ImageTextureStack stack;
	
	/** Displays the labels. */
	public CellSpheres cellLoc;

	/** The viewer, etc. */
	SimpleUniverse simpleU;
	
	/** "Rotation" object. */
    TransformGroup objRotate;

    /** Current time. */
    int time;
    
    /** Preferences about brightness, etc. */
    ImagePrefs prefs;

	private void createSceneGraph() {
	
		objRoot = new BranchGroup();
	    objRotate = new TransformGroup();
	    objRotate.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
	    objRotate.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
	    
	    objRoot.addChild(objRotate);
		
		// FIXME dimensions shouldn't be hardwired
		int[] numVoxels = {712, 512, 60};
		float[] resolution = { 0.087f, 0.087f, 0.504f };
		stack = new ImageTextureStack(numVoxels, resolution);
		
		Node a = stack.get();
		objRotate.addChild(a);

		cellLoc = new CellSpheres();
		objRotate.addChild(cellLoc.getNode());
		
		// add a sphere (for testing)
		/*
		CellSphere cs = new CellSphere();
		cs.setLocation(new Vector3f(4,4,4), new Vector3d(0.1,0.1,0.1));
				objRotate.addChild(cs.getNode());
			*/	
		
		addMouseBehaviors(objRotate);

		// Let Java 3D perform optimizations on this scene graph.
	    objRoot.compile();
	}
	
	/** Adds a simple set of mouse behaviors. */
	private void addMouseBehaviors(TransformGroup tg) {
		
        tg.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);	
        BoundingSphere bs = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 1e9);

        // rotation
        MouseRotate myMouseRotate = new MouseRotate();
//        myMouseRotate.setFactor(5);
        myMouseRotate.setTransformGroup(tg);
        myMouseRotate.setSchedulingBounds(bs);
        tg.addChild(myMouseRotate);

        // adding translation
        MouseTranslate myMouseTranslate = new MouseTranslate();
        myMouseTranslate.setFactor(2);
        myMouseTranslate.setTransformGroup(tg);
        myMouseTranslate.setSchedulingBounds(bs);

        tg.addChild(myMouseTranslate);
		
        // ... and zooming
        MouseWheelZoom myMouseWheelZoom = new MouseWheelZoom();
        myMouseWheelZoom.setFactor(3);
        myMouseWheelZoom.setTransformGroup(tg);
        myMouseWheelZoom.setSchedulingBounds(bs);
        tg.addChild(myMouseWheelZoom);
		
	}

	public void actionPerformed(ActionEvent a) {
		if (a.getActionCommand().equals("Update"))
			updateImages();
	}
	
	/** Updates the image stack, possibly changing the time point,
	 * and image settings.
	 * Currently, not very optimized. */
	private void updateImages() {
		System.err.print("in updateImages()...");
		ImagePrefs.Prefs p = prefs.getPrefs();
		if (p == null)
			return;
		SortedMap<Integer, BufferedImage> red = redImages.getImagesAtTime(p.time);
		SortedMap<Integer, BufferedImage> green = greenImages.getImagesAtTime(p.time);

		// get parameters
		/*
		float redBrightness = new Float(prefs.redBrightness.getText()).floatValue();
		float redAlpha = new Float(prefs.redalpha.getText()).floatValue();
		float greenBrightness = new Float(prefs.greenBrightness.getText()).floatValue();
		float greenAlpha = new Float(prefs.greenAlpha.getText()).floatValue();
*/
		// clear the buffer to be displayed
		stack.clearBuffer();
		
		float[] cRed = {p.redBrightness,0f,0f,p.redAlpha};
		stack.addToBuffer(red, cRed);
		float[] cGreen = {0f,p.greenBrightness,0f,p.greenAlpha};
		stack.addToBuffer(green, cGreen);
		
		// having recomputed the images, "flip" the buffer to show the updated stuff
		stack.showBuffer();
		
		// update label locations (if they're known)
		if (cellTree != null)
			cellLoc.setLocations(cellTree, p.time);
		
		System.err.println("finished update");
	}
	
	/** Constructor.
	 * 
	 * @param imagePathname  base 
	 */
	public ImagePlaneViewer(File imagePath) throws FileNotFoundException, IOException {
		
		String baseName = imagePath.getName();
		
		// set up to read the images (but don't read them yet)
	    redImages = new ImageSet(imagePath + "/tifR/" + baseName);
        greenImages = new ImageSet(imagePath + "/tif/" + baseName);

        // read the image annotation (if it's not present, leave this null,
        // as we'd still like to be able to view images in that case)
        String annotFile = imagePath + "/dats/CD" + baseName + ".csv";
        try {
        	CDFileParser parser =
        		new CDFileParser(new FileReader(annotFile));
        	cellTree = parser.parseFile();
        }
        catch (Exception e) {
        	System.err.println("can't open annotation file " + annotFile);
        }
        
        setLayout(new BorderLayout());
        GraphicsConfiguration config =
           SimpleUniverse.getPreferredConfiguration();

        Canvas3D canvas3D = new Canvas3D(config);
        add("Center", canvas3D);

        prefs = new ImagePrefs();
        add("East", prefs);
        prefs.updateButton.addActionListener(this);
        
        createSceneGraph();

        // SimpleUniverse is a Convenience Utility class
        simpleU = new SimpleUniverse(canvas3D);
        
        View v = simpleU.getViewer().getView();
        System.out.println("front clip distance = " + v.getFrontClipDistance());
        System.out.println("back clip distance = " + v.getBackClipDistance());
        // FIXME: update these distances?
//        v.setFrontClipDistance(0.1);
        v.setBackClipDistance(100);
        
        // move the ViewPlatform back a bit
        simpleU.getViewingPlatform().setNominalViewingTransform();     
        
        simpleU.addBranchGraph(objRoot);
	}
	
	/** Views a stack of image planes.
	 * @param args
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException {
		// if no image set is given, use a default
		String baseDir = "/gpfs/fs0/u/twal/images//20120426_JIM122_L2/";
//		String baseDir = "/var/tmp/data/image/20130529_nhr-67_JIM205_L3/";
		if (args.length == 2)
			baseDir = args[1];
		
		ImagePlaneViewer ipv = new ImagePlaneViewer(new File(baseDir));
        Frame frame = new MainFrame(ipv, 400, 300);
        frame.setVisible(true);
        ipv.updateImages();
	}

}
