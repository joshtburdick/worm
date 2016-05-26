package worm.image3d;

import java.awt.*;
import java.util.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import worm.tree.*;

/** A set of spheres and/or labels (pre-allocated to hopefully reduce
 * memory usage.) */
class CellSpheres {

	/** The group node. */
	TransformGroup tg;
	
	/** The individual spheres (not currently used.) */
	Vector<CellSphere> sphere = new Vector<CellSphere>(1000);
	
	/** ... and labels. */
	Vector<Label> label = new Vector<Label>(1000);
	
	/** Constructor. */
	public CellSpheres() {
		/*
		for(int i=0; i<sphere.size(); i++) {
			CellSphere s = new CellSphere();
			tg.addChild(s.getNode());
		}
		*/
		for(int i=0; i<label.size(); i++) {
			Label l = new Label(100f);
			tg.addChild(l.getNode());
		}
	}
	
	/** Returns the node to add to the scene graph. */
	public Node getNode() {
		return tg;		
	}
	
	/** Sets sphere locations to locations in a lineage tree
	 * at a particular time.
	 * @param c  root of the lineage tree
	 * @param time  time point to include */
	public void setLocations(HashMap<String, Cell> tree, double time) {

		// first, make all spheres invisible (??? is this inefficient?)
//		for(CellSphere s : sphere)
//			s.setVisible(false);

		for(Label l : label)
			l.setLocationAndText(new Vector3f(0,0,0), "");
		
		int i = 0;
		for (Cell c : tree.values()) {
			setLocations1(c, time, i);
			if (i >= label.size())
				return;
		}
	}
	
	/** Function which does the recursion for the above. */
	private void setLocations1(Cell c, double time, int i) {

		// is this time in the current cell?
		Nucleus n = c.nucleus.get(time);
		
		// if so, set this label
		if (n != null) {
			Vector3f loc = new Vector3f((float) n.x, (float) n.y, (float) n.z);
			System.out.println("set location and text : " + c.name + " " +
					n.x + " " + n.y + " " + n.z);
			label.get(i).setLocationAndText(loc, c.name);
		}
	}
	
}
