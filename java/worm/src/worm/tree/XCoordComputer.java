package worm.tree;

import java.util.HashMap;

/** Computes x coordinates for the tree.
 * Specifically, for each cell, this computes a
 * three-element array, containing
 * - the minimum x-coordinate for this subtree,
 * - this cell's x-coordinate, and
 * - the maximum x-coordinate for this subtree 
 * ??? this may end up being just this cell's x-coordinate */
class XCoordComputer {

	private Cell root;
	
	/** Maximum time to include. */
	private double maxTime;
	
	/** Maximum x coordinate so far. */
	private double maxX = 0;
	
	/** The mapping from Cell to x coordinate. */
	private HashMap<Cell, double[]> cellX
		= new HashMap<Cell, double[]>();
	
	/** Constructor.
	 * 
	 * @param root  the root of the tree
	 */
	XCoordComputer(Cell root) {
		this.root = root; this.maxTime = 1000;
	}
	
	XCoordComputer(Cell root, double maxTime) {
		this.root = root; this.maxTime = maxTime;
	}
	
	/** Computes the mapping from Cell to x coordinate. */
	public HashMap<Cell, double[]> getCellX() {
		updateCellX(root);
		return cellX;
	}
	
	/** Computes mapping from Cell to x coordinate. */
	private void updateCellX(Cell c) {
//		System.out.println(c.name + " " + c.t);

		// if the cell is after the max. time, skip it
		if (c.nucleus.firstKey() > maxTime)
			return;
		
		// if this is a leaf, or we've reached maxTime, allocate a column
		if (c.children.isEmpty() || c.nucleus.lastKey() >= maxTime) {
			// System.out.println(c.name + " " + maxX);
			double a = maxX++;
			double[] x = {a, a, a};
			cellX.put(c, x);
		}
		// otherwise...
		else {
			// do preorder traversal of children
			for(Cell child : c.children)
				updateCellX(child);
		
			// compute min, mean, and max of actual locations
			/* ??? this isn't working, but theoretically should
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;
			for(Cell child : c.children) {
				double x = cellX.get(child)[1];
				if (x < min) x = min;
				if (x > max) x = max;
			}
			*/
			
			// assumes children are sorted. XXX maybe not always true
			double min = cellX.get(c.children.elementAt(0))[0];
			double max = cellX.get(c.children.elementAt(1))[1];
			double[] x = { min, (min+max)/2.0, max };
			cellX.put(c, x);
		}
	}
}
