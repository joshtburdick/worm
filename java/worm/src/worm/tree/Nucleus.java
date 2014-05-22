package worm.tree;

import java.util.*;

/** The expression in a cell, over some range of time. */
public class Nucleus {
	
	/** Unique integer index for this cell (possibly used if
	 * expression gets represented by a matrix.) */
	int index;

	/** The coordinates for this cell, and its diameter.
	 * By convention, the expression data is only valid
	 * if size > 0. */
	public double x, y, z, size, t;
	
	/** The expression data. */
	public double[] expr;
	
	/** The children of this cell, if it divides at this point. */
	public Vector<Nucleus> children = new Vector<Nucleus>();

	/** Constructor for a cell with no location or expression data. */
	public Nucleus() {
		this.expr = new double[5];
	}
}
