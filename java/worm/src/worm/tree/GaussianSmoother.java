package worm.tree;

import java.util.*;
import cern.colt.list.DoubleArrayList;
import cern.jet.random.*;
import cern.jet.stat.*;

/** Gaussian smooths an image signal.
 * Um, yeah, actually, this may just take the maximum,
 * or compute an EWMA.
 * Adapted from smoother in EmbryoDB. */
public class GaussianSmoother {

	/** The root of the tree whose expression we're smoothing. */
	private Cell root;
	
	/** Number of points to include. */
	public int windowSize = 5;
	
	/** Normal distribution used for weighting. */
	private Normal f;
	
	/** Alpha weight for EWMA. */
	private double alpha = 0.1;
	
	/** The intensity values being added. */
	private DoubleArrayList x = new DoubleArrayList();
	
	/** The weights for those intensities. */
	private DoubleArrayList w = new DoubleArrayList();

	/** Constructor. */
	public GaussianSmoother(Cell root) {
		this.root = root;
	}
	
	/** Smooths expression.
	 * @param dest  where to write the smoothed data to
	 * @param destIndex  index of the signal in
	 * 	which to write expression
	 */
	public void smoothMaximizing(HashMap<String, Cell> dest, int destIndex) {
		smoothMaximizing(dest, destIndex, root, 0);
	}
	
	private void smoothMaximizing(HashMap<String, Cell> dest, int destIndex,
			Cell c, double maxExprSoFar) {

		// get where we'll be writing data to
		TreeMap<Double, Nucleus> destN = dest.get(c.name).nucleus;
		
		// copy over data at each timepoint
		for(Nucleus n : c.nucleus.values()) {
			
			// ??? is this the right "original" channel to include?
			double unsmoothed = n.expr[3];
			
			// possibly update maximum, and write it
			if (unsmoothed > maxExprSoFar)
				maxExprSoFar = unsmoothed;
			destN.get(n.t).expr[destIndex] = unsmoothed;
		}
		
		// recursive call
		for(Cell ch : c.children)
			smoothMaximizing(dest, destIndex, ch, maxExprSoFar);
	}
	
	/** Smooths expression using an exponentially weighted moving average.
	 * @param dest  where to write the smoothed data to
	 * @param destIndex  index of the signal in
	 * 	which to write expression
	 */
	public void smoothEWMA(HashMap<String, Cell> dest, int destIndex) {
		smoothEWMA(dest, destIndex, root, 0);
	}
	
	private void smoothEWMA(HashMap<String, Cell> dest, int destIndex,
			Cell c, double exprSmoothedSoFar) {
// System.out.println("smoothing cell " + c.name);
		// get where we'll be writing data to
		TreeMap<Double, Nucleus> destN = dest.get(c.name).nucleus;
		
		// copy over data at each timepoint
		for(Nucleus n : c.nucleus.values()) {
			
			// ??? is this the right "original" channel to include?
			double unsmoothed = n.expr[3];
			
			// update smoothed expression
			exprSmoothedSoFar = alpha * unsmoothed + (1-alpha) * exprSmoothedSoFar;
			
			// possibly update maximum, and write it
//			if (unsmoothed > maxExprSoFar)
//				maxExprSoFar = unsmoothed;
			destN.get(n.t).expr[destIndex] = exprSmoothedSoFar;
		}
		
		// recursive call
		for(Cell ch : c.children)
			smoothEWMA(dest, destIndex, ch, exprSmoothedSoFar);
	}

}
