package org.murrlab.image;

import java.util.*;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;

import org.rhwlab.tree.Cell;
import org.rhwlab.tree.CellData;

/** Smooths a channel of an image. */
public class MedianSmoother {

	/** Number of adjacent numbers to include in the median. */
	public int k;
	
	/** Root of the tree. */
	Cell root;
	
	/** List of Cells currently being median-smoothed. */
	private Deque<Double> x;
	
	/** Stores medians as they're being computed. */
	private HashMap<CellData, Double> median = new HashMap<CellData, Double>();
	
	/** Used for computing the median. */
	private DoubleArrayList a;
	
	/** Constructor.
	 * @param root  the root of the lineage. */
	public MedianSmoother(Cell root, int k) {
		this.root = root;
		this.k = k;
		a = new DoubleArrayList();
	}
	
	/** Sets everything to its median. */
	public void computeMedians() {
		walkTree();
		setMedians();
	}
	
	/** Gets the median of the k cells before a given Cell. */
	private double getMedian() {
		int i = 0;
		for(Double x1 : x) {
			a.set(i, x1);    // ??? use setQuick?
			i++;
		}
		a.setSize(i);
		return Descriptive.median(a);
	}
	
	/** Walks the tree, computing median of each number. */
	private void walkTree() {
		for(Enumeration e = root.depthFirstEnumeration(); e.hasMoreElements(); ) {
			Cell c = (Cell) e.nextElement();
			Vector<CellData> v = c.getCellData();
			for(int i=0; i<v.size(); i++) {
				CellData cd = (CellData) v.elementAt(i);
				median.put(cd, getMedian(c, i));
			}
		}
	}
	
	/** Walks the tree, computing lists of k previous numbers. */
	private void walkTree1(Cell c, int i) {
		
//		processCellData(c.getCellData());
		
		// process more data in this node?
		if (i+1 < c.getChildCount())
			walkTree1(c, i+1);
		else {
			// visit child nodes
			for(int j=0; j<c.getChildCount(); j++) {
				Cell ch = (Cell) c.getChildAt(j);
				walkTree1(ch, 0);
			}
		}
	}
	
	/** Deprecated. */
	private void processCellData(CellData cd) {
		Double oldLast = null;
		if (x.size() >= k)
			oldLast = x.removeLast();
		x.addFirst(cd.iNucleus.rwcorr3 + 0.0);
		
		// compute median of current cell
		median.put(cd, getMedian());

		x.removeFirst();
		if (oldLast != null)
			x.addLast(oldLast);
	}
	
	/** Gets the median for a given cell. */
	private double getMedian(Cell c, int i) {
		a.setSize(k);
		Vector<CellData> v = c.getCellData();
		int j;
		for(j = 0; j < k; j++) {
//			System.out.println("c = " + c + "   i = " + i);
			CellData cd = v.elementAt(i);
			a.set(j, cd.iNucleus.rweight);
			i--;
			
			// possibly add in data from the previous cell
			if (i < 0) {
				c = (Cell) c.getParent();
				if (c == null)
					break;
				else {
					v = c.getCellData();
					if (v == null || v.size() == 0)
						break;
					i = v.size() - 1;
				}
			}
		}
		
		a.setSize(j);
		return Descriptive.median(a);
	}
	
	/** Sets the medians for all cells. */
	private void setMedians() {
		for(Enumeration<Cell> e = root.depthFirstEnumeration(); e.hasMoreElements(); ) {
			Cell c = e.nextElement();
			Vector<CellData> v = c.getCellData();
			for(int i=0; i<v.size(); i++) {
				CellData cd = (CellData) v.elementAt(i);
				cd.iNucleus.rweight = (int) median.get(cd).intValue();
			}
		}
	}
}
