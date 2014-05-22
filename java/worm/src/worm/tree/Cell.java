package worm.tree;

import java.util.*;

/** Represents a cell (and indeed, a lineage tree.) */
public class Cell {

	/** The name of this cell. */
	public String name;
	
	/** The nuclei in this cell, indexed by time. */
	public TreeMap<Double, Nucleus> nucleus =
			new TreeMap<Double, Nucleus>();
	
	/** The daughters of this cell. */
	public Vector<Cell> children = new Vector<Cell>();
	
	public Cell(String name) {
		this.name = name;
	}
	
}
