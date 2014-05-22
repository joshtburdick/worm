package worm.tree;

import java.awt.Color;
import java.util.HashMap;
import java.util.TreeMap;

public class TreePanelSignal {

	/** The color to use at max. intensity when drawing this. */
	Color color = Color.RED;
	
	/** Intensities which should be "black" and "color", respectively. */
	float loIntensity = 0f, hiIntensity = 2400f;
	
	/** The cells, indexed by name, then time. */
	HashMap<String, Cell> cell;

	/** Constructor. */
	TreePanelSignal(HashMap<String, Cell> cell) {
		this.cell = cell;
	}
}
