package worm.tree;

import java.awt.*;
import java.awt.color.*;
import java.awt.geom.*;
import java.awt.print.*;
import javax.swing.JPanel;
import javax.vecmath.*;

import java.util.*;

/** Draws a lineage tree colored by expression.
 * Inspired by EmbryoDB tree-drawing code. */
public class TreePanel extends JPanel implements Printable {

	/** Information about the signals that will be shown in this
	 * panel. This includes colors and intensity thresholds.
	 * It also includes the expression data, but this component
	 * will ignore that, and just render the data in expr. */
	public Vector<TreePanelSignal> signal =
			new Vector<TreePanelSignal>();
	
	/** The expression of all signals that will be drawn.
	 * This should include any smoothing. */
	public HashMap<String, Cell> expr;
	
	/** The name of the root to redraw from. */
	public String root = "P0";
	
	/** The selected cell (if any.) */
	private Cell selectedCell;
	
	/** The last timepoint to include. */
	public float lastT = 100;
	
	/** The selected timepoint (if any.) */
	// currently not used
	// private double selectedTime = -1;
	
	/** Used for drawing. */
	private Line2D.Float line = new Line2D.Float();
	
	private float[] col;
	
	/** The scale at which to draw this. */
	public float scale = 1;
	
	/** The spacing of leaf cells. */
	public float lineSpacing = 13;
	
	/** Width of lines. */
	public float lineWidth = 5;
	
	/** Should names be drawn? */
	public boolean showNames = true;
	
	/** Should we show the maximum intensity seen so far? */
	public boolean showMaxIntensitySoFar = false;

	/** Mapping from cell to x coordinate. */
	private HashMap<Cell, double[]> cellX;
	
	/** These are slightly different, to avoid a tiny
	 * overly-long line. */
	private Stroke horizontalStroke, verticalStroke;
	
	/** The transform from (leaf cell, time) to user-space x,y.
	 * Note that we transform these "by hand", rather than using
	 * Graphics2D.setTransform(), so that things like line widths
	 * and text aren't scaled. */
	AffineTransform at = new AffineTransform();
	
	Map<String, String> anatomyInfo = new HashMap<String, String>();
	
    /** The scale at which to print this. */
    double printingScale = 0.1f;
    
	/** Constructor. */
	public TreePanel() {
		col = new float[3];
		col[1] = 1;
		col[2] = 0;

		this.setBackground(new Color(250,250,250));    // now that we have three colors...
	}
	
	/** Constructor, in which this is initialized with some signals.
	 * @param signal  the signals to add */
	public TreePanel(Vector<TreePanelSignal> signal) {
		this.signal = signal;
		
		// FIXME this has incorrect sharing
		expr = signal.elementAt(0).cell;
		
		// add each signal
		for(int i=0; i<signal.size(); i++) {
			GaussianSmoother gs =
					new GaussianSmoother(signal.elementAt(i).cell.get("P0"));
			gs.smoothMaximizing(expr, i);
		}
	}
	
	public Dimension getPreferredSize() {
		// ??? base this on the actual dimensions?
		return new Dimension(9000,1000);
	}
	
	/** Draws this. */
	public void paintComponent(Graphics g1) {
		if (expr == null)
			return;
		
		super.paintComponent(g1);
		Graphics2D g = (Graphics2D) g1;
		
		// "allocate" x-coordinates
		XCoordComputer xc = new XCoordComputer(expr.get(root), lastT);
		cellX = xc.getCellX();
		
		// set up the transform, and some line widths
		at.setToScale(lineSpacing, 1);
		at.translate(1, 1);
		
		// these were different in an attempt to work around a fairly
		// miniscule problem in rendering; for now, ignoring this
		horizontalStroke = new BasicStroke(lineWidth,
			BasicStroke.CAP_SQUARE,
			BasicStroke.JOIN_BEVEL);
		verticalStroke = new BasicStroke(lineWidth,
				BasicStroke.CAP_SQUARE,
				BasicStroke.JOIN_BEVEL);
		
		redrawSignal(g, expr.get(root), new Color3f(Color.BLACK));
	}
	
	/** Defines the color in some cell at a particular time.
	 * XXX this may well need some optimizing. */
	public Color3f getColorAtCellTime(String c, double t) {
		
		// get the expression data for all signals
		Nucleus n = expr.get(c).nucleus.get(t);
		
		// initialize color to black
		Color3f color = new Color3f();
		
		// loop through the signals
		for(int i=0; i<signal.size(); i++) {
			TreePanelSignal s = signal.elementAt(i);
			
			// scale signal
			float signal = (((float) n.expr[i]) - s.loIntensity) /
				(s.hiIntensity - s.loIntensity);
			if (signal < 0) signal = 0;
			if (signal > 1) signal = 1;
		
			Color3f c1 = new Color3f(s.color);
			c1.scale(signal);
			color.add(c1);
		}
		
		// ensure that this is within bounds
		color.clamp(0, 1);
		
		return color;
	}
	
	/** Redraws this, rooted from a given cell.
	 * @param g  Graphics object to draw onto */
	public void redrawSignal(Graphics2D g, Cell c, Color3f maxColorSoFar) {
		Color3f maxColor1 = new Color3f(maxColorSoFar);

		// if this is after cutoff time, then return
		if (c.nucleus.firstKey().doubleValue() > lastT)
			return;
		
		// if this is a leaf cell, or truncated, stop and draw label
		if (showNames && (c.children.isEmpty() ||
				c.nucleus.lastKey().doubleValue() >= lastT)) {

			// compute y location
			float labelT = c.nucleus.lastKey().floatValue();
			if (labelT > lastT)
				labelT = lastT;
			
			Point2D p = new Point2D.Float((float) cellX.get(c)[1], labelT);
			at.transform(p, p);
			
			AffineTransform t = g.getTransform();

			g.translate(p.getX() - lineWidth, p.getY() + lineWidth);
			g.rotate(Math.PI / 2);
			g.setColor(Color.BLACK);
			g.drawString(c.name, 0, 0);  // s.cellX.get(c).floatValue(), (float) c.t);
			
			g.setTransform(t);
		}

		// FIXME: if this is the selected cell, draw selection
		if (c == selectedCell) {
		
		}
		
		double previousT = -1 ;
		Color3f previousColor = new Color3f();   // to avoid writing each segment
		
		Point2D	line1 = new Point2D.Float(), line2 = new Point2D.Float();
		
		g.setStroke(verticalStroke);
		
		Color3f col = new Color3f();
		
		for (Map.Entry<Double, Nucleus> e: c.nucleus.entrySet()) {
			if (previousT > 0 && previousT <= lastT) {
				
				// set color
				col.set(getColorAtCellTime(c.name, e.getKey()));
				setToElementwiseMax(maxColor1, col);

				if (showMaxIntensitySoFar)
					col.set(maxColor1);
		
				// ??? only actually draw this if the color has changed, or if
				// this is the last timepoint for this cell ?
				if (true) {
				
					// draw vertical line
					// XXX the ends are, like, one pixel too long in PDF output
					// for now, not worrying about this
					g.setColor(col.get());
					line1.setLocation(cellX.get(c)[1], previousT);
					line2.setLocation(cellX.get(c)[1], e.getValue().t);
					at.transform(line1, line1);
					at.transform(line2, line2);
					line.setLine(line1, line2);		
					g.draw(line);
				}
			}
			
			previousT = e.getKey();
			previousColor.set(col);
		}
		
		// draw everything below here
		for (Cell child: c.children) {
		
			// possibly draw horizontal lines
			if (c.nucleus.lastKey() < lastT) {
				g.setStroke(horizontalStroke);
				g.setColor(col.get());
				line1.setLocation(cellX.get(c)[1], c.nucleus.lastKey());
				line2.setLocation(cellX.get(child)[1], c.nucleus.lastKey());
				at.transform(line1, line1);
				at.transform(line2, line2);
				line.setLine(line1, line2);
				g.draw(line);
			}
			
			// recursive call
			redrawSignal(g, child, maxColor1);
		}		
	}

	/** Utility to compute maximum.
	 * @param a, b  Color3f objects
	 *
	 * Side effects: sets a to element-wise max of a and b.
	 */
	private void setToElementwiseMax(Color3f a, Color3f b) {
		if (b.x > a.x) a.x = b.x;
		if (b.y > a.y) a.y = b.y;
		if (b.z > a.z) a.z = b.z;
	}
	
	/** Allows printing. */
	public int print(Graphics g, PageFormat pf, int page)
	    throws PrinterException {
	    if (page > 0) {
	        return NO_SUCH_PAGE;
	    }
//	    pf.setOrientation(PageFormat.PORTRAIT);

	    Graphics2D g2d = (Graphics2D)g;
/*
	    g2d.translate(pf.getImageableX(), pf.getImageableY());
	    g2d.scale(0.05f,0.05f);
*/
	    // attempting to scale this to fit the printable area
	    g2d.translate(pf.getImageableX(), pf.getImageableY());
//	    g2d.rotate(Math.PI/2, pf.getImageableX() / 2, pf.getImageableY() /2 );
	    g2d.scale(printingScale, printingScale);
	    printAll(g);
	
	    return PAGE_EXISTS;
	}
}
