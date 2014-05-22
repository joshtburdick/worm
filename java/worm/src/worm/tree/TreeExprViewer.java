package worm.tree;

import java.awt.*;
import javax.swing.*;
import java.io.*;
import java.util.*;

/** Views expression in a tree. */
public class TreeExprViewer extends JFrame {

	private JScrollPane sp;
	
	private TreePanel tp;
	
	private TreeExprViewerPrefs prefs;
	
	/** Constructor. */
	public TreeExprViewer() throws FileNotFoundException, IOException {
		super("Tree expression viewer");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		sp = new JScrollPane();
		// following advice of
		// http://stackoverflow.com/questions/12781179/moving-jscrollpane-horizontally-results-in-blured-text
		sp.getViewport().setScrollMode(JViewport.SIMPLE_SCROLL_MODE);
		add(sp);
		
		// image for testing
		// XXX FIXME this should totally have, like, a user interface or something
		String[] files = {
//				"/home/jburdick/data/image/CSV///SCD20091117_tag-185_b_6_L1.csv",
//				"/gpfs/fs0/u/twal/images//20120331_RW10871_L3/dats/SCD20120331_RW10871_L3.csv"
//	 "/gpfs/fs0/u/jlrichar/images/20110408_RW10425_pha-4_L2/dats/SCD20110408_RW10425_pha-4_L2.csv"  // pha-4
       "/gpfs/fs0/u/twal/images//20120331_RW10871_L3/dats/CD20120331_RW10871_L3.csv"
//	 "/gpfs/fs0/u/twal/images//20111011_L1/dats/SCD20111011_L1.csv",   // hlh-16
//				"/gpfs/fs0/u/twal/images//20120331_RW10871_L3/dats/SCD20120331_RW10871_L3.csv", // ceh-6
//	 "/gpfs/fs0/u/twal/images//20120426_JIM122_L2/dats/SCD20120426_JIM122_L2.csv"  // ceh-26
//		"/gpfs/fs0/u/azach/images//20120714_JIM136_L3/dats/SCD20120714_JIM136_L3.csv"  // JIM-136
		};
		
		tp = createPanelFromFilenames(files);
		// tp = new TreePanel();   // XXX
		prefs = new TreeExprViewerPrefs(tp);
		add(prefs, BorderLayout.NORTH);
		sp.getViewport().add(tp);
		
		/*
		CDFileParser p = new CDFileParser(new FileReader(s));
		HashMap<String, Cell> h = p.parseFile();
		
		TreePanelSignal si = new TreePanelSignal(h);

		tp.signal.add(si);
		GaussianSmoother gs = new GaussianSmoother(tp.signal.elementAt(0).cell.get("P0"));
		
		// XXX these objects shouldn't overlap
		tp.expr = si.cell;
		gs.smoothMaximizing(tp.expr, 0);
		*/
	}
	
	/** Creates a TreePanelSignal based on several movies, with default colors. */
	private TreePanel createPanelFromFilenames(String[] files)
	throws FileNotFoundException, IOException {
		Color[] colors = { Color.RED, Color.GREEN, Color.BLUE, Color.ORANGE, Color.CYAN };
		Vector<TreePanelSignal> s = new Vector<TreePanelSignal>();
		
		for(int i=0; i<files.length; i++) {
System.out.println("reading file " + i + ": " + files[i]);
			CDFileParser p = new CDFileParser(new FileReader(files[i]));
			TreePanelSignal si = new TreePanelSignal(p.parseFile());
			si.color = colors[i];
			System.out.println(si.color.toString());
			s.add(si);
		}
		
		return new TreePanel(s);
	}
	
	public Dimension getMinimumSize() {
		return new Dimension(900,500);
	}
	
	/** Main function. */
	public static void main(String[] args)
			throws FileNotFoundException, IOException {
		
		TreeExprViewer t = new TreeExprViewer();
		t.setSize(t.getMinimumSize());
//		t.pack();    // XXX this sets the image too large
		t.setVisible(true);
	}

}
