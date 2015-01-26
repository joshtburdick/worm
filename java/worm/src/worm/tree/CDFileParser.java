package worm.tree;

import java.io.*;
import java.util.*;

/** Parses a CD file (currently an SCD file.) */
public class CDFileParser {

	private BufferedReader r;
	
	/** Cells indexed first by name, then by time. */
	private HashMap<String, Cell> cells =
			new HashMap<String, Cell>();
	
	/** Line number. */
	int lineNumber = 0;
	
	/** Constructor. */
	public CDFileParser(Reader r) {
		this.r = new BufferedReader(r);
	}
	
	/** Parses a file. */
	public HashMap<String, Cell> parseFile()
			throws IOException {

		if (!checkHeader())
			return null;
		
		// read in all the lines
		String line = r.readLine();
		while (line != null) {
			parseLine(line);
			line = r.readLine();
			lineNumber++;
		}
		
		// add various links
		linkCells();
		
		return cells;
	}
	
	/** Checks the header. */
	public boolean checkHeader() throws IOException {
		String line = r.readLine();
		return (line.startsWith("cellTime,cell,time,none,global,local,blot,cross,z,x,y,size,gweight"));
	}
	
	/** Parses one line, and adds a Cell. */
	private void parseLine(String line) {
		String[] s = line.split(",");   // ??? use a Pattern here?
		assert( s.length >= 12);
		
		String cell = s[1];
		double time = Double.parseDouble(s[2]);
		
		// create Cell, if it's not already present
		if (!cells.containsKey(cell)) {
//			System.out.print(cell + " ");
			cells.put(cell, new Cell(cell));
		}
		
		// new Nucleus containing all this data
		Nucleus n = new Nucleus();
		n.index = lineNumber;
		n.t = time;
		n.x = Double.parseDouble(s[9]);
		n.y = Double.parseDouble(s[10]);
		n.z = Double.parseDouble(s[8]);
		n.size = Double.parseDouble(s[11]);
		n.expr = new double[5];
		for(int i=0; i<5; i++)
			n.expr[i] = Double.parseDouble( s[ i+3 ]);

		cells.get(cell).nucleus.put(time, n);
//		System.out.println(cell + " " + h.get(cell).nucleus.size());
	}
	
	/** Links up cells with their descendants. */
	private void linkCells() throws IOException {
		
		// go through the list of canonical lineage
		// relationships, and add links in it
		linkTriples();
		fakeTreeTop();
	}
	
	/** Links up lineage triples. */
	private void linkTriples() throws IOException {
		InputStream is = this.getClass().getResourceAsStream("lineageTriples.txt");
		BufferedReader r = new BufferedReader(new InputStreamReader(is));
		
		String line = r.readLine();
//		System.out.println(line);
		
		line = r.readLine();
		while (line != null) {
			String s[] = line.split("\t");
			String parent = s[0];
			String left = s[1];
			String right = s[2];
//			System.out.println(parent + " " + left + " " + right);
			// technically, I think these should always be present...
//			System.out.println(cells.get(parent));
//			System.out.println(cells.get(parent).children);
			
			if (cells.containsKey(s[0]) && cells.containsKey(s[1]))
				cells.get(s[0]).children.add(cells.get(s[1]));
			if (cells.containsKey(s[0]) && cells.containsKey(s[2]))
				cells.get(s[0]).children.add(cells.get(s[2]));
			
			line = r.readLine();
		}
	}
	
	/** Ensures that the top of the tree is present. */
	private void fakeTreeTop() {
		ensurePresent("ABa", "ABal", "ABar");
		ensurePresent("ABp", "ABpl", "ABpr");
		ensurePresent("EMS", "MS", "E");
		ensurePresent("P2", "C", "P3");
		ensurePresent("AB", "ABa", "ABp");
		ensurePresent("P1", "EMS", "P2");
		ensurePresent("P0", "AB", "P1");
	}
	
	/** Checks that a link is present, and "fakes" a cell at time 0 if it isn't. */
	private void ensurePresent(String parent, String left, String right) {

		// if the parent is present, we're done
		if (cells.containsKey(parent)) {
			System.out.println( parent + " is already present" );
			return;
		}
		
		// otherwise, create a "fake" cell
		Cell c = new Cell(parent);
		// ... with a nucleus, with all numbers set to 0
		Nucleus n = new Nucleus();
		n.expr = new double[5];
		c.nucleus.put(0d, n);

		c.children.add(cells.get(left));
		c.children.add(cells.get(right));
		cells.put(parent, c);
	}
}
