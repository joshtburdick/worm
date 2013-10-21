package edu.murrlab.image.translate;

import java.io.*;

import loci.common.DataTools;

/** Extracts the raw XML from a Leica .lif file.
 * Heavily based on the loci.formats.in.LIFReader class.
 * @author jburdick
 */
public class LeicaLIFXMLReader {

	/** Stream we're reading from. */
	DataInputStream is;
	
	private static final byte LIF_MAGIC_BYTE = 0x70;
	
	private static final byte LIF_MEMORY_BYTE = 0x2a;
	
	public LeicaLIFXMLReader(File file) throws IOException {
		this.is = new DataInputStream(new FileInputStream(file));
	}
	
	/** Parses the XML. */
	public String getXML() throws IOException {

		// read the header
		byte checkOne = is.readByte();
		is.skipBytes(2);
		byte checkTwo = is.readByte();
		if (checkOne != LIF_MAGIC_BYTE && checkTwo != LIF_MAGIC_BYTE) {
			throw new IllegalArgumentException("not a valid Leica LIF file");
		}
		is.skipBytes(4);

		// read and parse the XML description
		if (is.readByte() != LIF_MEMORY_BYTE) {
			throw new IllegalArgumentException("Invalid XML description");
		}

		// number of Unicode characters in the XML block
		int nc = Integer.reverseBytes(is.readInt());
//		System.err.println("nc = " + nc);
		
		// read in that many bytes
		byte[] buf = new byte[ 2 * nc ];
		is.readFully(buf);

		// as far as I know, this just strips out every other byte,
		// to convert from Unicode
		String xml = DataTools.stripString(new String(buf));

		return xml;
	}
}
