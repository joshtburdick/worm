package org.rhwlab.acetree;
import java.io.IOException;

/*
 * Created on Jan 14, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/**
 * @author biowolp
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AnnotInfo {

    /**
     * 
     */
    public String iName;
    public int iX;
    public int iY;
        
    public AnnotInfo(String name, int x, int y) {
            iName = name;
            iX = x;
            iY = y;
    }
    
    public String toString() {
        String s = iName + ", " + iX + ", " + iY;
        return s;
    }
    
    protected void finalize() throws IOException {
        //System.out.println("AnnotInfo: " + iName);
        //this = null;
    }

}
