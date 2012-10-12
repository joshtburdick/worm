import java.util.Vector;
import java.util.Enumeration;
import java.util.Hashtable;
import org.rhwlab.snight.Nucleus;
import ij.gui.PlotWindow;

import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.analyze.PlotData;
import java.text.DecimalFormat;
import org.rhwlab.utils.C;
import java.io.File;
import gov.noaa.pmel.sgt.dm.SGTData;
import gov.noaa.pmel.sgt.dm.SGTLine;
import gov.noaa.pmel.sgt.dm.SGTMetaData;
import gov.noaa.pmel.sgt.dm.SimpleLine;
import gov.noaa.pmel.sgt.swing.JPlotLayout;
import gov.noaa.pmel.util.Domain;
import gov.noaa.pmel.util.Range2D;
import java.awt.BorderLayout;
import javax.swing.JFrame;
import java.awt.Dimension;



public class TestPlugin2 extends org.rhwlab.analyze.AnalysisTemplate {
	
    private double     iXA; // for angle()
    private double     iYA;
    private double     iZA;

    public void test1() {
		append("MinimumPlugin.test1 entered");
		append("iTextField1=" + iTextField1.getText());
		append("iTextField2=" + iTextField2.getText());
	}
	
	public void test2() {
		append("MinimumPlugin.test1 entered");
		append("iTextField3=" + iTextField3.getText());
		initialize();
		test3();
	}
	
    private void test3() {
	    plotAngleByIncrements();
    }
    
    private void plotAngleByIncrements() {
        String yLabel = "angle";
        String xLabel = "time";
        File f = new File(iNucleiMgr.getConfig().iConfigFileName);
        String title = f.getName();
        PlotData pd = angleByIncrements(iNucleiMgr);
        pd.showMe();
        JPlotLayout plotLayout = plotlayout(title, xLabel, yLabel, pd.xValues, pd.yValues);
        JFrame frame = new JFrame("test");
        frame.getContentPane().setLayout(new BorderLayout());
        frame.getContentPane().add(plotLayout, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);
    }
    
    private PlotData angleByIncrements(NucleiMgr nucMgr) {
        //append("angle");
        initialize();
        double f = 180./Math.PI;
        double zfac = nucMgr.getZPixRes();
        Vector nucleirecord = nucMgr.getNucleiRecord();
        double cumAng = 0;
        double avgAng = 0;
        int first = 1;
        int last = nuclei_record.size() - 1;
        Vector nuclei1 = (Vector)nucleirecord.elementAt(first);
        double [] xValues = new double[last - first];
        double [] yValues = new double[last - first];
        DecimalFormat fmt = new DecimalFormat("####.##");
        
        //for (int i=1; i < nuclei_record.length - 1; i++) {
        //System.out.println("angle: " + first + C.CS + last);
        for (int i=first; i < last; i++) {
            if (i > 50) break;
            getAvgs(nuclei1);
            double xa = iXA;
            double ya = iYA;
            double za = iZA;
            Vector nuclei2 = (Vector)nucleirecord.elementAt(i);
            getAvgs(nuclei2);
            //Vector pairs = new Vector();
            Nucleus n1 = null;
            Nucleus n2 = null;
            String s = null;
            String CS = C.CS;
            double r = 0;
            double raprod = 0;
            
            Enumeration e = nuclei1.elements();
            while (e.hasMoreElements()) {
                n1 = (Nucleus)e.nextElement();
                if (n1.status == Nucleus.NILLI) continue;
                //if (!(n1.identity.indexOf("E") == 0 || n1.identity.indexOf("MS") == 0)) continue;
                if (n1.successor1 <= 0 || n1.successor2 > 0) continue;
                //System.out.println("n1: " + n1);
                n2 = (Nucleus)nuclei2.elementAt(n1.successor1 - 1);
                if (!n1.identity.equals(n2.identity)) continue;
                
                //System.out.println("n2: " + n2);
                double dy1 = n1.y - ya;
                double dz1 = n1.z * nucMgr.getZPixRes() - za;
                double dy2 = n2.y - iYA;
                double dz2 = n2.z * nucMgr.getZPixRes() - iZA;
                double [] c1 = new double[2];
                double [] c2 = new double[2];
                c1 = complex(dy1, dz1);
                c2 = complex(dy2, dz2);
                double a1 = c1[1]; //angrad(dy1, dz1);
                double a2 = c2[1]; //angrad(dy2, dz2);
                double r1 = c1[0];
                
                double da = a2 - a1;
                if (da > Math.PI * 1.5) da -= 2 * Math.PI;
                else if (da < - 1.5 * Math.PI) da += 2 * Math.PI;
                if (i < 8) debugShow(n1.identity, ya, za, dy1, dz1, dy2, dz2, a1*f, a2*f, da*f);
                
                r1 = r1*r1;
                r += r1; //1;
                raprod +=  r1 * da;
                
                if (!n1.identity.equals(n2.identity)) {
                    System.out.println("bad data");
                }
            }
            //s = String.valueOf(raprod/r * f);
            //System.out.println();
            
            if (r > 0) avgAng = raprod / r * f;
            cumAng += avgAng;
            //System.out.println("time: " + i + CS + cumAng + CS + avgAng + CS + raprod + CS + r);
            append("time, cumAng, avgAng, " + i + CS + fmt.format(cumAng) + CS + fmt.format(avgAng)
                    + CS + fmt.format(raprod) + CS + fmt.format(r));
            nuclei1 = nuclei2;
            yValues[i - 1] = (double)cumAng;
            xValues[i - 1] = (double)i;
        }
        
        //return new PlotData(xValues, yValues);
        double [] xx = new double[50];
        System.arraycopy(xValues, 0, xx, 0, 50);
        double [] yy = new double[50];
        System.arraycopy(yValues, 0, yy, 0, 50);
        return new PlotData(xx, yy);
    }
    
    private double [] complex(double y, double z) {
        double [] ra = new double[2];
        ra[1] = angrad(y,z);
        ra[0] = Math.sqrt(y * y + z * z);
        return ra;
    }
    
    private double angrad(double y, double z) {
        double angr = Math.atan(y/z);
        if (y < 0.) {
            if (z < 0.) angr += Math.PI;
            else angr += 2. * Math.PI;
        } else {
            if (z < 0.) angr += Math.PI; 
        }
        return angr;
    }
    
    // computes the center of gravity of the nuclei at a time point
    // the center is stored in iXA, iYA, iZA
    private void getAvgs(Vector nuclei) {
        Enumeration e = nuclei.elements();
        Nucleus n = null;
        int count = 0;
        iXA = 0;
        iYA = 0;
        iZA = 0;
        while (e.hasMoreElements()) {
            n = (Nucleus)e.nextElement();
            if (n.status == Nucleus.NILLI) continue;
            count++;
            iXA += n.x;
            iYA += n.y;
            iZA += n.z * iNucleiMgr.getZPixRes();
        }
        if (count > 0) {
            iXA /= count;
            iYA /= count;
            iZA /= count;
        }
    }
    
    private void debugShow(String cell, double ya, double za, double dy1, double dz1, 
            double dy2, double dz2, double a1, double a2, double da) {
        DecimalFormat fmt = new DecimalFormat("####.##");
        String s = cell;
        s += C.CS + fmt.format(ya);
        s += C.CS + fmt.format(za);
        s += C.CS + fmt.format(dy1);
        s += C.CS + fmt.format(dz1);
        s += C.CS + fmt.format(dy2);
        s += C.CS + fmt.format(dz2);
        s += C.CS + fmt.format(a1);
        s += C.CS + fmt.format(a2);
        s += C.CS + fmt.format(da);
        append(s);
    }
    
    private JPlotLayout plotlayout(String title, String xLabel, String yLabel, double [] xValues, double [] yValues) {
        JPlotLayout layout_ = new JPlotLayout(false, false, false, "bogus", null, false);
        layout_.setSize(new Dimension(240,160));
        layout_.setMouseEventsEnabled(false);
        layout_.setBatch(true);
        layout_.setTitles(title, "", "");
        layout_.setTitleHeightP(0.4, 0.4);
        SimpleLine data = new SimpleLine((double [])xValues, yValues, "CumAngle ");
        SGTMetaData meta = new SGTMetaData(xLabel, "", false, false);
        data.setXMetaData(meta);
        meta = new SGTMetaData(yLabel, "", false, false);
        data.setYMetaData(meta);
        layout_.addData(data, "cum angle");
        Domain d = layout_.getRange();
        //System.out.println(title + C.CS  + d.getYRange());
        /*
        d.setYRange(new Range2D(-80, 10));
        try {
            layout_.setRange(d);
        } catch(Exception e) {
            e.printStackTrace();
        }
        */
        
        layout_.setBatch(false);
        return layout_;
        
    }
    
    

}
