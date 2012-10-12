import java.util.Vector;
import java.util.Enumeration;
import java.util.Hashtable;
import org.rhwlab.snight.Nucleus;
import ij.gui.PlotWindow;


public class TestPlugin extends org.rhwlab.analyze.AnalysisTemplate {
	
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
        append("test3 entered");
        float [] xValues = new float[nuclei_record.size()];
        float [] yValues = new float[nuclei_record.size()];
        String yLabel = "cell count";
        String xLabel = "time";
        String title = "cell count vs time";
        for (int i=0; i < nuclei_record.size(); i++) {
            Vector nuclei = (Vector)nuclei_record.elementAt(i);
            yValues[i] = countNuclei(nuclei);
            xValues[i] = i;
        }
        
        
        PlotWindow pw = new PlotWindow(title, xLabel, yLabel, xValues, yValues);
        pw.setLimits(0, 250, 0, 500);
        pw.draw();
    }
    
    private int countNuclei(Vector nuclei) {
        int count = 0;
        Enumeration e = nuclei.elements();
        while(e.hasMoreElements()) {
            Nucleus n = (Nucleus)e.nextElement();
            if (n.status > 0) count++;
        }
        return count;
    }
    

}
