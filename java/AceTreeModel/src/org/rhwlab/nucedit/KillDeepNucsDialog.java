package org.rhwlab.nucedit;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DecimalFormat;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.Border;

import org.rhwlab.acetree.AceTree;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.Nucleus;

public class KillDeepNucsDialog  extends JDialog implements ActionListener {

	AceTree			iAceTree;
	NucleiMgr		iNucleiMgr;
	int				iZLim;
	int				iCount;
	JLabel			iEstimatedCount;
	JLabel			iZLimLabel;
	JCheckBox		iFlippedImages;


	public KillDeepNucsDialog(AceTree aceTree, Frame owner, boolean modal) {
        super(owner, modal);
        iAceTree = aceTree;
        iNucleiMgr = iAceTree.getNucleiMgr();
        setTitle("Kill Deep Nucs");

        JDialog dialog = this;
        iZLim = 27;
        //estimateNucs(false);
        //iEstimatedCount = new JLabel(String.valueOf(iCount));

        JPanel pWhole = new JPanel();

        pWhole.setLayout(new BoxLayout(pWhole, BoxLayout.PAGE_AXIS));
        Border blackline = BorderFactory.createLineBorder(Color.black);
        //Border empty = BorderFactory.createEmptyBorder();
        Border topBorder = BorderFactory.createEmptyBorder(10,0,0,0);
        Border botBorder = BorderFactory.createEmptyBorder(0,0,10,0);
        JPanel p = new JPanel();
        iFlippedImages = new JCheckBox("Flipped images?");
        p.setLayout(new GridLayout(1,0));
        p.add(iFlippedImages);
        pWhole.add(p);
        estimateNucs(false);
        iEstimatedCount = new JLabel(String.valueOf(iCount));

        p = new JPanel();
        p.setLayout(new GridLayout(1,0));
        p.add(new JLabel("z limit:"));
        iZLimLabel = new JLabel(String.valueOf(iZLim));
        p.add(iZLimLabel);
        JButton plus = new JButton("+");
        JButton minus = new JButton("-");
        plus.addActionListener(this);
        minus.addActionListener(this);
        p.add(plus);
        p.add(minus);
        pWhole.add(p);

        p = new JPanel();
        p.setLayout(new GridLayout(1,0));
        JButton estimate = new JButton("estimate");
        p.add(estimate);
        estimate.addActionListener(this);
        pWhole.add(p);


        p = new JPanel();
        p.setLayout(new GridLayout(1,0));
        p.add(new JLabel("Estimated count: "));
        p.add(iEstimatedCount);
        pWhole.add(p);

        p = new JPanel();
        p.setLayout(new GridLayout(1,0));
        JButton killem = new JButton("Kill'em");
        killem.addActionListener(this);
        p.add(killem);
        pWhole.add(p);

        pWhole.setOpaque(true); //content panes must be opaque
        dialog.setContentPane(pWhole);
        //dialog.setSize(new Dimension(200, 400));
        dialog.setLocationRelativeTo(owner);
        dialog.pack();
        dialog.setVisible(true);


	}

	private void estimateNucs(boolean implement) {
		println("estimateNucs, " + iZLim + CS + iCount + CS + implement);
		iCount = 0;
        Vector nucRec = (Vector)iNucleiMgr.getNucleiRecord();
        for (int i=0; i < nucRec.size(); i++) {
        	Vector nuclei = (Vector)nucRec.get(i);
        	for (int j=0; j < nuclei.size(); j++) {
        		Nucleus n = (Nucleus)nuclei.get(j);
        		if (n.status == Nucleus.NILLI) continue;
        		if (n.z < iZLim && !iFlippedImages.isSelected()) continue;
        		if (n.z > iZLim && iFlippedImages.isSelected()) continue;
        		//println("killDeepNucs, " + i + CS + n);
        		if (implement) {
        			n.status = Nucleus.NILLI;
        		} else iCount++;
        	}
        }
        if (implement) {
            iAceTree.clearTree();
            iAceTree.buildTree(true);

        }
		//println("estimateNucs, " + iZLim + CS + iCount);
	}

	public void actionPerformed(ActionEvent e) {
		String c = e.getActionCommand();
		if (c.equals("+")) iZLim++;
		else if (c.equals("-")) iZLim--;
		else if (c.equals("estimate")) {
			estimateNucs(false);
			iEstimatedCount.setText(String.valueOf(iCount));
		} else if (c.equals("Kill'em")) {
			estimateNucs(true);
			iEstimatedCount.setText(String.valueOf(iCount));
		}

		iZLimLabel.setText(String.valueOf(iZLim));

	}



	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

    private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ";
    private static final String TAB = "\t";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt4(double d) {return DF4.format(d);}
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF0.format(d);}

}
