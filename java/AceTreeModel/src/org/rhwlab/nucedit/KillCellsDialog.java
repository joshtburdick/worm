/*
 * Created on Apr 8, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.rhwlab.nucedit;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.Border;

import org.rhwlab.acetree.AceTree;
import org.rhwlab.acetree.NucUtils;
import org.rhwlab.snight.Identity;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.AncesTree;
import org.rhwlab.tree.Cell;
import org.rhwlab.utils.Log;

/**
 * @author biowolp
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class KillCellsDialog extends GenericDialog {
    JTextField 			iCellToKill;
    JTextField 			iKillTime;
    JRadioButton 		iOneCell;
    JRadioButton 		iSeveralCells;
    JLabel 				iCellCount;
    JSpinner       		iCellsToKill;
    Vector 				iCandidateCells;
    int 				iTime;
    String 				iCellName;
    EditLog 			iEditLog;
    int					iNumCellsToEnd;
    JButton				iUseAll;
    SpinnerNumberModel	iSNModel;
    //NucleiMgr iNucleiMgr;


    /**
     * @param aceTree
     * @param owner
     * @param modal
     */
    public KillCellsDialog(AceTree aceTree, Frame owner, boolean modal,
            Cell cell, int time, EditLog editLog) {
        super(aceTree, owner, modal);
        Border blackline = BorderFactory.createLineBorder(Color.black);
        iEditLog = iNucleiMgr.getEditLog();
        iCandidateCells = new Vector();
        iTime = time;
        iCellName = cell.getName();
        setTitle(TITLE);
        iContentPanel.setLayout(new BoxLayout(iContentPanel, BoxLayout.PAGE_AXIS));
        JPanel p = new JPanel();
        p.setLayout(new GridLayout(1,0));
        p.setBorder(blackline);
        JLabel test = new JLabel(CELLTOKILL);
        p.add(test);
        iCellToKill = new JTextField();
        iCellToKill.setColumns(12);
        //iCellToKill.setMaximumSize(new Dimension(20, 200));
        iCellToKill.setText(cell.getName());
        p.add(iCellToKill);
        iContentPanel.add(p);
        p = new JPanel();
        p.setLayout(new GridLayout(1,0));
        p.setBorder(blackline);
        test = new JLabel(KILLTIME);
        //iContentPanel.add(test);
        iKillTime = new JTextField();
        iKillTime.setColumns(12);
        iKillTime.setText(String.valueOf(time));
        //iContentPanel.add(iKillTime);
        p.add(test);
        p.add(iKillTime);
        iContentPanel.add(p);

        JPanel choices = new JPanel();
        choices.setLayout(new GridLayout(1,0));
        choices.setBorder(blackline);

        /*
        choices.setLayout(new GridLayout(0,1));
        iOneCell = new JRadioButton(ONECELL);
        iOneCell.setSelected(true);
        iOneCell.addActionListener(this);
        iSeveralCells = new JRadioButton(SEVERAL);
        iSeveralCells.addActionListener(this);
        ButtonGroup group = new ButtonGroup();
        group.add(iOneCell);
        group.add(iSeveralCells);
        choices.add(iOneCell);
        choices.add(iSeveralCells);
        */

        // put starting cell on the list and then look for more
        // up to the point where the name changes
        //System.out.println("KillCellsDialog: " + iNucleiMgr + CS + iTime);
        Vector nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(iTime - 1);
        iCandidateCells = new Vector();
        iCandidateCells.add(nuclei);
        collectCandidateCells();
        int nmax = iCandidateCells.size();
        iNumCellsToEnd = nmax;
        //System.out.println("at time: " + iTime);
        //System.out.println("killable: " + nmax);
        iCellCount = new JLabel(NUMBER);
        choices.add(iCellCount);
        iSNModel = new SpinnerNumberModel(1, 0, nmax, 1);
        iCellsToKill = new JSpinner(iSNModel);
        choices.add(iCellsToKill);

        iUseAll = new JButton("useAll");
        iUseAll.addActionListener(this);
        choices.add(iUseAll);

        iContentPanel.add(choices);


        setSize(new Dimension(WIDTH, HEIGHT));
        setVisible(true);

    }

    private void collectCandidateCells() {
        int k = iTime; // the first time is already in there
        int count = 1;
        boolean found = true;
        String name = iCellName;
        Nucleus n = null;
        Vector nuclei = null;
        Vector nuclei_record = iNucleiMgr.getNucleiRecord();
        while (found && (count < MAX_KILLS) && k < nuclei_record.size()) {
            nuclei = (Vector)nuclei_record.elementAt(k++);
            found = false;
            for (int j=0; j < nuclei.size(); j++) {
                n = (Nucleus)nuclei.elementAt(j);
                if (!n.identity.equals(name)) continue;
                if (n.status <= 0) break;
                iCandidateCells.add(nuclei);
                //n.status = Nucleus.NILLI;
                found = true;
            }
            count++;
        }

    }

    public void actionPerformed(ActionEvent e) {
        Object o = e.getSource();
        if (o == iUseAll) {
        	Object oo = iSNModel.getValue();
        	println("actionPerformed, " + oo);
        	iSNModel.setValue(new Integer(iNumCellsToEnd));
        	println("actionPerformed, " + oo);
        	return;
        }

        Integer xx = (Integer)iCellsToKill.getValue();
        int k = xx.intValue();
        if (k == 0) return;
        //System.out.println("actionPerformed: " + x);
        // locate the cell to start at after a rebuild
        String indexString = iNucleiMgr.getIndex(iCellName, iTime);
        //String indexString = " (" + index + ") ";
        iEditLog.appendx("KILLING " + iCellName + indexString + " at time " + iTime);
        //int k = iCandidateCells.size();
        if (k == 1) iEditLog.append(".");
        else {
            iEditLog.appendx(" and following ");
            if (k == 2) iEditLog.append("time." + iNucleiMgr.getIndex(iCellName, iTime + 1));
            else {
                int km = k - 1;
                String s = km + " times ";
                for (int i = 0; i < km; i++) {
                    s += iNucleiMgr.getIndex(iCellName, iTime + i + 1);

                }
                iEditLog.append(s);
            }
        }
        iNucleiMgr.makeBackupNucleiRecord();
        Nucleus predecessorNuc = null;
        if (iTime > 1) {
            Nucleus nc = iNucleiMgr.getCurrentCellData(iCellName, iTime);
            Vector predNuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(iTime - 2);
            if (nc.predecessor > 0) {
                predecessorNuc = (Nucleus)predNuclei.elementAt(nc.predecessor - 1);
            }
        }
        int namingMethod = AceTree.getAceTree(null).getNucleiMgr().getIdentity().getNamingMethod();
        for (int i=0; i < k; i++) {
            Vector nuclei = (Vector)iCandidateCells.elementAt(i);
            Nucleus n = null;
            //Vector indices = new Vector();
            for (int j=0; j < nuclei.size(); j++) {
                n = (Nucleus)nuclei.elementAt(j);
                if (!n.identity.equals(iCellName)) continue;
                n.status = Nucleus.NILLI;
                n.identity = "";
                n.assignedID = "";
                //indices.add(new Integer(j));
            }
            /*
            if (namingMethod == Identity.MANUAL) {
                int lastp = indices.size();
                for (int j = lastp; j > 0; j--) {
                    int kn = ((Integer)indices.elementAt(j - 1)).intValue();
                    nuclei.remove(kn);
                }
            }
            */
        }
        Cell c = null;
        int strTime = iTime - 1;
        if (predecessorNuc != null) {
            AncesTree ances = iAceTree.getAncesTree();
            Hashtable h = ances.getCellsByName();
            c = (Cell)h.get(predecessorNuc.identity);
        }

        if (o == iApplyAndRebuild) {
            iAceTree.clearTree();
            iAceTree.buildTree(true);
            iEditLog.setModified(true);
            AncesTree ances = iAceTree.getAncesTree();
            Hashtable h = ances.getCellsByName();
            if (c != null) c = (Cell)h.get(c.getName());
            System.out.println("killCellsDialog.actionPerformed: " + c + CS + strTime);
            if (c != null) iAceTree.setStartingCell(c, strTime);
            dispose();
        } else super.actionPerformed(e);

    }

    private final static int
    WIDTH = 300
   ,HEIGHT = 300
   ,MAX_KILLS = 999
   ;

    public final static String
         TITLE = "Kill Cells"
        ,ONECELL = "one cell"
        ,SEVERAL = "several cells"
        ,KILLABLE = " are killable"
        ,CELLTOKILL = "kill cell"
        ,KILLTIME = "start time"
        ,NUMBER = "number"
        ,SPACER = "     "
        ;

    public static void main(String[] args) {
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
