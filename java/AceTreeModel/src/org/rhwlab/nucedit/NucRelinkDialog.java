/*
 * Created on Mar 29, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.rhwlab.nucedit;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.KeyboardFocusManager;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.border.Border;

import org.rhwlab.acetree.AceTree;
import org.rhwlab.acetree.NucUtils;
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
public class NucRelinkDialog extends JDialog implements ActionListener {
    private AceTree 			iAceTree;
    public static NucleiMgr 	iNucleiMgr;
    private JLabel			iRelinkTime;
    private JLabel 			iRelinkNuc;
    private JLabel 			iLinkTime;
    private JLabel 			iLinkNuc;
    //private JButton 			iDoit;
    private JButton 			iApplyAndRebuild;
    private JButton 			iApplyOnly;
    private JButton 			iRelinkButton;
    private JButton 			iLinkButton;
    private JButton 			iLinkRootButton;
    private EditLog 			iEditLog;
    //private Log     iDLog;


    public NucRelinkDialog(AceTree aceTree, Frame owner, boolean modal,
        Cell cell, int time) {
        super(owner, modal);
        iAceTree = aceTree;
        iNucleiMgr = iAceTree.getNucleiMgr();
        iEditLog = iAceTree.getEditLog();
        //iDLog = iAceTree.getDebugLog();
        setTitle(TITLE);

        JDialog dialog = this;
        //dialog.setFocusable(false);
        addWindowListener(new WindowEventHandler());

        JPanel pWhole = new JPanel();


        pWhole.setLayout(new BoxLayout(pWhole, BoxLayout.PAGE_AXIS));
        Border blackline = BorderFactory.createLineBorder(Color.black);
        //Border empty = BorderFactory.createEmptyBorder();
        Border topBorder = BorderFactory.createEmptyBorder(10,0,0,0);
        Border botBorder = BorderFactory.createEmptyBorder(0,0,10,0);

        JPanel labelAtTop = new JPanel();
        labelAtTop.setLayout(new GridLayout(2,1)); //labelAtTop, BoxLayout.PAGE_AXIS));
        labelAtTop.setBorder(blackline);
        JLabel topLab = new JLabel(LATER);
        Font f = topLab.getFont();
        //iDLog.append("NucRelinkDialog: " + f.getName() + CS + f.getStyle() + CS + f.getSize());
        int size = (int)(f.getSize() * 1.3);
        //iDLog.append("NucRelinkDialog: " + size);

        Font f2 = new Font(f.getName(), f.getStyle(), size);
        topLab.setFont(f2);
        labelAtTop.add(topLab);
        JLabel botLab = new JLabel(EARLIER);
        botLab.setFont(f2);
        labelAtTop.add(botLab);
        pWhole.add(labelAtTop);


        //earlier time
        JPanel pStr = new JPanel();
        pStr.setLayout(new BoxLayout(pStr, BoxLayout.PAGE_AXIS));
        pStr.setBorder(blackline);
        JPanel s = new JPanel();
        s.setLayout(new FlowLayout());
        JLabel label = new JLabel(LINKTIME);
        s.add(label);

        iLinkTime = new JLabel(FIVE);
        iLinkTime.setBorder(blackline);
        //iLinkTime.setColumns(5);
        s.add(iLinkTime);
        s.setBorder(topBorder);
        pStr.add(s);

        s = new JPanel();
        s.setLayout(new FlowLayout());
        label = new JLabel(LINKNUC);
        s.add(label);
        iLinkNuc = new JLabel(TWELVE);
        iLinkNuc.setBorder(blackline);
        //iLinkNuc.setColumns(12);
        s.add(iLinkNuc);
        s.setBorder(botBorder);
        pStr.add(s);

        iLinkButton = new JButton(SETEARLYCELL);
        iLinkButton.addActionListener(this);
        iLinkRootButton = new JButton(SETROOTCELL);
        iLinkRootButton.addActionListener(this);
        s = new JPanel();
        s.setLayout(new GridLayout(0,1));
        s.setBorder(topBorder);
        s.add(iLinkButton);
        s.add(iLinkRootButton);
        pStr.add(s);
        pWhole.add(pStr);

        // later time
        JPanel pEnd = new JPanel();
        pEnd.setLayout(new BoxLayout(pEnd, BoxLayout.PAGE_AXIS));
        pEnd.setBorder(blackline);
        s = new JPanel();
        s.setLayout(new FlowLayout());
        label = new JLabel(RELINKTIME);
        s.add(label);

        //iRelinkTime = new JTextField();
        //iRelinkTime.setColumns(5);
        iRelinkTime = new JLabel(FIVE);
        iRelinkTime.setBorder(blackline);
        //iRelinkTime.setColumns(5);

        s.setBorder(topBorder);
        s.add(iRelinkTime);
        pEnd.add(s);

        s = new JPanel();
        s.setLayout(new FlowLayout());
        label = new JLabel(RELINKNUC);
        s.add(label);
        iRelinkNuc = new JLabel(TWELVE);
        iRelinkNuc.setBorder(blackline);


        //iRelinkNuc.setColumns(12);
        //iRelinkNuc.setText(cell.getName());
        s.add(iRelinkNuc);
        s.setBorder(botBorder);
        pEnd.add(s);

        iRelinkButton = new JButton(SETLATECELL);
        iRelinkButton.addActionListener(this);
        s = new JPanel();
        s.setLayout(new GridLayout(1,1));
        s.setBorder(topBorder);
        s.add(iRelinkButton);

        pEnd.add(s);
        pWhole.add(pEnd);


        JPanel xp = new JPanel();
        xp.setBorder(blackline);
        xp.setLayout(new BoxLayout(xp, BoxLayout.PAGE_AXIS));
        Border b = BorderFactory.createEmptyBorder(10,0,10,0);

        s = new JPanel();
        s.setLayout(new GridLayout(3,1));
        s.add(new JLabel(""));

        iApplyAndRebuild = new JButton(APPLYANDREBUILD);
        iApplyAndRebuild.addActionListener(this);
        s.add(iApplyAndRebuild);

        iApplyOnly = new JButton(APPLYONLY);
        iApplyOnly.addActionListener(this);
        s.add(iApplyOnly);
        pWhole.add(s);

        pWhole.setOpaque(true); //content panes must be opaque
        dialog.setContentPane(pWhole);
        dialog.setSize(new Dimension(200, 400));
        dialog.setLocationRelativeTo(owner);
        dialog.setVisible(true);

        addKeyListener(new MyKeyListener());
        //setFocusableFalse();
        setDefaultButtonBehavior(iLinkButton, "SPACE");
        JButton [] jba = new JButton[4];
        jba[0] = iLinkButton;
        jba[1] = iRelinkButton;
        jba[2] = iApplyOnly;
        jba[3] = iApplyAndRebuild;
        //setKeyBehavior(jba, "ENTER");
        //iLinkButton.requestFocus();
    }

    private void setFocusableFalse() {
        iLinkTime.setFocusable(false);
        iLinkNuc.setFocusable(false);
        iLinkButton.setFocusable(false);
        iLinkRootButton.setFocusable(false);
        iRelinkTime.setFocusable(false);
        iRelinkNuc.setFocusable(false);
        iRelinkButton.setFocusable(false);
        iApplyAndRebuild.setFocusable(false);
        iApplyOnly.setFocusable(false);

    }

    void setDefaultButtonBehavior(JButton jb, String key) {
        String s = key;
        Action home = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	println("focus to AceTree");
            	iAceTree.requestFocus();
            }
        };
        //jb.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(s), "pressed");
        jb.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(KeyEvent.VK_END, 0, false), "pressed");
        jb.getActionMap().put("pressed", home );


    }

    void setKeyBehavior(JButton [] jba, String key) {
        String s = key;
        Action home = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	Component compFocusOwner = KeyboardFocusManager.getCurrentKeyboardFocusManager().getFocusOwner();
            	if (compFocusOwner instanceof JButton) {
            		println("its a button");
            		((JButton)compFocusOwner).doClick();
            	}
            }
        };
        for (int i=0; i < jba.length; i++) {
        	jba[i].getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(s), s);
        	jba[i].getActionMap().put(s, home );
        }


    }


    public void actionPerformed(ActionEvent e) {
        //System.out.println("NucRelinkDialog.actionPerformed");
        Object o = e.getSource();
        String cmd = e.getActionCommand();
        if (cmd.equals(SETEARLYCELL)) {
            int time = iAceTree.getImageTime() + iAceTree.getTimeInc();
            iLinkTime.setText(String.valueOf(time));
            iLinkNuc.setText(iAceTree.getCurrentCell().getName());
        } else if (o == iLinkRootButton) {
            iLinkTime.setText(String.valueOf((iAceTree.getNucleiMgr()).getStartingIndex()));
            iLinkNuc.setText(AceTree.ROOTNAME);
        } else if (cmd.equals(SETLATECELL)) {
            int time = iAceTree.getImageTime() + iAceTree.getTimeInc();
            iRelinkTime.setText(String.valueOf(time));
            iRelinkNuc.setText(iAceTree.getCurrentCell().getName());

        } else if (cmd.equals(APPLYANDREBUILD) || cmd.equals(APPLYONLY)) {
            int endTime;
            try {
                endTime = Integer.parseInt(iRelinkTime.getText());
            } catch(NumberFormatException nfe) {
                showMessage("invalid relink time, aborting");
                return;
            }
            int strTime;
            try {
                strTime = Integer.parseInt(iLinkTime.getText());
            } catch(NumberFormatException nfe) {
                showMessage("invalid link time, aborting");
                return;
            }

            if (strTime > 1 && endTime <= strTime) {
                showMessage("endTime is not greater than start time, aborting");
                return;
            }

            String endCellName = iRelinkNuc.getText();
            String strCellName = iLinkNuc.getText();
            boolean b = checkCellValidities(endCellName, endTime, strCellName, strTime);
            //println("actionPerformed: checkCellValidities returned: " + b);

            if (!b) return;
            /*
            boolean endGood = NucUtils.isValidCell(endCellName, endTime);
            boolean strGood = NucUtils.isValidCell(strCellName, strTime);
            if (!(endGood && strGood)) {
                String s0 = "";
                String s1 = "";
                s0 = "invalid cell: ";
                if (!endGood) s0 = "invalid cell: " + endCellName + CS + endTime;
                if (!strGood) s1 = "invalid cell: " + strCellName + CS + strTime;
                JOptionPane pane = new JOptionPane(s0 + "\n" + s1);
                JDialog dialog = pane.createDialog(iAceTree, "About AceTree");
                dialog.setModal(true);
                dialog.show();
                //dispose();
                return;
            }
            */

            StringBuffer sb = new StringBuffer("RELINKING: ");
            sb.append(endTime);
            sb.append(CS + endCellName);
            sb.append(iNucleiMgr.getIndex(endCellName, endTime));
            sb.append(CS + strTime);
            sb.append(CS + strCellName);
            sb.append(iNucleiMgr.getIndex(strCellName, strTime));
            //System.out.println(sb.toString());
            iEditLog.append(sb.toString());

            iNucleiMgr.makeBackupNucleiRecord();
            createAndAddCells(endCellName, endTime, strCellName, strTime);
            //System.out.println("returned from createAndAddCells");
            if (cmd.equals(APPLYANDREBUILD)) {
                //println("\n\nNucRelinkDialog.actionPerformed: applyAndRebuild");
                iAceTree.clearTree();
                iAceTree.buildTree(true);
                AncesTree ances = iAceTree.getAncesTree();
                Hashtable h = ances.getCellsByName();
                Cell c = (Cell)h.get(strCellName);
                iAceTree.setStartingCell(c, strTime);
                iEditLog.setModified(true);
                //dispose();
                iRelinkNuc.setText(FIVE);
                iRelinkTime.setText(TWELVE);
                iLinkNuc.setText(strCellName);
                char x = endCellName.charAt(0);
                if (x != '_' && x != 'N') iLinkNuc.setText(endCellName);
                iLinkTime.setText(String.valueOf(endTime));
            }

        }
        iAceTree.requestFocus();
    }

    private boolean checkCellValidities(String endCellName, int endTime, String strCellName, int strTime) {
        if (strCellName.equals(AceTree.ROOTNAME)) return true;
        Nucleus nEnd = iNucleiMgr.getCurrentCellData(endCellName, endTime);
        Nucleus nStr = iNucleiMgr.getCurrentCellData(strCellName, strTime);
        if (nEnd == null || nStr == null) {
            String s  = "";
            String s0 = "";
            String s1 = "";
            //s0 = "invalid cell: ";
            if (nEnd == null) s0 = "invalid cell: " + endCellName + CS + endTime + NL;
            if (nStr == null) s1 = "invalid cell: " + strCellName + CS + strTime + NL;
            showMessage(s0 + s1);
            return false;
        }
        if (nStr.successor2 > 0) {
            String s = "Cell " + strCellName + " already has 2 successors\n";
            s = s + "cannot complete relink.";
            showMessage(s);
            return false;
        }
        return true;
    }

    private void showMessage(String s) {
        JOptionPane pane = new JOptionPane(s);
        JDialog dialog = pane.createDialog(iAceTree, "About AceTree");
        dialog.setModal(true);
        dialog.show();
    }

    public static void createAndAddCells(String endCellName, int endTime, String strCellName, int strTime) {
        // access nucleus record of end and start cells
        println("createAndAddCells, " + endCellName + CS + endTime + CS + strCellName + CS + strTime);
        Nucleus nEnd = getNucleus(endCellName, endTime);
        if (strCellName.equals(AceTree.ROOTNAME)) {
            nEnd.predecessor = Nucleus.NILLI;
            return;
        }
        //System.out.println("endCell: " + endCellName + CS + endTime);
        //System.out.println("nEnd: " + nEnd);
        //System.out.println("startCell: " + strCellName + CS + strTime);
        Nucleus nStr = getNucleus(strCellName, strTime);
        //System.out.println("actionPerformed: nStr: " + nStr);
        Vector nuclei_record = iNucleiMgr.getNucleiRecord();
        Vector nucleiAdd = null;
        Nucleus n = nStr;
        int predecessor = nStr.index;
        for (int k = strTime + 1; k < endTime; k++) {
            nucleiAdd = (Vector)nuclei_record.elementAt(k - 1);
            n = interpolateNucleus(nEnd, nStr, endTime, strTime, k);
            n.index = nucleiAdd.size() + 1;
            //n.snindex = n.index;
            n.predecessor = predecessor;
            predecessor = n.index;
            //System.out.println("adding: " + n);
            //iEditLog.append("adding: " + n.toString());
            nucleiAdd.add(n);
        }
        nEnd.predecessor = n.index;
        //System.out.print("nEnd: " + nEnd);
    }

    private static Nucleus getNucleus(String name, int time) {
        System.out.println("getNucleus, seeming: " + name + CS + time);
        Nucleus nRtn = null;
        Nucleus n = null;
        Vector nuclei_record = iNucleiMgr.getNucleiRecord();
        Vector nuclei = (Vector)nuclei_record.elementAt(time - 1);
        for (int j=0; j < nuclei.size(); j++) {
            n = (Nucleus)nuclei.elementAt(j);

            if (n.status > 0 && n.identity.equals(name)) {
                nRtn = n;
                break;
            }
        }
        return nRtn;
    }

    private static Nucleus interpolateNucleus(Nucleus nEnd, Nucleus nStr, int endTime, int strTime, int midTime) {
        Nucleus n = nStr.copy();
        int deltaT = endTime - strTime;
        int deltaM = midTime - strTime;
        n.x = (nEnd.x - nStr.x)*deltaM/deltaT + nStr.x;
        n.y = (nEnd.y - nStr.y)*deltaM/deltaT + nStr.y;
        n.z = (nEnd.z - nStr.z)*deltaM/deltaT + nStr.z;
        n.size = (nEnd.size - nStr.size)*deltaM/deltaT + nStr.size;
        return n;
    }

    public class MyKeyListener extends KeyAdapter{
        public void keyPressed(KeyEvent ke){
          char i = ke.getKeyChar();
          String str = Character.toString(i);
          println("MyKeyListener, " + str);
          //iAceTree.requestFocus();
        }
      }

    private class WindowEventHandler extends WindowAdapter {
        public void windowClosing(WindowEvent e) {
        	iAceTree.iNucRelinkDialog = null;
            println("NucRelinkDialog.windowclosing");
        }

    }

    public final static String
        TITLE = "Nuclei Editor"
       ,RELINKTIME = "Later time"
       ,RELINKNUC = "Cell name"
       ,LINKTIME = "Earlier time"
       ,LINKNUC = "Cell name"
       ,DOIT = "Apply"
       ,CS = ", "
       ,NL = "\n"
       ,APPLYANDREBUILD = "apply/rebuild"
       ,APPLYONLY = "apply only"
       ,SETCURRENTCELL = "set current cell"
       ,SETEARLYCELL = "set early cell"
       ,SETLATECELL = "set late cell"
       ,SETROOTCELL = "set root cell"
       ,LATER = "link cell at later time"
       ,EARLIER = "to cell at earlier time"
       ,FIVE = "     "
       ,TWELVE = "            "
       ;

       private static void println(String s) {System.out.println(s);}
       private static void print(String s) {System.out.print(s);}
       private static final String TAB = "\t";
       private static final DecimalFormat DF0 = new DecimalFormat("####");
       private static final DecimalFormat DF1 = new DecimalFormat("####.#");
       private static final DecimalFormat DF4 = new DecimalFormat("####.####");
       private static String fmt4(double d) {return DF4.format(d);}
       private static String fmt1(double d) {return DF1.format(d);}
       private static String fmt0(double d) {return DF0.format(d);}
}
