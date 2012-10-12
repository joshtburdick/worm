package org.rhwlab.nucedit;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.KeyStroke;

import org.rhwlab.acetree.AceTree;
import org.rhwlab.acetree.NucUtils;
//import org.rhwlab.image.EditImage3;
import org.rhwlab.image.ImageWindow;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.Cell;
import org.rhwlab.utils.C;

public class AddOneDialog extends JDialog implements ActionListener, WindowFocusListener {
    ImageWindow             iParent;
    private AceTree         iAceTree;
    int                     iImageTime;
    int                     iTimeInc;
    int                     iImagePlane;
    int                     iPlaneInc;
    int                     iPrevTime;
    Cell                    iCurrentCell;
    Nucleus                 iNucleus;
    protected int           iNucSize;
    private Cell            iCurrentCellSave;
    private Nucleus         iNucleusCopy;
    private Nucleus         iNucleusActive;
    private boolean         iTimeChanged;
    private int             iTimeSave;



    private JRadioButton    iAdjust;
    private JRadioButton    iAdd;
    //private JRadioButton    iAddSeries;

    private JButton         iRelink;
    private JButton         iKillCells;
    private JButton         iRebuildAndRename;
    private JButton         iRebuildOnly;
    private JButton         iUp;
    private JButton         iDown;
    private JButton         iLeft;
    private JButton         iRight;
    private JButton         iBig;
    private JButton         iSmall;
    private JButton         iIncZ;
    private JButton         iDecZ;
    private JButton         iUndo;
    private JButton         iTest;
    private JTextField      iName;
    private JTextField      iForceName;
    private JTextField      iX;
    private JTextField      iY;
    private JTextField      iZ;
    private JTextField      iD;
    private JButton         iSetN;
    private JButton         iForce;
    private JButton         iSetX;
    private JButton         iSetY;
    private JButton         iSetZ;
    private JButton         iSetD;
    private JButton			iDefault;


    public AddOneDialog(AceTree aceTree, Frame owner, boolean modal, Cell cell, int time) {
        super(owner, modal);
        addWindowListener(new WindowEventHandler());
        //println("AddOneDialog, " + cell.getName() + CS + time);
        setTitle("Add single nuclei");
        iAceTree = aceTree;
        iNucSize = 50;
        JDialog dialog = this;
        JPanel pWhole = new JPanel();
        iParent = (ImageWindow)owner;
        pWhole.setOpaque(true); //content panes must be opaque
        fillControlPanel(pWhole);

        iTimeSave = -1;
        updateCurrentInfo(true);
        updateTextFields();
        //iCurrentCellSave = iCurrentCell;
        //System.out.println("EditImage: " + iNucleus);
        if (iNucleus == null) iNucleus = new Nucleus();
        iNucleusCopy = iNucleus.copy();
        iNucleusActive = iNucleus;

        iTimeSave = iImageTime + iTimeInc;
        iCurrentCellSave = iCurrentCell;

        dialog.setContentPane(pWhole);
        dialog.setSize(new Dimension(310, 450));
        dialog.setLocationRelativeTo(owner);
        dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
        dialog.setVisible(true);
        setKeyboardActions();

        addWindowFocusListener(this);



    }

    private void setKeyboardActions() {
        Action home = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	println("setKeyboardActions");
            	iAceTree.requestFocus();
            }
        };
        KeyStroke key = null;
        key = KeyStroke.getKeyStroke("F2");
        key = KeyStroke.getKeyStroke(KeyEvent.VK_SPACE, 0, false);

        iDefault.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(key, "pressed");
        iDefault.getActionMap().put("pressed", home );
        iDefault.requestFocus();

    }

    private class WindowEventHandler extends WindowAdapter {
        public void windowClosing(WindowEvent e) {
        	iAceTree.iAddOneDialog = null;
            println("AddOneDialog.windowclosing");
        }

    }



    /* (non-Javadoc)
     * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
     */
    public void actionPerformed(ActionEvent e) {
        Object o = e.getSource();
        String cmd = e.getActionCommand();
        //println("EIDialog2.actionPerformed: " + o);
        if (o == iAdd) return;
        if (o == iRelink) {
            iAceTree.relinkNucleus();
            iAdjust.setSelected(true);
            iParent.refreshDisplay(null);
        } else if (o == iKillCells) {
                iAceTree.killCells();
                iAdjust.setSelected(true);
                iParent.refreshDisplay(null);

        } else if (cmd.equals(REBUILDANDRENAME)) {
                updateCurrentInfo(false);
                int time = iImageTime + iTimeInc;
                Cell c = iCurrentCell;
                iAceTree.clearTree();
                iAceTree.buildTree(true);
                //iEditLog.setModified(true);
                //System.out.println("actionPerformed: " + c + C.CS + time);
                if (c != null) iAceTree.setStartingCell(c, time);
                iAdjust.setSelected(true);

        } else if (o == iRebuildOnly) {
                updateCurrentInfo(false);
                int time = iImageTime + iTimeInc;
                Cell c = iCurrentCell;
                iAceTree.clearTree();
                iAceTree.buildTree(false);
                //iEditLog.setModified(true);
                System.out.println("actionPerformed: " + c + C.CS + time);
                if (c != null) iAceTree.setStartingCell(c, time);
                iAdjust.setSelected(true);

        } else {
            //println("AddOneDialog.actionPerformed:2 iCurrentCell: " + iCurrentCell);
            updateCurrentInfo(false);
            setKeypadEnabled(true);
            //println("EIDialog2.actionPerformed:3 iCurrentCell: " + iCurrentCell);
            Nucleus n = ImageWindow.cNucleiMgr.getNucleusFromHashkey(iCurrentCell.getHashKey(), iImageTime + iTimeInc);
            if (cmd.equals(UP)) n.y--;
            else if (cmd.equals(DOWN)) n.y++;
            else if (cmd.equals(LEFT)) n.x--;
            else if (cmd.equals(RIGHT)) n.x++;
            else if (cmd.equals(BIG)) {
                n.size += 2;
                iNucSize = n.size;
            }
            else if (cmd.equals(SMALL)) {
                n.size -= 2;
                iNucSize = n.size;
            } else if (cmd.equals(INCZ)) {
            	n.z += ZINCREMENT;
            	iZ.setText(String.valueOf(n.z));
            } else if (cmd.equals(DECZ)) {
            	n.z -= ZINCREMENT;
            	iZ.setText(String.valueOf(n.z));
            }
            else if (o == iX || o == iSetX) n.x = Integer.parseInt(iX.getText());
            else if (o == iY || o == iSetY) n.y = Integer.parseInt(iY.getText());
            else if (o == iZ || o == iSetZ) n.z = Float.parseFloat(iZ.getText());
            else if (o == iD || o == iSetD) n.size = Integer.parseInt(iD.getText());
            else if (o == iName) {
                String oldName = n.identity;
                n.identity = iName.getText();
                iCurrentCell.setName(n.identity);
                iParent.updateCellAnnotation(iCurrentCell, oldName, iImageTime + iTimeInc);
            } else if (o == iForce) {
                n.assignedID = iForceName.getText();
                n.identity = n.assignedID;
                iName.setText(n.assignedID);
                iForceName.setText("");
                // we need to change the hashkey in the nucleus object and cell object
                // assume we know the iCurrentCell at this point
                int time = iImageTime + iTimeInc;
                String hashKey = NucUtils.makeHashKey(time, n);
                //System.out.println("addCell: " + hashKey);
                n.setHashKey(hashKey);
                iCurrentCell.setHashKey(hashKey);
                iCurrentCell.setName(n.identity);


            } else if (o == iUndo) {
                undo();
                iParent.refreshDisplay(null);
            }
        //} else {
            //iAddSeriesAction = true;
            //setKeypadEnabled(false);
            else if (o == iUndo && iAdd.isSelected()) {
                undo();
                iParent.refreshDisplay(null);
            }
            iParent.refreshDisplay(null);
        }

    }

    private void undo() {
        //System.out.println("test undo");
        if (iAdd.isSelected()) {
            iNucleus.status = Nucleus.NILLI;
            //Log editLog = iAceTree.getEditLog();
            //editLog.append("CELLADDITIONCANCELED at " + (iImageTime + iTimeInc));
            //editLog.append("now: " + iNucleus);
        }
        //System.out.println("test2: " + iNucleus);
        //iNucleus = iNucleusCopy.copy();
        if(iAdjust.isSelected()) {
            iNucleus.x = iNucleusCopy.x;
            iNucleus.y = iNucleusCopy.y;
            iNucleus.z = iNucleusCopy.z;
            iNucleus.size = iNucleusCopy.size;
            iNucleus.identity = iNucleusCopy.identity;
        }
        //System.out.println("test3: " + iNucleus);
        //refreshDisplay(null);

    }

    public void processMouseEvent(MouseEvent e) {
        //println("AddOneDialog.processMouseEvent: " + e);
        int button = e.getButton();
        if (button == 1) {
            if (iAdd.isSelected()) {
                addCell(e.getX(), e.getY(), false);
            }
        } else if (button == 3) {
            updateCurrentInfo(false);
            Nucleus n = ImageWindow.cNucleiMgr.findClosestNucleus(e.getX(), e.getY(), iImagePlane + iPlaneInc, iImageTime + iTimeInc);
            if (n == null) {
                System.out.println("cant find closest nucleus");
                return;
            }
            Cell c = (Cell)iAceTree.getAncesTree().getCells().get(n.hashKey);

            //System.out.println("mouseClicked1: " + c + C.CS + iCurrentCell
            //        + C.CS + iImagePlane + C.CS + iPlaneInc);
            iAceTree.setCurrentCell(c, iImageTime + iTimeInc, AceTree.RIGHTCLICKONEDITIMAGE);
            iNucleus = n;
            updateTextFields();

        } else if (button == 2) {
            addCell(e.getX(), e.getY(), false);

        }
    }

    protected void addCell(int x, int y, boolean continuation) {
        //System.out.println("addCell: " + x + C.CS + y);
        updateCurrentInfo(false);
        int time = iImageTime + iTimeInc;
        Vector nuclei = (Vector)ImageWindow.cNucleiMgr.getNucleiRecord().elementAt(time - 1);
        Nucleus n = new Nucleus();
        n.index = nuclei.size() + 1;
        String hashKey = NucUtils.makeHashKey(time, n);
        //System.out.println("addCell: " + hashKey);
        n.setHashKey(hashKey);
        n.status = 1;
        n.x = x;
        n.y = y;
        n.z = iImagePlane + iPlaneInc;
        n.size = iNucSize;
        n.identity = "_" + hashKey;
        n.predecessor = -1;
        n.successor1 = -1;
        n.successor2 = -1;
        nuclei.add(n);
        iNucleus = n;

        Cell c = new Cell(n.identity, time);
        c.setHashKey(hashKey);
        iAceTree.getAncesTree().getCells().put(hashKey, c);
        iAceTree.setShowCentroids(true);
        iAceTree.setShowAnnotations(true);

        c.setParameters(time, time, n);
        Cell root = iAceTree.getAncesTree().getRoot();
        c.setParent(root);
        //iAceTree.setCurrentCell(c, time, AceTree.CONTROLCALLBACK);
//        iAceTree.setCurrentCell(c, iImageTime + iTimeInc, AceTree.RIGHTCLICKONEDITIMAGE);
        iAceTree.setCurrentCell(c, time, AceTree.RIGHTCLICKONEDITIMAGE);
        iParent.addAnnotation(x, y, true);



        iAdjust.setSelected(true);
        setKeypadEnabled(true);
        iName.setText(n.identity);


        iParent.refreshDisplay(null);
        //System.out.println("addCell: " + iCurrentCell);
        //iAceTree.clearTree();
        //iAceTree.buildTree(true);
    }


    protected void updateCurrentInfo(boolean detectChange) {
        //System.out.println("EditImage.updateCurrentInfo called: " + new GregorianCalendar().getTime());
        iImageTime = iAceTree.getImageTime();
        iImagePlane = iAceTree.getImagePlane();
        iTimeInc = iAceTree.getTimeInc();
        iPlaneInc = iAceTree.getPlaneInc();
        iCurrentCell = iAceTree.getCurrentCell();
        Vector nuclei = ImageWindow.cNucleiMgr.getNuclei(iImageTime + iTimeInc - 1);
        iNucleus = NucUtils.getCurrentCellNucleus(nuclei, iCurrentCell);
        //println("updateCurrentInfo: " + iCurrentCell + CS + iNucleus);
        //if (!detectChange) return;
        //System.out.println("updateCurrentInfo: detect change steps implemented");
    }

    protected void updateTextFields() {
        updateCurrentInfo(false);
        //System.out.println("updateTextFields: " + iName + C.CS + iNucleus);
        if (iNucleus == null) return; //workaround
        iName.setText(iNucleus.identity);
        iForceName.setText(iNucleus.assignedID);
        iX.setText(String.valueOf(iNucleus.x));
        iY.setText(String.valueOf(iNucleus.y));
        iZ.setText(String.valueOf(iNucleus.z));
        iD.setText(String.valueOf(iNucleus.size));
    }



    protected void fillControlPanel(JPanel pp) {
        addChoices(pp);
        addKeypad(pp);
        addTextFields(pp);

    }
    protected void addKeypad(JPanel mp) {
        JPanel p = new JPanel();
        iLeft = new JButton(LEFT);
        iRight = new JButton(RIGHT);
        iUp = new JButton(UP);
        iDown = new JButton(DOWN);
        //iUndo = new JButton(UNDO);
        iTest = new JButton(TEST);
        iBig = new JButton(BIG);
        iSmall = new JButton(SMALL);
        iIncZ = new JButton(INCZ);
        iDecZ = new JButton(DECZ);

        //iShowC = new JButton(SHOWC);
        iLeft.addActionListener(this);
        iRight.addActionListener(this);
        iUp.addActionListener(this);
        iDown.addActionListener(this);
        //iUndo.addActionListener(this);
        iTest.addActionListener(this);
        iBig.addActionListener(this);
        iSmall.addActionListener(this);
        iIncZ.addActionListener(this);
        iDecZ.addActionListener(this);
        //iHome.addActionListener(this);
        p.setLayout(new GridLayout(3,3));
        p.setBorder(BorderFactory.createLineBorder(Color.white));
        p.add(iBig);
        p.add(iUp);
        p.add(iSmall);
        p.add(iLeft);
        p.add(new JButton());
        p.add(iRight);
        p.add(iIncZ);
        p.add(iDown);
        //p.add(iUndo);
        p.add(iDecZ);
        mp.add(p);
        //setKeypadEnabled(false);
        iDefault = iSmall;

    }

    protected void setKeypadEnabled(boolean b) {
        iUp.setEnabled(b);
        iDown.setEnabled(b);
        iLeft.setEnabled(b);
        iRight.setEnabled(b);
        iBig.setEnabled(b);
        iSmall.setEnabled(b);
        iName.setEnabled(b);
        iX.setEnabled(b);
        iY.setEnabled(b);
        iZ.setEnabled(b);
        iD.setEnabled(b);
    }

    protected void addChoices(JPanel mp) {
        JRadioButton rb = null;
        ButtonGroup bg = new ButtonGroup();
        JPanel rp = new JPanel(new GridLayout(0, 1));
        //rp.setBorder(BorderFactory.createEmptyBorder(20,20,20,20));
        //rp.setBorder(BorderFactory.createLineBorder(Color.black));

        iAdd = new JRadioButton(ADD);
        iAdd.addActionListener(this);
        bg.add(iAdd);
        rp.add(iAdd);

        //iAddSeries = new JRadioButton(ADDSERIES);
        //iAddSeries.addActionListener(this);
        //bg.add(iAddSeries);
        //rp.add(iAddSeries);

        iAdjust = new JRadioButton(ADJUST);
        iAdjust.addActionListener(this);
        iAdjust.setSelected(true);
        bg.add(iAdjust);
        rp.add(iAdjust);

        mp.add(rp);

        rp = new JPanel(new GridLayout(0, 1));
        //rp.setBorder(BorderFactory.createEmptyBorder(20,20,20,20));
        rp.setBorder(BorderFactory.createLineBorder(Color.white));
        rp.setPreferredSize(new Dimension(350, 100));

        iRelink= new JButton(RELINK);
        iRelink.addActionListener(this);
        //bg.add(iRelink);
        rp.add(iRelink);

        iKillCells= new JButton(KILLCELLS);
        iKillCells.addActionListener(this);
        //bg.add(iKillCells);
        rp.add(iKillCells);

        iRebuildOnly= new JButton(REBUILDONLY);
        iRebuildOnly.addActionListener(this);
        //bg.add(iRebuildOnly);
        rp.add(iRebuildOnly);

        iRebuildAndRename= new JButton(REBUILDANDRENAME);
        iRebuildAndRename.addActionListener(this);
        //bg.add(iRebuildAndRename);
        rp.add(iRebuildAndRename);

        mp.add(rp);
    }

    protected void addTextFields(JPanel mp) {
        JPanel p = new JPanel();
        //p.setBorder(BorderFactory.createEmptyBorder(20,20,20,20));
        p.setLayout(new GridLayout(0,3));
        iName = new JTextField(NAME, 8);
        iForceName = new JTextField("");
        iX = new JTextField(X);
        iY = new JTextField(Y);
        iZ = new JTextField(Z);
        iD = new JTextField(D);
        iName.addActionListener(this);
        iX.addActionListener(this);
        iY.addActionListener(this);
        iZ.addActionListener(this);
        iD.addActionListener(this);
        JLabel l = null;
        l = new JLabel(NAME);
        p.add(l);
        p.add(iName);
        iSetN = new JButton(SETN);
        iSetN.addActionListener(this);
        p.add(iSetN);
        l = new JLabel(FORCENAME);
        p.add(l);
        p.add(iForceName);
        iForce = new JButton(FORCE);
        iForce.addActionListener(this);
        p.add(iForce);
        l = new JLabel(X);
        p.add(l);
        p.add(iX);
        iSetX = new JButton(SETX);
        iSetX.addActionListener(this);
        p.add(iSetX);
        l = new JLabel(Y);
        p.add(l);
        p.add(iY);
        iSetY = new JButton(SETY);
        iSetY.addActionListener(this);
        p.add(iSetY);
        l = new JLabel(Z);
        p.add(l);
        p.add(iZ);
        iSetZ = new JButton(SETZ);
        iSetZ.addActionListener(this);
        p.add(iSetZ);
        l = new JLabel(D);
        p.add(l);
        p.add(iD);
        iSetD = new JButton(SETD);
        iSetD.addActionListener(this);
        p.add(iSetD);
        mp.add(p);

    }

    public void processWindowEvent(WindowEvent e){
        //println("processWindowEvent: " + e);
        int id = e.getID();
        if (id == WindowEvent.WINDOW_CLOSING) {
            iParent.parentNotifyDialogClosing(this);
            iAceTree.iAddOneDialog = null;
            dispose();
            //System.exit(0);
        }

    }



    public static final String
    UP = "UP"
   ,DOWN = "DOWN"
   ,LEFT = "LEFT"
   ,RIGHT = "RIGHT"
   ,TEST = "TEST"
   ,UNDO = "UNDO"
   ,BIG = "BIG"
   ,SMALL = "SMALL"
   ,INCZ = "INCZ"
   ,DECZ = "DECZ"
   ,ADJUST = "adjustCell"
   ,ADD = "addCell"
   ,ADDSERIES = "addCellSeries"
   ,REBUILDANDRENAME =  "rebuildAndRename"
   ,RELINK = "relink"
   ,KILLCELLS = "killCells"
   ,REBUILDONLY = "rebuildOnly"
   ,NAME = "Name"
   ,FORCENAME = "ForceName"
   ,X = "X"
   ,Y = "Y"
   ,Z = "Z"
   ,D = "D"
   ,FORCE = "Force"
   ,SETN = "SetN"
   ,SETX = "SetX"
   ,SETY = "SetY"
   ,SETZ = "SetZ"
   ,SETD = "SetD"

   ;

    static final double ZINCREMENT = 0.5;



    public static void main(String[] args) {
        //EIDialog2 eid2 = new EIDialog2(null, null, false, null, 0);

    }
    protected static void println(String s) {System.out.println(s);}
    protected static final String CS = ", ";


	public void windowGainedFocus(WindowEvent e) {
		println("AddOneDialog.windowGainedFocus, ");
		iParent.iDialog = this;

	}

	public void windowLostFocus(WindowEvent e) {
		println("AddOneDialog.windowLostFocus, ");

	}


}
