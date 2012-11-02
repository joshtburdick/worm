/*
 * Copyright 2005 University of Washington Genome Sciences
 * All rights reserved
 */
package org.rhwlab.acetree;

import org.murrlab.image.MedianSmoother;

import org.rhwlab.analyze.DeathAndDivisionLog;
import org.rhwlab.help.AceTreeHelp;
import org.rhwlab.help.TestWindow;
import org.rhwlab.image.CellMovementImage;
//import org.rhwlab.image.EditImage;
//import org.rhwlab.image.EditImage3;
import org.rhwlab.image.Image3D;
import org.rhwlab.image.Image3D2;
import org.rhwlab.image.Image3D2Z;
import org.rhwlab.image.ImageAllCentroids;
import org.rhwlab.image.ZipImage;
import org.rhwlab.image.ImageWindow;
//import org.rhwlab.nucedit.AddNucToRoot;
import org.rhwlab.manifest.ManifestX;
import org.rhwlab.nucedit.EditLog;
import org.rhwlab.nucedit.KillCellsDialog;
//import org.rhwlab.nucedit.NucAddDialog;
import org.rhwlab.nucedit.AddOneDialog;
import org.rhwlab.nucedit.DeathsAdjacencies;
//import org.rhwlab.nucedit.EIDialog2;
import org.rhwlab.nucedit.EditTraverse;
import org.rhwlab.nucedit.Juvenesence;
import org.rhwlab.nucedit.KillDeepNucsDialog;
import org.rhwlab.nucedit.Lazarus;
import org.rhwlab.nucedit.NucEditDialog;
import org.rhwlab.nucedit.NucRelinkDialog;
import org.rhwlab.nucedit.Orientation;
import org.rhwlab.nucedit.SetEndTimeDialog;
import org.rhwlab.nucedit.Siamese;
import org.rhwlab.nucedit.Zafer1;
import org.rhwlab.snight.Config;
import org.rhwlab.snight.Identity;
import org.rhwlab.snight.NucZipper;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.snight.Parameters;
import org.rhwlab.tree.AncesTree;
import org.rhwlab.tree.CanonicalTree;
import org.rhwlab.tree.Cell;
import org.rhwlab.tree.CellData;
import org.rhwlab.tree.Newick;
import org.rhwlab.tree.SulstonTree;
import org.rhwlab.tree.VTree;
import org.rhwlab.utils.AuxFrame;
import org.rhwlab.utils.C;
import org.rhwlab.utils.EUtils;
import org.rhwlab.utils.Log;

//import org.rhwlab.utils.*;
//import org.rhwlab.snight.*;
import ij.IJ;
import ij.ImagePlus;

import java.awt.AWTEvent;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Event;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.KeyboardFocusManager;
import java.awt.Toolkit;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.Calendar;
import java.util.Collections;
import java.util.Enumeration;
import java.util.GregorianCalendar;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ActionMap;
import javax.swing.BoxLayout;
import javax.swing.Icon;
import javax.swing.InputMap;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.JTree;
import javax.swing.KeyStroke;
import javax.swing.event.MouseInputAdapter;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;
/**
 *
 * The main class for the embryo data and analysis user interface.
 * <p>
 * The major supporting classes of the application are:
 * <br><b>AncesTree</b> - builds the tree from the "nuclei" analysis files
 * <br><b>Cell</b> - the node cell of the tree from DefaultMutableTreeNode
 * <br><b>ZipImage</b> - handles display of the selected image
 * <br><b>ZipNuclei</b> - opens the zip data file and enables access to the analysis
 * files therein. Also supplies an instance of class ZipFile used elsewhere.
 * <br><b>NucUtils</b> - class holding many static utility functions closely related
 * </p><p>
 * The following additional supporting classes are involved:
 * <br>AnnotInfo - small data structure
 * <br>ControlCallback - interfaced defined so UI elements can be in a separate class
 * <br>EUtils - Designed to hold simple static utility functions not specific to this app.
 * <br>InputCtrl - a control now used to create on panel in the Our_Tree3 UI
 * <br>MouseHandler - captures mouse movements on the image window
 * to the app.
 * <br>SpringUtilities - class from java tutorial needed by InputCtrl
 *
 *
 * @author biowolp
 * @version 1.0 January 18, 2005
 */
public class AceTree extends JPanel
            implements /*TreeSelectionListener, PlugIn, */
            ActionListener, ControlCallback, Runnable {

    protected static AceTree  iAceTree;

    private Hashtable   iNucleiMgrHash;

    private String      iConfigFileName;
    private JTree       iTree;
    private String []   iImgSuffix;
    private Cell        iRoot;
    private JTextPane   iText;
    private JTextPane   iText2;
    private JTextPane   iText3;
    private JFrame      iMainFrame;
    private Vector      iTempV;
    private String      iFilePath;
    private boolean     iRootEstablished;
    private int         iImageTime;
    private int         iImagePlane;
    public  ImageWindow iImgWin;
    private boolean     iImgWinSet;
    private NucleiMgr   iNucleiMgr;
    //private Parameters  iParameters;
    public AncesTree   iAncesTree;
    private String      iOrientation;

    public AceMenuBar  iAceMenuBar;
    private EditLog     iEditLog;
    private Log         iDDLog;
    private Log         iDLog;
    //private boolean     iEditLogInitialized;

    private WindowEventHandler  iWinEvtHandler;
    private FileInputStream     iFstream;
    private BufferedReader      iBReader;

    private JButton     iCopy;
    private JButton     iShow;
    private JButton     iClear;
    private JButton     iNext;
    private JButton     iPrev;
    private JButton     iUp;
    private JButton     iDown;
    private JButton     iHome;
    private JButton		iDefault;
    private JButton     iShowC;
    private JButton     iTrack;
    private JButton     iSister;
    private JButton     iColorToggle;
    private JLabel      iLabel;
    private int         iTimeInc;
    private int         iPlaneInc;
    private Cell        iCurrentCell;
    private boolean     iCurrentCellPresent;
    private int         iCurrentCellXloc;
    private int         iCurrentCellYloc;
    private float       iCurrentCellZloc;
    private InputCtrl   iInputCtrl;
    //private String      iAxis;

    // key run parameters
    private String      iZipFileName; // a full path to the zip with nuclei and parameters
    private String      iZipNucDir;   // subdirectory in above zip for nuclei
    private String      iZipTifFilePath; // a full path to the zip file with tifs in it
    private String      iTifPrefix;      // leading part of image file names also parameters file
    private int         iStartingIndex;
    private int         iEndingIndex;
    private int         iNamingMethod;
    private Hashtable   iConfigHash;
    public Hashtable    iCellsByName;

    private double      iZPixRes;
    private int         iPlaneEnd;
    private int         iPlaneStart;

    private boolean     iShowAnnotations;
    private boolean     iShowAnnotationsSave;
    private boolean     iShowCentroids;
    //private Vector      iAnnotsShown;
    public Integer      iTrackPosition;
    public Integer      iTrackPositionSave;
    private boolean     iIgnoreValueChanged;

    private Image3D     iImage3D;
    private Image3D2    iImage3D2;
    private Image3D2Z    iImage3D2Z;
    private Object      iDispProps3D;
    private Object      iDispProps3D2;
    private Object      iDispProps3D2Z;
//    private EditImage   iEditImage;
//    private EditImage   iEditImage2;
    //public EditImage3   iEditImage3;
    public boolean		iEditTools;
    private CellMovementImage iCellMovementImage;
    private boolean     iCallSaveImage;
    private int		iUseZip;

    private CanonicalTree   iCanonicalTree;
    private PlayerControl   iPlayerControl;
    private EditTraverse    iEditTraverse;

    private boolean     iDebugTest;
    private int         iColor;

    public	NucRelinkDialog			iNucRelinkDialog;
    public	AddOneDialog			iAddOneDialog;


    /**
     * The only constructor defined for this class.
     * Instantiated in the main program.
     */
    private AceTree() {
        this(null);
        //System.out.println("AceTree empty constructor");
        //this("config.dat");
    }

    private AceTree(String configFileName) {
        super();
        iAceTree = this;
        iMainFrame = new JFrame(TITLE);
        iAceMenuBar = new AceMenuBar(this);
        iTempV = new Vector();
        iConfigFileName = configFileName;
        //System.out.println("AceTree constructor using config file: " + iConfigFileName);
        //NucUtils.setConfigFileName(iConfigFileName);
        iNucleiMgrHash = new Hashtable();
        iRootEstablished = false;
        iImageTime = 0;
        iImagePlane = 0;
        iTimeInc = 0;
        iPlaneInc = 0;
        iCurrentCell = null;
        iCurrentCellXloc = 0;
        iCurrentCellYloc = 0;
        setShowAnnotations(false);
        iShowCentroids = false;
        iInputCtrl = null;
        iRoot = new Cell(ROOTNAME);
        iTree = new JTree(iRoot);
        iTree.addMouseListener(new TreeMouseAdapter());
        iEditLog = new EditLog("EditLog");
        iDLog = new Log("Debug Log");

        setKeyboardActions();
        displayTree();
        iTrackPosition = ImageWindow.ANTERIOR;
        iTrackPositionSave = iTrackPosition;
        iDebugTest = false;
        iCanonicalTree = CanonicalTree.getCanonicalTree();
        iColor = 0;
        iTree.addMouseListener(new TreeMouseAdapter());
        if (iConfigFileName != null) {
            bringUpSeriesUI(iConfigFileName);
        }
        iTree.setCellRenderer(new AceTreeCellRenderer(this));
        setDefaultKeyboardActions();
    }

    public PlayerControl getPlayerControl() {
    	return iPlayerControl;
    }

    private void testRoot() {
        System.out.println("testRoot: " + iRoot.getLeafCount());
    }

    public synchronized static AceTree getAceTree(String configFileName) {
        if (iAceTree == null) {
            //System.out.println("AceTree.getAceTree making a new AceTree: " + configFileName);
            if (configFileName != null) {
                iAceTree = new AceTree(configFileName);

            } else iAceTree = new AceTree();

        }
        return iAceTree;
    }

    public synchronized static AceTree getAceTree(Object configFileName) {
        if (iAceTree == null) {
            //System.out.println("AceTree.getAceTree making a new AceTree: " + configFileName);
            if (configFileName != null) {
                String [] configs = ((String [])configFileName);
                if (configs.length == 1) {
                    iAceTree = new AceTree(configs[0]);
                } else {
                    iAceTree = new AceTree(configs[1]);
                    for (int i=2; i < configs.length; i++) {
                        iAceTree.setConfigFileName(configs[i]);
                        iAceTree.bringUpSeriesUI(configs[i]);
                    }
                    iAceTree.setConfigFileName(configs[0]);
                    iAceTree.bringUpSeriesUI(configs[0]);

                }
            } else iAceTree = new AceTree();

        }
        return iAceTree;
    }

    private class AceTreeCellRenderer extends DefaultTreeCellRenderer {
    	AceTree iAceTree;
    	public AceTreeCellRenderer(AceTree aceTree) {
    		super();
    		iAceTree = aceTree;
    	}

    	public Font getFont() {
    		Font f = super.getFont();
    		if (iAceTree != null) {
    			Cell x = iAceTree.getCellByName(getText());
    			if (x == null) {
    				f = f.deriveFont(Font.PLAIN);
    			} else if (x.getTime() < 1) {
    		    	f = f.deriveFont(Font.ITALIC);
    		    } else {
    			    f = f.deriveFont(Font.BOLD);
    		    }
    		}
    		return f;
    	}

    	public Icon getLeafIcon() {
    		return null;
    	}

    	public Icon getOpenIcon() {
    		return null;
    	}

    	public Icon getClosedIcon() {
    		return null;
    	}

    }


    public void bringUpSeriesUI(String configFileName) {
        System.out.println("\n\nbringUpSeriesUI: " + configFileName);
        System.gc();
        // check to see if the series is already in the hash
        String shortName = Config.getShortName(configFileName);
        NucleiMgr nucMgr = (NucleiMgr)iNucleiMgrHash.get(shortName);
        if (nucMgr == null) {
            // if not in hash then make sure there is such a file before proceeding
            try {
                FileInputStream fis = new FileInputStream(configFileName);
                fis.close();
            } catch(Exception fnfe) {
                new AceTreeHelp("/org/rhwlab/help/messages/ConfigError.html", 200, 200);
                return;
            }

            int k = bringUpSeriesData(configFileName);
            if (k != 0) return; //problem finding the zipNuclei
        }
        iNucleiMgr = (NucleiMgr)iNucleiMgrHash.get(shortName);
        if (iNucleiMgr == null) {
            System.out.println(HELPMSG + configFileName);
            System.exit(1);
        }
        iEditLog = iNucleiMgr.getEditLog();
        iNucleiMgr.sendStaticParametersToImageWindow();
        ImageWindow.setNucleiMgr(iNucleiMgr);
        //clearTree();
        setConfigFileName(configFileName);
        grabConfigStuff();
        iPlaneEnd = iNucleiMgr.getPlaneEnd();
        iPlaneStart = iNucleiMgr.getPlaneStart();

        //clearTree();
        //iTree.updateUI();


        buildTree(false);
        setShowAnnotations(true);

    }

    public void bringUpSeriesUI(Config config) {
    	String configFileName = config.iConfigFileName;
        System.out.println("\n\nbringUpSeriesUI: " + configFileName);
        System.gc();
        // check to see if the series is already in the hash
        String shortName = Config.getShortName(configFileName);
        NucleiMgr nucMgr = (NucleiMgr)iNucleiMgrHash.get(shortName);
        if (nucMgr == null) {
            // if not in hash then make sure there is such a file before proceeding

            int k = bringUpSeriesData(config);
            if (k != 0) return; //problem finding the zipNuclei
        }
        iNucleiMgr = (NucleiMgr)iNucleiMgrHash.get(shortName);
        if (iNucleiMgr == null) {
            System.out.println(HELPMSG + configFileName);
            System.exit(1);
        }
        iEditLog = iNucleiMgr.getEditLog();
        iNucleiMgr.sendStaticParametersToImageWindow();
        ImageWindow.setNucleiMgr(iNucleiMgr);
        //clearTree();
        setConfigFileName(configFileName);
        grabConfigStuff();
        iPlaneEnd = iNucleiMgr.getPlaneEnd();
        iPlaneStart = iNucleiMgr.getPlaneStart();

        //clearTree();
        //iTree.updateUI();
        buildTree(false);
        setShowAnnotations(true);

    }

    public int bringUpSeriesData(String configFileName) {
        System.out.println("bringUpSeriesData: " + configFileName);
        File fx = new File(configFileName);
        String ss = TITLE + ": " + fx.getName();
        iMainFrame.setTitle(ss);

        // this is the only place where we construct a NucleiMgr
        NucleiMgr nucMgr = new NucleiMgr(configFileName);
        if (!nucMgr.iGoodNucleiMgr) {
            return -1;
        }
        nucMgr.processNuclei(true, nucMgr.getConfig().iNamingMethod);
        String config = nucMgr.getConfig().getShortName();
        println("bringUpSeriesData, " + config);
        if (!iNucleiMgrHash.containsKey(config)) {
            iNucleiMgrHash.put(config, nucMgr);
            iAceMenuBar.addToRecent(config);
        }
        System.gc();
        return 0;
    }

    public int bringUpSeriesData(Config config) {
    	String configFileName = config.iConfigFileName;
        System.out.println("bringUpSeriesData: " + configFileName);
        File fx = new File(configFileName);
        String ss = TITLE + ": " + fx.getName();
        iMainFrame.setTitle(ss);

        // this is the only place where we construct a NucleiMgr
        NucleiMgr nucMgr = new NucleiMgr(config);
        if (!nucMgr.iGoodNucleiMgr) {
            return -1;
        }
        nucMgr.processNuclei(true, nucMgr.getConfig().iNamingMethod);
        String configName = nucMgr.getConfig().getShortName();
        if (!iNucleiMgrHash.containsKey(configName)) {
            iNucleiMgrHash.put(configName, nucMgr);
            iAceMenuBar.addToRecent(configName);
        }
        System.gc();
        return 0;
    }


    public void openSeveralConfigs(String configList) {
        String sr = null;
        boolean first = true;
        //iConfigFiles = new Vector();
        try {
            FileInputStream fis = new FileInputStream(configList);
            BufferedReader br = new BufferedReader(new InputStreamReader(fis));
            sr = br.readLine();
            while (sr != null && sr.length() > 2) {
                if (sr.indexOf("#") != 0) {
                    String [] sa = sr.split(" ");
                    sr = sa[0];
                    System.out.println("\n\n***config file: " + sr);
                    bringUpSeriesData(sr);
                    if (first) {
                        bringUpSeriesUI(sr);
                        first = false;
                    }
                }
                sr = br.readLine();
            }
            br.close();
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }


    }

    public void removeRecent(String item) {
        iNucleiMgrHash.remove(item);
        System.gc();
    }

    public void clearAll() {
        Enumeration e = iNucleiMgrHash.keys();
        while (e.hasMoreElements()) {
            NucleiMgr nm = (NucleiMgr)iNucleiMgrHash.get(e.nextElement());
            nm = null;
            System.gc();

        }
        iNucleiMgrHash = new Hashtable();
        System.gc();
    }

    private void grabConfigStuff() {
        Config c = iNucleiMgr.getConfig();
        iTifPrefix = c.iTifPrefix;
        iStartingIndex = c.iStartingIndex;
        iEndingIndex = c.iEndingIndex;
        iUseZip = c.iUseZip;

    }


    public void clearTree() {
        //new Throwable().printStackTrace();
        if (iAncesTree == null) {
            return;
        }
        Cell root = iRoot;
        if (root == null) return;
        int m = 0;
        int count = 0;
        Cell cc = null;
        while (iRoot.getChildCount() > 1) {
            Cell c = (Cell)iRoot.getFirstLeaf();
            while ((cc = (Cell)c.getNextLeaf()) != null) {
                c.removeFromParent();
                c = null;
                c = cc;
                count++;
            }
            m++;
        }
        println("clearTree: removed: " + count + CS + m);
        Runtime runtime = Runtime.getRuntime();
        runtime.gc();
        println("clearTree: memory: " + runtime.freeMemory() + CS + runtime.totalMemory() + CS + runtime.maxMemory());
        if (root != null) root.removeAllChildren();
        Hashtable x = iAncesTree.getCells();
        if (x != null) x.clear();
        iTree.updateUI();
    }

    private void clearTreeTest() {
        if (iAncesTree == null) {
            return;
        }
        Cell root = iAncesTree.getRoot();
        root = iRoot;
        System.out.println("\nAceTree.clearTree entered with root=" + root);
        int count = 0;
        /*
        println("clearTreeTest: leafcount: " + iRoot.getLeafCount());
        for (int i=0; i < 3000; i++) {
            int k = iRoot.getLeafCount();
            //println("clearTreeTest: " + i + CS + k);
            if (k == 1) break;
            Cell c = (Cell)iRoot.getFirstLeaf();
            //println("clearTreeTest: " + i + CS + k + CS + c.getName());
            c.removeFromParent();
            //c = null;
        }
        */
        int m = 0;
        Cell cc = null;
        while (iRoot.getChildCount() > 1) {
            //int k = iRoot.getChildCount();
            //println("clearTreeTest: childcount: " + k);
            Cell c = (Cell)iRoot.getFirstLeaf();
            while ((cc = (Cell)c.getNextLeaf()) != null) {
                c.removeFromParent();
                c = null;
                c = cc;
                count++;
            }
            m++;
        }
        println("clearTreeTest: removed: " + count + CS + m);


        if (root != null) root.removeAllChildren();
        Hashtable x = iAncesTree.getCells();
        if (x != null) x.clear();
        //iNucleiMgr.clearAllHashkeys(); //************ BOGUS
        iTree.updateUI();

    }

    void reviewNuclei() {
    	Vector nr = iNucleiMgr.getNucleiRecord();
    	for (int i=189; i < 195; i++) {
    		Vector nuclei = (Vector)nr.get(i);
    		for (int j=0; j < nuclei.size(); j++) {
    			Nucleus n = (Nucleus)nr.get(j);
    			println("reviewNuclei, " + i + CS + j );
    		}
    	}
    }

    public void buildTree(boolean doIdentity) {
        iShowAnnotationsSave = iShowAnnotations;
        setShowAnnotations(false);
        iShowCentroids = false;
        iShowC.setText(SHOWC);
        if (doIdentity) iNucleiMgr.processNuclei(doIdentity, iNamingMethod);

        if (iEditLog != null) {
            //iEditLog.append(new GregorianCalendar().getTime().toString());
            iEditLog.append("buildTree(" + doIdentity +
                ") start = " + iStartingIndex + " end = " + iEndingIndex
                + iEditLog.getTime());
        }

        grabConfigStuff();
        System.out.println("StartingIndex: " + iStartingIndex);
        System.out.println("EndingIndexX: " + iEndingIndex);
        Cell.setEndingIndexS(iEndingIndex);
        iAncesTree = iNucleiMgr.getAncesTree();
        iCellsByName = iAncesTree.getCellsByName();
        //Cell P = (Cell)iCellsByName.get("P");
        //int kkk = P.getChildCount();
        //println("AceTree.buildTree, 1, " + kkk + CS + P.getName());



        updateRoot(iAncesTree.getRootCells());

        iCellsByName = iAncesTree.getCellsByName();
        //Cell PP = (Cell)iCellsByName.get("P");
        //int kk = PP.getChildCount();
        //println("AceTree.buildTree, 2, " + kk + CS + PP.getName());
        //iAxis = getAxis();

        //System.out.println("buildTree: " + iRoot + CS + iRoot.getChildCount());
        int k = 0;
        Cell c = walkUpToAGoodCell();

        iAceMenuBar.setEditEnabled(true);
        iAceMenuBar.setEnabled(true);

        // 20050808 added in response to detected bug related to killCells
        iTree.updateUI();
        setTreeSelectionMode();
        setStartingCell(c, 1);
        iTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setOpenIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setClosedIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setLeafIcon(null);

        if (iEditTraverse != null) iEditTraverse.buildNotification();
        setShowAnnotations(iShowAnnotationsSave);

    }

    private Cell walkUpToAGoodCell() {
        Cell c = null;
        if (iRoot.getChildCount() <= 1) return iRoot;
        // assume the first child is P0
        // look for a real cell off of P0
        c = (Cell)iRoot.getChildAt(0);
        while (c.getTime() < 0 && c.getChildCount() > 0) {
            //System.out.println("buildTree: " + c + CS + c.getChildCount() + CS + k);
            //c = (Cell)iRoot.getChildAt(++k);
            c = (Cell)c.getChildAt(0);
            //System.out.println("buildTree: " + c);
        }
        // if you don't find one, go back to the root and look
        // for a Nuc or something
        if (c.getTime() < 0) {
            for (int i=1; i < iRoot.getChildCount(); i++) {
                c = (Cell)iRoot.getChildAt(i);
                if (c.getTime() > 0) break;
            }

        }
        return c;

    }

    /*
    private String getAxis() {
        String axis = "adl";
        Identity id = iNucleiMgr.getIdentity();
        if (id.getParameters().dvInit < 0) axis = "avr";
        println("\ngetAxis: " + axis);
        //return axis;
        return "sam";
    }
    */

    public void restoreTree(String shortName) {
        System.out.println("\n\nAceTree.restoreTree called: " + shortName);
        iMainFrame.setTitle(TITLE + ": " + shortName);
        //new Throwable().printStackTrace();
        NucleiMgr nucMgr = (NucleiMgr)iNucleiMgrHash.get(shortName);
        if (nucMgr == null) {
            System.out.println("SORRY: " + shortName + " is not hashed");
            return;
        }
        iNucleiMgr = nucMgr;
        iEditLog = iNucleiMgr.getEditLog();
        grabConfigStuff();
        iPlaneEnd = iNucleiMgr.getPlaneEnd();
        iPlaneStart = iNucleiMgr.getPlaneStart();
        //NucUtils.setNucleiMgr(iNucleiMgr);
        ImageWindow.setNucleiMgr(iNucleiMgr);
        ImageWindow.setStaticParameters(iNucleiMgr.getConfig().iZipTifFilePath,
                iNucleiMgr.getConfig().iTifPrefix, iNucleiMgr.getConfig().iUseZip);

        System.out.println("StartingIndex: " + iStartingIndex);
        System.out.println("EndingIndex: " + iEndingIndex);
        Cell.setEndingIndexS(iEndingIndex);
        iAncesTree = iNucleiMgr.getAncesTree();
        updateRoot(iAncesTree.getRootCells());
        iCellsByName = iAncesTree.getCellsByName();
        setShowAnnotations(false);
        iShow.setText(SHOW);
        iShowCentroids = false;
        iShowC.setText(SHOWC);
        Cell.setEndingIndexS(iEndingIndex); // what does this do?

        Cell c = walkUpToAGoodCell();
        iTree.updateUI();
        setTreeSelectionMode();
        setStartingCell(c, 1);
        iTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setOpenIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setClosedIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setLeafIcon(null);
        //iAxis = getAxis();
    }


    private void updateRoot(Vector rootCells) {
        Cell PP = (Cell)iCellsByName.get("P");
        int kk = PP.getChildCount();
        //println("AceTree.updateRoot, 1, " + kk + CS + PP.getName());
        //System.out.println("\n#######updateRoot in: " + iRoot.showStuff());
        iRoot.removeAllChildren();

        //PP = (Cell)iCellsByName.get("P");
        //kk = PP.getChildCount();
       // println("AceTree.updateRoot, 2, " + kk + CS + PP.getName());


        // struggled with what must be a bug in DefaultMutableTreeNode
        // until I broke out the collecting of children
        // from the adding of them to a different root
        Vector v = new Vector();
        Enumeration e = rootCells.elements();
        while (e.hasMoreElements()) {
            Cell c = (Cell)e.nextElement();
            v.add(c);
        }
        //PP = (Cell)iCellsByName.get("P");
        //kk = PP.getChildCount();
        //println("AceTree.updateRoot, 3, " + kk + CS + PP.getName() + CS + iCellsByName.size());

        for (int i=0; i < v.size(); i++) {
            Cell cc = (Cell)v.elementAt(i);
            //println("AceTree.updateRoot, " + i + CS + cc.getName() + CS + ((Cell)cc.getParent()).getName());
            cc.removeFromParent();
            iRoot.add(cc);
        }
        //iCellsByName.put("P", PP);
        iCellsByName.remove("P");
        iCellsByName.put("P", iRoot);


        //PP = (Cell)iCellsByName.get("P");
        //kk = PP.getChildCount();
        //println("AceTree.updateRoot, 4, " + kk + CS + PP.getName() + CS + iCellsByName.size());

		iRoot.setEndTime(1);
		//System.out.println("\n#######updateRoot out: " + iRoot.showStuff());
		//Cell xx = (Cell)iAncesTree.getCellsByName().put(iRoot.getName(), iRoot);
		//System.out.println("\n#######updateRoot out2: " + xx.showStuff());
    }

    public void setStartingCell(Cell c, int time) {
        // seem to need to exercise iAncesTree to start things off well
        System.out.println("setStartingCell, cell, time: " + c + CS + time);


	    //new Throwable().printStackTrace();
        if (c != iRoot) {
            if (c == null) c = (Cell)iRoot.getChildAt(0);
            while (c.getChildCount() > 0 && c.getTime() < 1) {
                //println("setStartingCell while loop: " + c + CS + c.getTime());
                c = (Cell)c.getChildAt(0);
            }


            //c.showParameters();
            time = Math.max(time, c.getTime());
            time = Math.min(time, c.getEndTime());
            iImageTime = time;
            iTimeInc = 0;
            getTimeAndPlane(c);
            getCurrentCellParameters();
            showTreeCell(iCurrentCell);
        } else {
            iImageTime = 1;
            iTimeInc = 0;
            iImagePlane = 15;
            iPlaneInc = 0;
        }
        handleCellSelectionChange(c, time - iImageTime); // this will bring up an image
        if (!c.getName().equals("P") && iRoot.getChildCount() > 0) {
            //setShowAnnotations(true);
            iShowCentroids = true;
            iShowC.setText(HIDEC);
            addMainAnnotation();
        }
        iAceMenuBar.setClearEnabled(true);
        //System.out.println("setStartingCell -iImgWin: " + iImgWin);
        iImgWin.refreshDisplay(null);
    }

    public void run() {
        println("AceTree.run: entered");
        try {
            Thread.sleep(5);
        } catch(Exception e) {

        }
        nextTime();
        updateDisplay();
        try {
            Thread.sleep(5);
        } catch(Exception e) {

        }
        prevTime();
        updateDisplay();
    }

    public void expandTree() {
        Cell c = (Cell)iRoot.getFirstLeaf();
        while (c != null) {
            showTreeCell(c);
            //System.out.println(c);
            //TreeNode [] tna = c.getPath();
            //TreePath tp = new TreePath(tna);
            //iTree.makeVisible(tp);
            c = (Cell)c.getNextLeaf();
        }
    }

    private void setTreeSelectionMode() {
        iTree.getSelectionModel().setSelectionMode
        (TreeSelectionModel.SINGLE_TREE_SELECTION);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setOpenIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setClosedIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setLeafIcon(null);
    }

    private void displayTree() {
        iTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setOpenIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setClosedIcon(null);
        ((DefaultTreeCellRenderer)(iTree.getCellRenderer())).setLeafIcon(null);
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        JPanel treev = new JPanel();
        treev.setLayout(new BorderLayout());
        JScrollPane treeView = new JScrollPane(iTree);
        treeView.setPreferredSize(new Dimension(WIDTH,HEIGHT200));
        treeView.setMaximumSize(new Dimension(Integer.MAX_VALUE,Integer.MAX_VALUE));
        treev.add(treeView);
        add(treev);

        JPanel textv = new JPanel();
        textv.setLayout(new BorderLayout());
        textv.setMaximumSize(new Dimension(Integer.MAX_VALUE,Integer.MAX_VALUE));
        textv.setPreferredSize(new Dimension(WIDTH,HEIGHT100));
        iText = new JTextPane();
        iText.setEditable(false);
        JScrollPane textView = new JScrollPane(iText);
        textView.setPreferredSize(new Dimension(WIDTH,HEIGHT100));
        textv.add(textView);
        add(textv);
        iInputCtrl = new InputCtrl(this);
        add(iInputCtrl);
        iPlayerControl = new PlayerControl(this);
        add(iPlayerControl);
        JPanel pad = createPad();
        add(pad);
        //iTree.addTreeSelectionListener(this);

        JPanel textv2 = new JPanel();
        textv2.setLayout(new BorderLayout());
        iText2 = new JTextPane();
        iText2.setPreferredSize(new Dimension(WIDTH/2, HEIGHT30));
        iText2.setEditable(false);
        textv2.add(iText2, BorderLayout.WEST);
        iText3 = new JTextPane();
        iText3.setPreferredSize(new Dimension(WIDTH/2, HEIGHT30));
        iText3.setEditable(false);
        textv2.add(iText3, BorderLayout.EAST);


        add(textv2);
    }

    private JPanel createPad() {
        JPanel p = new JPanel();
        iNext = new JButton(NEXTT);
        iPrev = new JButton(PREV);
        iUp = new JButton(UP);
        iDown = new JButton(DOWN);
        iHome = new JButton(HOME);
        iDefault = iHome;
    	iMainFrame.getRootPane().setDefaultButton(iDefault);

        //iShowC = new JButton(SHOWC);
        iNext.addActionListener(this);
        iPrev.addActionListener(this);
        iUp.addActionListener(this);
        iDown.addActionListener(this);
        iHome.addActionListener(this);
        p.setLayout(new GridLayout(4,3));

        iShow = new JButton(SHOW);
        iShowC = new JButton(SHOWC);
        iClear = new JButton(CLEAR);
        iShow.addActionListener(this);
        iShowC.addActionListener(this);
        iClear.addActionListener(this);
        //JButton x1 = new JButton("");
        //JButton x3 = new JButton("");

        iCopy = new JButton(COPY);
        iCopy.addActionListener(this);
        iTrack = new JButton(TRACK);
        iTrack.addActionListener(this);
        iSister = new JButton(SISTER);
        iSister.addActionListener(this);
        iColorToggle = new JButton(COLORTOGGLE);
        iColorToggle.addActionListener(this);
        //iEdit = new JButton(EDIT);
        //iEdit.addActionListener(this);

        p.setPreferredSize(new Dimension(WIDTH, HEIGHT100));
        p.add(iShow);
        p.add(iUp);
        p.add(iClear);
        p.add(iPrev);
        p.add(iHome);
        p.add(iNext);
        p.add(iShowC);
        p.add(iDown);
        p.add(iCopy);
        p.add(iTrack);
        p.add(iSister);
        p.add(iColorToggle);
        //p.add(iEdit);
        return p;
    }


    private void setSpecialKeyboardActions() {
    	KeyStroke key = null;
    	String xxx = null;
    	final AceTree aceTree = this;

    	xxx = "ctrl_left";
    	Action ctrl_left = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_left");
    			if (iAddOneDialog != null) {
    				ActionEvent ae = new ActionEvent(aceTree, 1, "LEFT");
    				iAddOneDialog.actionPerformed(ae);
    			}
    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_left );

    	xxx = "ctrl_left_a";
    	Action ctrl_left_a = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_left");
    			if (iAddOneDialog != null) {
    				ActionEvent ae = new ActionEvent(aceTree, 1, "LEFT");
    				iAddOneDialog.actionPerformed(ae);
    			}
    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_left );

    	xxx = "ctrl_right";
        Action ctrl_right = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_right");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "RIGHT");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_right );

    	xxx = "ctrl_right_d";
        Action ctrl_right_d = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_right");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "RIGHT");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_right );

        xxx = "ctrl_up";
    	Action ctrl_up = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_up");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "UP");
				//ActionEvent ae = new ActionEvent(aceTree, 1, "DOWN");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_up );

        xxx = "ctrl_up_w";
    	Action ctrl_up_w = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_up");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "UP");
				//ActionEvent ae = new ActionEvent(aceTree, 1, "DOWN");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_W, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_up );

    	xxx = "ctrl_down";
        Action ctrl_down = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_down");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "DOWN");
				//ActionEvent ae = new ActionEvent(aceTree, 1, "UP");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_down );

    	xxx = "ctrl_down_s";
        Action ctrl_down_s = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, ctrl_down");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "DOWN");
				//ActionEvent ae = new ActionEvent(aceTree, 1, "UP");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, ctrl_down );

    	xxx = "shift_left";
    	Action shift_left = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_left");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "BIG");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_LEFT, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_left );

    	xxx = "shift_left_a";
    	Action shift_left_a = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_left");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "BIG");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_A, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_left );

    	xxx = "shift_right";
        Action shift_right = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_right");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "SMALL");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_RIGHT, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_right );

    	xxx = "shift_right_d";
        Action shift_right_d = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_right");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, "SMALL");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_D, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_right );

        xxx = "shift_up";
    	Action shift_up = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_up");
    			if (iAddOneDialog == null) return;
				//ActionEvent ae = new ActionEvent(aceTree, 1, "INCZ");
				ActionEvent ae = new ActionEvent(aceTree, 1, "DECZ");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_UP, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_up );

        xxx = "shift_up_w";
    	Action shift_up_w = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_up");
    			if (iAddOneDialog == null) return;
				//ActionEvent ae = new ActionEvent(aceTree, 1, "INCZ");
				ActionEvent ae = new ActionEvent(aceTree, 1, "DECZ");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_W, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_up );

    	xxx = "shift_down";
        Action shift_down = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_down");
    			if (iAddOneDialog == null) return;
				//ActionEvent ae = new ActionEvent(aceTree, 1, "DECZ");
				ActionEvent ae = new ActionEvent(aceTree, 1, "INCZ");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_down );

    	xxx = "shift_down_s";
        Action shift_down_s = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, shift_down");
    			if (iAddOneDialog == null) return;
				//ActionEvent ae = new ActionEvent(aceTree, 1, "DECZ");
				ActionEvent ae = new ActionEvent(aceTree, 1, "INCZ");
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.SHIFT_MASK, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, shift_down );

        // this one is a delete cell special
    	xxx = "DELETE";
        Action DELETE = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, DELETE");
    			if (iAddOneDialog == null) return;
    			killCell(0);
    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_DELETE, 0, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, DELETE );

        xxx = "F5";
        Action F5 = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, F5");
    			if (iAddOneDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, AddOneDialog.REBUILDANDRENAME);
				iAddOneDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_F5, 0, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, F5 );

        // these go to the NucRelinkDialog

        xxx = "F1";
        Action F1 = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, F1");
    			if (iNucRelinkDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, NucRelinkDialog.SETEARLYCELL);
				iNucRelinkDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_F1, 0, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, F1 );

    	xxx = "F2";
        Action F2 = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, F2");
    			if (iNucRelinkDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, NucRelinkDialog.SETLATECELL);
				iNucRelinkDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_F2, 0, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, F2 );

    	xxx = "F3";
        Action F3 = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, F3");
    			if (iNucRelinkDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, NucRelinkDialog.APPLYONLY);
				iNucRelinkDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_F3, 0, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, F3 );

    	xxx = "F4";
        Action F4 = new AbstractAction() {
    		public void actionPerformed(ActionEvent e) {
    			//println("AceTree.setSpecialKeyBoardActions, F4");
    			if (iNucRelinkDialog == null) return;
				ActionEvent ae = new ActionEvent(aceTree, 1, NucRelinkDialog.APPLYANDREBUILD);
				iNucRelinkDialog.actionPerformed(ae);

    		}
    	};
        key = KeyStroke.getKeyStroke(KeyEvent.VK_F4, 0, false);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(key, xxx);
        getActionMap().put(xxx, F4 );







    }

    private void setKeyboardActions() {

    	setSpecialKeyboardActions();

    	String actionKey = "";
    	KeyStroke stroke = null;
    	InputMap inputMap = null;
    	ActionMap actionMap = null;


        String s = "PAGE_UP";
        Action PageUp = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                if (iShow.getText().equals(HIDE)) {
                    setShowAnnotations(false);
                } else {
                    setShowAnnotations(true);
                }
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(KeyStroke.getKeyStroke("PAGE_UP"), "PAGE_UP");
        getActionMap().put(s, PageUp );

        s = "PAGE_DOWN";
        Action PageDn = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                if (iShowC.getText().equals(HIDEC)) {
                    iShowCentroids = false;
                    iShowC.setText(SHOWC);
                } else {
                    iShowCentroids = true;
                    iShowC.setText(HIDEC);
                }
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).
            put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, PageDn );

        s = "END";
        Action end = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                copyImage();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).
            put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, end );

        s = "HOME";
        Action home = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                iTimeInc = iPlaneInc = 0;
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).
            put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, home );

        s = "UP";
        Action up = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                //System.out.println("up key pressed");
                incPlane(-1);
                iTrackPosition = ImageWindow.NONE;
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).
            put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, up );

        actionKey = "w_up";
        stroke = KeyStroke.getKeyStroke("typed w");
        inputMap = this.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        inputMap.put(stroke, actionKey);
        actionMap = this.getActionMap();
        actionMap.put(actionKey, up);

        s = "DOWN";
        Action down = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                //System.out.println("down key pressed");
                incPlane(1);
                iTrackPosition = ImageWindow.NONE;
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).
            put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, down );

        actionKey = "s_down";
        stroke = KeyStroke.getKeyStroke("typed s");
        inputMap = this.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        inputMap.put(stroke, actionKey);
        actionMap = this.getActionMap();
        actionMap.put(actionKey, down);


        s = "LEFT";
        Action left = new AbstractAction("LEFT") {
            public void actionPerformed(ActionEvent e) {
            	AWTEvent awt = (AWTEvent)e;
                //System.out.println("AceTree--left key pressed, " + e);
                //System.out.println("AceTree--left key pressed, awt, " + awt);
                prevTime();
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, left );

        actionKey = "a_left";
        stroke = KeyStroke.getKeyStroke("typed a");
        inputMap = this.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        inputMap.put(stroke, actionKey);
        actionMap = this.getActionMap();
        actionMap.put(actionKey, left);



        s = "RIGHT";
        Action right = new AbstractAction(s) {
            public void actionPerformed(ActionEvent e) {
                //System.out.println("right key pressed");
                if (nextTime()) updateDisplay();
            }
        };

        //AceTreeActions right = new AceTreeActions("RIGHT", 12345);
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, right );

        actionKey = "d_right";
        stroke = KeyStroke.getKeyStroke("typed d");
        inputMap = this.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        inputMap.put(stroke, actionKey);
        actionMap = this.getActionMap();
        actionMap.put(actionKey, right);

        s = "ENTER";
        Action get = new AbstractAction(s) {
            public void actionPerformed(ActionEvent e) {
                iInputCtrl.getIt();
                updateDisplay();
            }
        };
        getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).
            put(KeyStroke.getKeyStroke(s), s);
        getActionMap().put(s, get );

    }

    class AceTreeActions extends AbstractAction {
    	int iID;

    	public AceTreeActions(String name, int id) {
    		super(name);
    		iID = id;

    	}

		public void actionPerformed(ActionEvent e) {
			println("AceTreeActions.actionPerformed, " + e);

		}

    }

    private void setDefaultKeyboardActions() {
        String s = "F2";
        Action home = new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
            	println("setDefaultKeyboardActions, ");
            	//iAceTree.requestFocus();
            	Component compFocusOwner =
                    KeyboardFocusManager.getCurrentKeyboardFocusManager().getFocusOwner();
            	Window windowFocusOwner =
                    KeyboardFocusManager.getCurrentKeyboardFocusManager().getFocusedWindow();
            	if (compFocusOwner instanceof JButton) {
            		println("its a button");
            		//((JButton)compFocusOwner).doClick();
            	}
            	println("setKeyboardActions, " + compFocusOwner);
            }
        };
        iDefault.getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).
            put(KeyStroke.getKeyStroke(s), s);
        iDefault.getActionMap().put(s, home );

    }


    /////////////////////////////////////////////////////////////

    private void handleCellSelectionChange(Cell c, int timeInc) {
        //System.out.println("handleCellSelectionChange: " + c + CS + timeInc);
        if (c == null) return;
        //iAnnotsShown.clear();
        getTimeAndPlane(c);
        iTimeInc = timeInc;
        iPlaneInc = 0;
        //println("handleCellSelectionChange:2 " + iImageTime + CS + iImagePlane);
        if (iImageTime < 1 || iImagePlane < 1) return;
        updateDisplay();
    }

    private void getTimeAndPlane(Cell c) {
        if (c == null) return;
        if (c == iRoot) {
            iImageTime = 1;
            iImagePlane = 15;
        } else {
            iImageTime = c.getTime();
            iImagePlane = (int)((double)c.getPlane() + HALFROUND);
        }
        iTimeInc = 0;
        iPlaneInc = 0;
        iCurrentCell = c;
    }


    public void updateDisplay() {
    	//println("updateDisplay:1 " + System.currentTimeMillis());
        if (iDebugTest) println("updateDisplay:1 " + System.currentTimeMillis());
        if ((iImageTime + iTimeInc) < iStartingIndex) return;
        if ((iImagePlane + iPlaneInc) <= 0) return;
        getCurrentCellParameters();
        handleImage();
        if (iCallSaveImage) {
            iCallSaveImage = false;
            iImgWin.saveImageIfEnabled();
        }
        String s = makeDisplayText();
        iText.setText(s);
        if (iDebugTest) {
            debugTest(false);
            println("updateDisplay:2 " + System.currentTimeMillis());
        }
    }

    public void handleImage() {
        String cfile = makeImageName();
        ImagePlus ip = null;
        if (cfile == null) {
            IJ.error("no image available");
            iImgWin.makeImage(null);
            return;
        } else {
            if (iImgWin != null) {
                try {

                    ip = iImgWin.refreshDisplay(iTifPrefix + cfile);
                } catch(Exception e) {
                    System.out.println("handleImage -- no image available: " + iTifPrefix + cfile);
                    System.out.println(e);
                    e.printStackTrace();
                    iPlayerControl.pause();
                }
            } else {
                System.out.println("\nhandleImage making new one: " + ip + CS + iTifPrefix + CS + cfile);
                ip = ImageWindow.makeImage(iTifPrefix + cfile);
                iImgWin = new ImageWindow(iTifPrefix + cfile, ip);
                iImgWin.setAceTree(this);
                //iImgWin.refreshDisplay(iTifPrefix + makeImageName(iCurrentCell);
                iImgWinSet = true;
            }
        }

        /*
        if (iEditImage != null) {
            try {
                ip = iEditImage.refreshDisplay(iTifPrefix + cfile);
            } catch(Exception e) {
                System.out.println("handleImage -- no image available: " + iTifPrefix + cfile);
                System.out.println(e);
                e.printStackTrace();
                //iPlayerControl.pause();
            }
        }
        */

        /*
        if (iEditImage3 != null) {
            try {
                ip = iEditImage3.refreshDisplay(iTifPrefix + cfile);
            } catch(Exception e) {
                System.out.println("handleImage -- no image available: " + iTifPrefix + cfile);
                System.out.println(e);
                e.printStackTrace();
                //iPlayerControl.pause();
            }
        }
        //iImgWin.requestFocus();
         * */

    }

    public void addMainAnnotation() {
        //System.out.println("addMainAnnotation: " + iCurrentCellXloc + CS + iCurrentCellYloc);
        if (iCurrentCellXloc <= 0) return;
        iImgWin.addAnnotation(iCurrentCellXloc, iCurrentCellYloc, true);
        //if (iEditImage != null) iEditImage.addAnnotation(iCurrentCellXloc, iCurrentCellYloc, true);
        //if (iEditImage3 != null) iEditImage3.addAnnotation(iCurrentCellXloc, iCurrentCellYloc, true);
    }

    private String makeImageName() {
        // typical name: t001-p15.tif
        // to be augmented later to something like: images/050405-t001-p15.tif
        // which specifies a path and prefix for the set
        StringBuffer name = new StringBuffer("t");
        name.append(EUtils.makePaddedInt(iImageTime + iTimeInc));
        name.append("-p");
        String p = EUtils.makePaddedInt(iImagePlane + iPlaneInc, 3);
        if (p.startsWith("0")) p = p.substring(1, p.length());
        name.append(p);
        switch(iUseZip) {
        case 0:
        case 1:
        case 3:
            name.append(".tif");
            break;
        default:
            name.append(".zip");
        }
        return(name.toString());
    }

    public String makeImageName(int time, int plane) {
        // typical name: t001-p15.tif
        // to be augmented later to something like: images/050405-t001-p15.tif
        // which specifies a path and prefix for the set
        StringBuffer name = new StringBuffer("t");
        name.append(EUtils.makePaddedInt(time));
        name.append("-p");
        String p = EUtils.makePaddedInt(plane, 3);
        if (p.startsWith("0")) p = p.substring(1, p.length());
        name.append(p);

        switch(iUseZip) {
        case 0:
        case 1:
            name.append(".tif");
            break;
        default:
            name.append(".zip");
        }
        return(name.toString());
    }


    private int trackCellPlane() {
        if (iTrackPosition != ImageWindow.NONE) {
            iPlaneInc = 0;
            return (int)(iCurrentCellZloc + 0.5);
        }
        else {
            return iImagePlane;
        }
    }


    private void getCurrentCellParameters() {
    	//System.out.println("getCurrentCellParameters: " + iImageTime + CS + iTimeInc);
        if (iCurrentCell == null) return;
        int time = iImageTime + iTimeInc;
        if (time == 0) {
            time = 1;
            iImageTime = 1;
        }
        //Vector nuclei = iNucleiMgr.getNucleiRecord()[time - 1];
        Nucleus n = null;
        try {
            Vector nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(time - 1);
            n = NucUtils.getCurrentCellNucleus(nuclei, iCurrentCell);
        } catch(Exception e) {
            System.out.println("AceTree.getCurrentCellParameters error at time=" + time);
        }
        //System.out.println("getCurrentCellParameters: " + time + CS + iCurrentCell + CS + n);
        iCurrentCellXloc = -1;
        iCurrentCellYloc = -1;
        iCurrentCellZloc = -1;
        iCurrentCellPresent = false;
        if (n != null) {
            iCurrentCellXloc = n.x;
            iCurrentCellYloc = n.y;
            iCurrentCellZloc = n.z;
            iImagePlane = trackCellPlane();
            iCurrentCellPresent = true;
        }
    }

    private String getRedDataFromCell(int time) {
        Vector vcd = iCurrentCell.getCellData();
        int item = time - iCurrentCell.getTime();
        CellData cd = (CellData)vcd.elementAt(item);
        String s = ", rweight: " + cd.iNucleus.rweight;
        return s;
    }

    private String makeDisplayText() {
        int time = iImageTime + iTimeInc;
        Vector nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(time - 1);

        String name = "";
        if (iCurrentCell != null) name = iCurrentCell.getName();
        StringBuffer sb2 = new StringBuffer();
        //System.out.println("makeDisplayText: " + name);
        if (iCurrentCell == null) iCurrentCellPresent = false;
        if (iCurrentCellPresent) {
            sb2.append(name + " is one of ");
            sb2.append(NucUtils.countLiveCells(nuclei) + " cells at time " + (iImageTime + iTimeInc));
            Nucleus n = NucUtils.getCurrentCellNucleus(nuclei, iCurrentCell);
	        if (n != null) {
            sb2.append("\nlocation: " + iCurrentCellXloc + ", " + iCurrentCellYloc + ", " + n.z);
	        //sb2.append(CS + iAxis);
            double d = iNucleiMgr.nucDiameter(n,
                    (double)(iImagePlane + iPlaneInc));
            String sd = (new DecimalFormat("###.#")).format(d);
            sb2.append("\nsize: " + n.size + " displayed diameter: " + sd);
            sb2.append("\ncurrent index: " + n.index);
            //sb2.append(getRedDataFromCell(time));
            sb2.append(" weightg/r: " + n.weight);
            sb2.append(", " + n.rweight);
            sb2.append("\nstart=" + iCurrentCell.getTime());
            sb2.append(", end=" + iCurrentCell.getEnd());
            sb2.append(", fate=" + iCurrentCell.getFate());
            String track;
            switch(iTrackPosition.intValue()) {
                case 1:
                    track = "\ntrack anterior";
                    break;
                case 2:
                    track = "\ntrack posterior";
                    break;
                default:
                    track = "\nnot tracking";
                    break;

            }
            //if (iTrackPosition != ImageWindow.NONE) track = "\ntracking";
            sb2.append(track);
            }
        } else sb2.append(name + " not present");
        //System.out.println("makeDisplayText: " + iCurrentCell.getRedDataString());
        return sb2.toString();

    }

    /**
     * Called from AceMenuBar on file open action
     * @param name the name of the config file we opened
     */
    public void setConfigFileName(String name) {
        //System.out.println("AceTree.setConfigFileName: " + name);
        //new Throwable().printStackTrace();
        iConfigFileName = name;
        //NucUtils.setConfigFileName(name);
    }

    /**
     * Called from AceMenuBar on quickopen action
     * @param
     */
    public void quickOpen() {
        println("AceTree.quickOpen");
        new QuickOpen();

    }


    public String getConfigFileName() {
        //System.out.println("AceTree.getConfigFileName: " + iConfigFileName);
        return iConfigFileName;
    }
////////////////////////////////////////////////////////////////////
////////////  image handling ///////////////////////////////////////
////////////////////////////////////////////////////////////////////


    private void copyImage() {
        ImagePlus ip = iImgWin.getImagePlus();
        String s = iTifPrefix + makeImageName() + Math.random();
        ImageWindow iw = new ImageWindow(s, ip);
        iw.setAceTree(this);
        iw.refreshDisplay(s);
        iw.setLocation(iw.getX() + XINC, iw.getY() + YINC);
    }

/*
    public void editImage() {
        ImagePlus ip = iImgWin.getImagePlus();
        String s1 = iTifPrefix + makeImageName();
        String s = s1 + Math.random();
        //System.out.println("editImage: " + s);
        iEditImage = new EditImage(s, ip);
        iEditImage.setLocation(iEditImage.getX() + XINC, iEditImage.getY() + YINC);
        if (iEditImage != null) {
            try {
                ip = iEditImage.refreshDisplay(s1);
            } catch(Exception e) {
                System.out.println("editImage -- no image available: " + s1);
                System.out.println(e);
                e.printStackTrace();
                //iPlayerControl.pause();
            }
        }
    }
*/
    public void editImage3() {
    	/*
        ImagePlus ip = iImgWin.getImagePlus();
        String s1 = iTifPrefix + makeImageName();
        String s = s1 + Math.random();
        //System.out.println("editImage: " + s);
        iEditImage3 = new EditImage3(s, ip);
        iEditImage3.setLocation(iImgWin.getX() + XINC, iImgWin.getY() + YINC);
        if (iEditImage3 != null) {
            try {
                ip = iEditImage3.refreshDisplay(s1);
            } catch(Exception e) {
                System.out.println("editImage -- no image available: " + s1);
                System.out.println(e);
                e.printStackTrace();
                //iPlayerControl.pause();
            }
        }
        */
    }

    public void editTools() {
    	//editImage3();
    	if (iAddOneDialog != null) return;
    	iEditTools = true;
    	iAddOneDialog = new AddOneDialog(this, iImgWin, false, iCurrentCell, iImageTime);
    	//iAddOneDialog = new AddOneDialog(this, iEditImage3, false, iCurrentCell, iImageTime);
    	iImgWin.iDialog = iAddOneDialog;
    	iImgWin.setDialogsEnabled(true);
    	//iEditImage3.iDialog = iAddOneDialog;
    	relinkNucleus(); //brings up the NucRelinkDialog
    	//println("AceTree.editTools, ");
    }

    //public EditImage3 getEditImage3() {
	//return iEditImage3;
    //}

    public void cellMovementImage() {
        ImagePlus ip = iImgWin.getImagePlus();
        String s1 = iTifPrefix + makeImageName();
        String s = s1;
        //String s1 = iTifPrefix + makeImageName();
        //String s = s1 + Math.random();
        //System.out.println("cellMovementImage: " + s);
        iCellMovementImage = new CellMovementImage(s, ip);
        /*
        if (iEditImage != null) {
            try {
                ip = iCellMovementImage.refreshDisplay(s1);
            } catch(Exception e) {
                System.out.println("cellMovememtImage -- no image available: " + s1);
                System.out.println(e);
                e.printStackTrace();
                //iPlayerControl.pause();
            }
        }
        */

    }


    public void setEditImageNull(int which) {
        //println("setEditImageNull:");
        switch(which) {
            case 1:
                //iEditImage = null;
                break;
            case 3:
                //iEditImage3 = null;
                break;
            case 4:
                iCellMovementImage = null;
                break;

            default:
        }
    }

////////////////////////////////////////////////////////////////////
////////////image handling end ///////////////////////////////////
////////////////////////////////////////////////////////////////////

// introduced to permit right click on tree to select end time of cell
    private class TreeMouseAdapter extends MouseInputAdapter {
        public void mouseClicked(MouseEvent e) {
            int button = e.getButton();
            //System.out.println("TreeMouseAdapter.mouseClicked: " + button);

            Cell c = null;
            if (button == 2) return;
            else {
                int selRow = iTree.getRowForLocation(e.getX(), e.getY());
                TreePath selPath = iTree.getPathForLocation(e.getX(), e.getY());
                if (selPath == null) return;
                c = (Cell)selPath.getLastPathComponent();
                //if (c.getName().startsWith("AB")) {
                //	println("TreeMouseAdapter.mouseClicked: debug");
                //}

            }

            if (button == 1) {
                //int selRow = iTree.getRowForLocation(e.getX(), e.getY());
                //TreePath selPath = iTree.getPathForLocation(e.getX(), e.getY());
                //if (selPath == null) return;
                //Cell c = (Cell)selPath.getLastPathComponent();
                //Cell c = (Cell) iTree.getLastSelectedPathComponent();
                showTreeCell(c);
                if (c.getTime() < 0) {
                    if (iCurrentCell != null) c = iCurrentCell;
                    else return;

                }
                int time = c.getTime();
                //println("TreeMouseAdapter.mouseClicked: " + selPath + CS + c + CS + time);
                setCurrentCell(c, time, LEFTCLICKONTREE);

            }

            else if (button == 3) {
                //iIgnoreValueChanged = true;
                //int selRow = iTree.getRowForLocation(e.getX(), e.getY());
                //TreePath selPath = iTree.getPathForLocation(e.getX(), e.getY());
                //if (selPath == null) return;
                //Cell c = (Cell)selPath.getLastPathComponent();
                //Cell c = (Cell) iTree.getLastSelectedPathComponent();
                showTreeCell(c);
                if (c.getTime() < 0) {
                    if (iCurrentCell != null) c = iCurrentCell;
                    else return;

                }
                int time = c.getEnd();
                setCurrentCell(c, time, RIGHTCLICKONTREE);
            }
        }
    }


    /** Required by TreeSelectionListener interface.
     * */
    public void valueChanged(TreeSelectionEvent e) {
        // events are ignored if they are created programmatically
        // as part of the tracking code
        //System.out.println("valueChanged: " + iIgnoreValueChanged);
        /*
        if (iIgnoreValueChanged) {
            iIgnoreValueChanged = false;
            return;
        }
        int now = iImageTime + iTimeInc;
        Cell c = (Cell) iTree.getLastSelectedPathComponent();
        setCurrentCell(c, now, LEFTCLICKONTREE);
        */
    }


    public void imageUp() {
    	println("AceTree.imageUp, ");
        incPlane(-1);
        iTrackPosition = ImageWindow.NONE;
	iCallSaveImage = true;
        updateDisplay();
    }

    public void imageDown() {
        incPlane(1);
        iTrackPosition = ImageWindow.NONE;
	iCallSaveImage = true;
        updateDisplay();
    }


    public void actionPerformed(ActionEvent e) {
    	//println("AceTree.actionPerformed, " + e);
        boolean doUpdate = true;
        if (!iImgWinSet) return;
        iImgWin.setSpecialEffect(null);
        String cmd = e.getActionCommand();
        int inc = 0;
        if (e.getActionCommand().equals(NEXTT)) {
            doUpdate = nextTime(); //incTime(1);
        }
        else if (cmd.equals("F2")) println("AceTree.actionPerformed, F2");
        else if (e.getActionCommand().equals(PREV)) prevTime(); //incTime(-1);
        else if (e.getActionCommand().equals(UP)) {
            imageUp();
            return;
        }
        else if (e.getActionCommand().equals(DOWN)) {
            imageDown();
            return;
        }
        else if (e.getActionCommand().equals(HOME)) {
            getTimeAndPlane(iCurrentCell);
            if (iCurrentCell.isAnterior()) iTrackPosition = ImageWindow.ANTERIOR;
            else iTrackPosition = ImageWindow.POSTERIOR;
            //setTrack();
        }
        else if (e.getActionCommand().equals(SHOW)) {
            setShowAnnotations(true);
        }
        else if (e.getActionCommand().equals(HIDE)) {
            setShowAnnotations(false);
        }
        else if (e.getActionCommand().equals(SHOWC)) {
            iShowCentroids = true;
            iShowC.setText(HIDEC);
        }
        else if (e.getActionCommand().equals(HIDEC)) {
            iShowCentroids = false;
            iShowC.setText(SHOWC);
        }
        else if (e.getActionCommand().equals(CLEAR)) {
            setShowAnnotations(false);
            iImgWin.clearAnnotations();
            //if (iEditImage3 != null) iEditImage3.clearAnnotations();
        }
        else if (e.getActionCommand().equals(COPY)) {
            //copyImage();
            //John requested that we trashcan this
            //for (int i=0; i < 10; i++) {
            //    clearTree();
            //    buildTree(true);
            //}
            //debugTest(true);
            //nextTime();
        }

        else if (e.getActionCommand().equals(EDIT)) {
            editImage3();
        }

        else if (e.getSource() == iTrack) {
            setTrack();
        }
        else if (e.getSource() == iSister) {
            handleSisterRequest();
        }
        else if (e.getSource() == iColorToggle) {
            toggleColor();
        }
        if (doUpdate) updateDisplay();
    }


    private void toggleColor() {
        iColor = (iColor + 1) % 4;
    }

    public int getColor() {
        return iColor;
    }

    // handle track/no track button action
    private void setTrack() {
        if (iTrackPosition != ImageWindow.NONE) {
            iTrackPositionSave = iTrackPosition;
            iTrackPosition = ImageWindow.NONE;
        } else {
            iTrackPosition = iTrackPositionSave;
        }
    }

    public void forceTrackingOn() {
        iTrackPosition = ImageWindow.ANTERIOR;
        iTrackPositionSave = ImageWindow.POSTERIOR;
    }

    public void controlCallback(Vector v) {
        iImgWin.setSpecialEffect(null);
        String ctrl = (String)v.elementAt(0);
        boolean haveTime = false;
        boolean haveCellName = false;
        boolean haveCellIndex = false;
        if (ctrl.equals("InputCtrl1")) {
            //println("controlCallback: ");
        	//requestFocus();
            String time = ((String)v.elementAt(1)).trim();

            int requestedTime = -1;
            Vector v2 = null;
            try {
                requestedTime = Integer.parseInt(time);
                v2 = (Vector)iNucleiMgr.getNucleiRecord().elementAt(requestedTime - 1);
                haveTime = true;

            } catch(Exception e) {
                //System.out.println("bad image time: " + time);
                //return;
            }
            String cell = ((String)v.elementAt(2)).trim();
	    //System.out.println("controlCallback: " + cell + CS + time);
            boolean numeric = false;
            if (cell.length() > 0) {
                numeric = Character.isDigit(cell.charAt(0));
                if (numeric) haveCellIndex = true;
                else haveCellName = true;

            }
            //System.out.println("controlCallback: " + time + CS + cell + CS + numeric);
            boolean valid = haveCellName || haveTime;
            if (!valid) {
                println("controlCallback: invalid choice");
                return;
            }

            Cell c = null;
            if (numeric) {
                int index = Integer.parseInt(cell);
                try {
                    cell = ((Nucleus)v2.elementAt(index - 1)).identity;
                    c = (Cell)iAncesTree.getCellsByName().get(cell);
                    c = null;
                } catch(ArrayIndexOutOfBoundsException aiob) {
                    System.out.println("bad cell index: " + cell);
                    return;
                } catch(Exception e) {
                    System.out.println("ControlCallback: " + e);
                    return;
                }
            } else if (cell.length() > 0) {

                c = (Cell)iAncesTree.getCellsByName().get(cell);
            } else c = null;

            Cell csave = null;
            if (c != null) {
                if (!haveTime) requestedTime = c.getEnd();
            } else if (haveTime) {
                // try for a cell near the middle of the embryo
                for (int i=0; i < v2.size(); i++) {
                    Nucleus n = (Nucleus)v2.get(i);
                    c = null;
                    if (n.status > 0) {
                        c = (Cell)iAncesTree.getCellsByName().get(n.identity);
                        if (csave == null) csave = c;
                        if (n.z > 15 && n.z < 18) break;
                    }
                }
                if (c == null) c = csave;
            }
            setCurrentCell(c, requestedTime, CONTROLCALLBACK);
        }
    }

    public void controlCallback(Vector v, int oldversion) {
        iImgWin.setSpecialEffect(null);
        String ctrl = (String)v.elementAt(0);
        if (ctrl.equals("InputCtrl1")) {
            String time = ((String)v.elementAt(1)).trim();

            int requestedTime = -1;
            Vector v2 = null;
            try {
                requestedTime = Integer.parseInt(time);
                //v2 = (Vector)iNucleiMgr.getNucleiRecord()[requestedTime - 1];
                v2 = (Vector)iNucleiMgr.getNucleiRecord().elementAt(requestedTime - 1);

            } catch(Exception e) {
                System.out.println("bad image time: " + time);
                return;
            }
            String cell = ((String)v.elementAt(2)).trim();
            //System.out.println("controlCallback: " + time + CS + cell);
            boolean numeric = Character.isDigit(cell.charAt(0));
            Cell c = null;
            if (numeric) {
                int index = Integer.parseInt(cell);
                try {
                    cell = ((Nucleus)v2.elementAt(index - 1)).identity;
                    c = (Cell)iAncesTree.getCellsByName().get(cell);
                } catch(ArrayIndexOutOfBoundsException aiob) {
                    System.out.println("bad cell index: " + cell);
                    return;
                } catch(Exception e) {
                    System.out.println("ControlCallback: " + e);
                    return;
                }
            } else {
                c = (Cell)iAncesTree.getCellsByName().get(cell);
            }
            if (c == null) {
                System.out.println("bad cell name: " + cell);
                return;
            }
            //System.out.println("controlCallback: " + c.getName() + CS + requestedTime);
            setCurrentCell(c, requestedTime, CONTROLCALLBACK);
        }
    }

    public Cell getCellByName(String name) {
        if (iAncesTree == null) return null;
    	else return (Cell)iAncesTree.getCellsByName().get(name);
    }

    public Cell getCellByHash(String hashKey) {
        return (Cell)iAncesTree.getCells().get(hashKey);
    }

    public AncesTree getAncesTree() {
        return iAncesTree;
    }

    public JTree getJTree() {
        return iTree;
    }

    public NucleiMgr getNucleiMgr() {
        return iNucleiMgr;
    }

    private boolean nextTime() {
        //System.out.println("AceTree.nextTime: " + iImageTime + CS + iTimeInc + CS + iEndingIndex);
        //new Throwable().printStackTrace();
        if (iImageTime + iTimeInc == iEndingIndex) return false;
        iTimeInc++;
        if (iImage3D != null) iImage3D.insertContent(getImageTitle());
        if (iImage3D2 != null) iImage3D2.insertContent(getImageTitle());
        if (iImage3D2Z != null) iImage3D2Z.insertContent(getImageTitle());
        iCallSaveImage = true;
        int now = iImageTime + iTimeInc;
        int end = 999;
        if (iCurrentCell != null) end = iCurrentCell.getEnd();
        if (now <= end) return true; // we will call updateDisplay next
        if (iCurrentCell.getFateInt() == Cell.DIED) {
            iCurrentCellPresent = false;
            //System.out.println("cell death -- tracking cannot continue on this cell: " + iCurrentCell);
            //System.out.println("nextTime turning tracking off");
            iImageTime += iTimeInc;
            iTimeInc = 0;
            iTrackPosition = ImageWindow.NONE;
            return true;
        }
        // at this point we know that a cell division occurred in this transition
        // iCurrentCell will change as a side effect of doDaughterDisplayWork
        setCurrentCell(iCurrentCell, now, NEXTTIME);
        /*
        */
        return true;
    }

    public boolean prevTime() {
        if (iImageTime + iTimeInc <= iStartingIndex) return false;
        iTimeInc--;
        if (iImage3D != null) iImage3D.insertContent(getImageTitle());
        if (iImage3D2 != null) iImage3D2.insertContent(getImageTitle());
        if (iImage3D2Z != null) iImage3D2Z.insertContent(getImageTitle());
        iCallSaveImage = true;
        int now = iImageTime + iTimeInc;
        int start = 0;
        if (iCurrentCell != null) start = iCurrentCell.getTime();
        if (now >= start) return true;
        // a cell change occurs as we move to parent here
        //println("prevTime: " + iCurrentCell.getName() + CS + now);
        setCurrentCell(iCurrentCell, now, PREVTIME);
        return true;
    }


    private void showTreeCell(Cell c) {
        iIgnoreValueChanged = true;
        TreeNode [] tna = c.getPath();
        TreePath tp = new TreePath(tna);
        //iTree.expandPath(tp);
        iTree.makeVisible(tp);
        int row = iTree.getRowForPath(tp);
        iTree.setSelectionInterval(row,row);
        iTree.scrollRowToVisible(row);
        iTree.makeVisible(tp);
        //iIgnoreValueChanged = false;
    }

    private void makeDaughterDisplay(Cell c) {
        iTimeInc = 0;
        iPlaneInc = 0;
        getTimeAndPlane(c);
        if (iImageTime < 1 || iImagePlane < 1) return;
        doDaughterDisplayWork((Cell)c.getParent(), c);
    }

    private void doDaughterDisplayWork(Cell parent, Cell selectedDaughter) {
        //Cell parent = (Cell)c.getParent();
        //System.out.println("doDaughterDisplayWork: " + parent + CS + selectedDaughter);
        if (parent == null) System.out.println("*******NULL PARENT");
        if (!isTracking()) return;
        if (parent.getName() == ROOTNAME) return;
        if (iTimeInc != 0) return;
        int k = parent.getChildCount();
        if (k <= 1) return;
        Cell anteriorCell = (Cell)parent.getChildAt(0);
        Cell posteriorCell = (Cell)parent.getChildAt(1);
        //System.out.println("makeDaughterDisplay: " + anteriorCell + CS + posteriorCell);
        if (selectedDaughter != null) {
            if (selectedDaughter == anteriorCell) iTrackPosition = ImageWindow.ANTERIOR;
            else iTrackPosition = ImageWindow.POSTERIOR;
        }
        Cell save = iCurrentCell;
        if (iTrackPosition == ImageWindow.ANTERIOR) iCurrentCell = anteriorCell;
        else iCurrentCell = posteriorCell;
        if (iCurrentCell == null) {
            iCurrentCell = save;
            return;
        }

        //iNoTrack.setText(NOTRACK);

        //Vector nuclei = iNucleiMgr.getNucleiRecord()[iImageTime + iTimeInc - 1];
        Vector nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(iImageTime + iTimeInc - 1);
        String currentName = parent.getName();
        StringBuffer dummy = new StringBuffer();
        Nucleus anterior = NucUtils.getCurrentCellData(nuclei, anteriorCell.getName());
        Nucleus posterior = NucUtils.getCurrentCellData(nuclei, posteriorCell.getName());
        //System.out.println("makeDaughterDisplay: anterior: " + anterior);
        //System.out.println("makeDaughterDisplay: posterior: " + posterior);
        if (anterior != null && posterior != null) {
            makeAndSetSpecialEffects(anterior, posterior);
        }
    }

    private void makeAndSetSpecialEffects(Nucleus anterior, Nucleus posterior) {
        //System.out.println("makeAndSetSpecialEffects1: " + anterior);
        //new Throwable().printStackTrace();
        //System.out.println("makeAndSetSpecialEffects1: " + posterior);
        Object [] se = new Object[8];
        if (iTrackPosition == ImageWindow.ANTERIOR) {
            se[0] = new Integer(anterior.x);
            se[1] = new Integer(anterior.y);
            se[2] = new Integer((int)(anterior.z + HALFROUND));
            se[3] = new Integer(posterior.x);
            se[4] = new Integer(posterior.y);
            se[5] = new Integer(posterior.size/2);
            se[6] = new Integer((int)(posterior.z + HALFROUND));
            se[7] = posterior.identity;
        } else {
            se[0] = new Integer(posterior.x);
            se[1] = new Integer(posterior.y);
            se[2] = new Integer((int)(posterior.z + HALFROUND));
            se[3] = new Integer(anterior.x);
            se[4] = new Integer(anterior.y);
            se[5] = new Integer(anterior.size/2);
            se[6] = new Integer((int)(anterior.z + HALFROUND));
            se[7] = anterior.identity;
        }
        iImgWin.setSpecialEffect(se);

    }


    private void trackingActionsOnCurrentCellChange() {
        //System.out.println("trackingActionsOnCurrentCellChange");
        // set iImageTime and iTimeInc cleanly
        // set iImagePlane and iPlaneInc cleanly
        // assume initially that the transition was to a previous time
        int time = iImageTime + iTimeInc;
        iImageTime = iCurrentCell.getTime();
        iTimeInc = time - iImageTime;
        int plane = iImagePlane + iPlaneInc;
        iImagePlane = iCurrentCell.getPlane();
        iPlaneInc = plane - iImagePlane;
    }

    public void setCurrentCell(Cell c, int time, int source) {
        //println("AceTree.setCurrentCell: " + c + CS + time + CS + source);
        if (c == null) {
        	if (source == CONTROLCALLBACK) showSelectedCell(c, time);
        	return;
        }
        //println("AceTree.setCurrentCell: " + c.getName() + CS + time + CS + source);
        iImgWin.setSpecialEffect(null);
        //System.out.println("setCurrentCell: " + c + CS + time + CS + source);
        if (source != RIGHTCLICKONEDITIMAGE && !iCellsByName.containsKey(c.getName())) {
            //System.out.println("setCurrentCell:2 " + c + CS + time + CS + source);
            return;
        }
        if (source == RIGHTCLICKONIMAGE) {
            Cell old = iCurrentCell;
            iCurrentCell = c; //(Cell)iAncesTree.getCellsByName().get(cellName);
            trackingActionsOnCurrentCellChange();
            iAceTree.forceTrackingOn();
            showTreeCell(iCurrentCell);
            String s = makeDisplayText();
            iText.setText(s);
            //println("setCurrentCell:3 " + iCurrentCell + CS + old);
            iImgWin.updateCurrentCellAnnotation(iCurrentCell, old, -1);
            Cell parent = (Cell)c.getParent();
            if (iTimeInc == 0 && !parent.getName().equals("P0")) doDaughterDisplayWork((Cell)c.getParent(), c);
            else {
                if (iCurrentCell.isAnterior()) iTrackPosition = ImageWindow.ANTERIOR;
                else iTrackPosition = ImageWindow.POSTERIOR;
            }

            updateDisplay();
        } else if (source == RIGHTCLICKONEDITIMAGE) {
            //println("setCurrentCell:4 ");
            Cell old = iCurrentCell;
            iCurrentCell = c;
            trackingActionsOnCurrentCellChange();
            iAceTree.forceTrackingOn();
            showTreeCell(iCurrentCell);
            String s = "added cell in progress";
            iImgWin.updateCurrentCellAnnotation(iCurrentCell, old, -1);
            iText.setText(s);
            updateDisplay();

        } else if (source == LEFTCLICKONTREE) {
            showSelectedCell(c, time);
            updateDisplay();
            if (iImage3D != null) iImage3D.insertContent(getImageTitle());

        } else if (source == RIGHTCLICKONTREE) {
            //System.out.println("setCurrentCell RIGHTCLICKONTREE: " + c.getName() + CS + time);
            if (c.isAnterior()) iTrackPosition = ImageWindow.ANTERIOR;
            else iTrackPosition = ImageWindow.POSTERIOR;
            showSelectedCell(c, time);
            if (iImage3D != null) iImage3D.insertContent(getImageTitle());
        } else if (source == CONTROLCALLBACK) {
            showSelectedCell(c, time);
        } else if (source == NEXTTIME) {
            iImageTime = time;
            iTimeInc = 0;
            Cell currentCellSave = iCurrentCell;
            doDaughterDisplayWork(iCurrentCell, null);
            if (currentCellSave != iCurrentCell) {
                trackingActionsOnCurrentCellChange();
                iImgWin.updateCurrentCellAnnotation(iCurrentCell, currentCellSave, time);
            }
            showTreeCell(iCurrentCell);
        } else if (source == PREVTIME) {
            Vector nuclei1 = (Vector)iNucleiMgr.getNucleiRecord().elementAt(iImageTime + iTimeInc);
            Vector nuclei0 = (Vector)iNucleiMgr.getNucleiRecord().elementAt(iImageTime + iTimeInc - 1);
            Nucleus n = NucUtils.getParent(nuclei0, nuclei1, iCurrentCell.getName());
            Cell currentCellSave = iCurrentCell;
            if (n != null) {
                iCurrentCell = (Cell)iAncesTree.getCellsByName().get(n.identity);
                if (iCurrentCell == null) {

                    iCurrentCell = currentCellSave;
                    return;
                }
                if (currentCellSave != iCurrentCell) {
                    trackingActionsOnCurrentCellChange();
                    iImgWin.updateCurrentCellAnnotation(iCurrentCell, currentCellSave, time);
                }
                showTreeCell(iCurrentCell);
            } else {
                iTrackPosition = ImageWindow.NONE;
                iCurrentCell = null;
                showTreeCell(iRoot);
            }

        }
    }

    /**
     * does the work required by the cell selection control
     * @param c Cell the cell desired
     * @param requestedTime int the time index where it is to be shown
     * @param v2 vector of nuclei at this time point
     */
    private void showSelectedCell(Cell c, int requestedTime) {
        if (c == null) {
            iImageTime = requestedTime;
            iTimeInc = 0;
            iImagePlane = 15;
            iPlaneInc = 0;
            iCurrentCell = iRoot;
            showTreeCell(iCurrentCell);
            updateDisplay();

        	return;
        }

        String name = c.getName();
        Nucleus n = iNucleiMgr.getCurrentCellData(name, requestedTime);
        //if (n == null) {
        //    System.out.println("cell " + c + " not present at time " + requestedTime);
        //    return;
        //}
        if (n != null) {
            Cell old = iCurrentCell;
            iImageTime = c.getTime();
            iTimeInc = requestedTime - iImageTime;
            iImagePlane = (int)(n.z + HALFROUND);
            iPlaneInc = 0;
            iCurrentCell = c;
            //System.out.println("showSelectedCell: " + iCurrentCell + CS + c + CS + iImagePlane + CS + iPlaneInc);
            //if (iImageTime < 1 || iImagePlane < 1) return;
            if (iImageTime < 1) return;
            iCurrentCellPresent = true;
            if (iCurrentCell.isAnterior()) iTrackPosition = ImageWindow.ANTERIOR;
            else iTrackPosition = ImageWindow.POSTERIOR;
            showTreeCell(iCurrentCell);

        //if (iTimeInc == 0) makeDaughterDisplay(iCurrentCell);

            int baseTime = c.getTime(); //Integer.parseInt(sa[0]);
            iImgWin.updateCurrentCellAnnotation(iCurrentCell, old, -1);
            updateDisplay();
        }
        else {
            iImageTime = requestedTime;
            iTimeInc = 0;
            iImagePlane = 15;
            iPlaneInc = 0;
            iCurrentCell = c;
            showTreeCell(iCurrentCell);
            updateDisplay();

        }

    }



    private void incTime(int inc) {
        if (inc > 0) {
            if (iImageTime + iTimeInc < iEndingIndex) iTimeInc ++;
        }
        else if (iImageTime + iTimeInc > iStartingIndex) iTimeInc--;
    }

    private void incPlane(int inc) {
    	//println("AceTree.incPlane, " + inc);
        if (inc > 0) {
            if (iImagePlane + iPlaneInc < iPlaneEnd) iPlaneInc ++;
        }
        else if (iImagePlane + iPlaneInc > 1) iPlaneInc --;
    }

    private void handleSisterRequest() {
        //if (iImgWin.getSpecialEffect() == null)
        int k = iCurrentCell.getSiblingCount();
        if (k != 2) return;
        Cell parent = (Cell)iCurrentCell.getParent();
        Cell anteriorCell = (Cell)parent.getChildAt(0);
        Cell posteriorCell = (Cell)parent.getChildAt(1);
        int now = iImageTime + iTimeInc;
        Nucleus anterior = iNucleiMgr.getCurrentCellData(anteriorCell.getName(), now);
        Nucleus posterior = iNucleiMgr.getCurrentCellData(posteriorCell.getName(), now);
        if (anterior == null || posterior == null) {
            System.out.println("sister not present");
            return;
        }
        makeAndSetSpecialEffects(anterior, posterior);
    }

    /**
     * called from ImageWindow when user clicks a nucleus
     * @param e MouseEvent detected in ImageWindow
     */
    public void mouseMoved(MouseEvent e) {
        String s = POSITION + e.getX() + ", " + e.getY();
       iText2.setText(s);
    }

    public void cellAnnotated(String name) {
        iText3.setText(name);
    }
    ///////////////////// editing ///////////////////////////////////

    public void saveNuclei(File file) {
        System.out.println("saveNuclei: " + file);
        //iEditLog.showMe();
        NucZipper nz = new NucZipper(file, iNucleiMgr);
        nz = null;
        //iEditLog.setModified(false);
    }

    public void viewNuclei() {
        new NucEditDialog(this, iMainFrame, false);
    }

    //public void addNucleus() {
        //new NucAddDialog(false);
    //    new AddNucToRoot(false);
    //}

    public void relinkNucleus() {
        int time = iImageTime + iTimeInc;
        iNucRelinkDialog = new NucRelinkDialog(this, iMainFrame, false, iCurrentCell, time);
    }

    public void killCell(int x) {
    	println("killCell, ");
        //if (iTimeInc != 0 && iPlaneInc != 0) return;
        //Vector nuclei = iNucleiMgr.getNucleiRecord()[iImageTime + iTimeInc - 1];
        Vector nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(iImageTime + iTimeInc - 1);
        String name = iCurrentCell.getName();
        Nucleus n = null;
        for (int j=0; j < nuclei.size(); j++) {
            n = (Nucleus)nuclei.elementAt(j);
            if (!n.identity.equals(name)) continue;
            n.status = Nucleus.NILLI;
            break;
        }
        prevTime();
        updateDisplay();

        //clearTree();
        //buildTree(true);

    }

    public void killDeepNucs() {
    	new KillDeepNucsDialog(this, iMainFrame, true);
    }

    public void killDeepNucs(int zLim) {
        Vector nucRec = (Vector)iNucleiMgr.getNucleiRecord();
        for (int i=0; i < nucRec.size(); i++) {
        	Vector nuclei = (Vector)nucRec.get(i);
        	for (int j=0; j < nuclei.size(); j++) {
        		Nucleus n = (Nucleus)nuclei.get(j);
        		if (n.status == Nucleus.NILLI) continue;
        		if (n.z < zLim) continue;
        		println("killDeepNucs, " + i + CS + n);
        		n.status = Nucleus.NILLI;
        	}
        }
        clearTree();
        buildTree(true);

    }

    public void testWindow() {
    	new TestWindow(this, iMainFrame, false);
    }

    public void killCells() {
        int time = iImageTime + iTimeInc;
        new KillCellsDialog(this, iMainFrame, true, iCurrentCell, time, iEditLog);
    }

    public void pausePlayerControl() {
        iPlayerControl.pause();
    }

    public void setEndTime() {
        new SetEndTimeDialog(this, iMainFrame, true);
    }

    public void incrementEndTime() {
        setEndingIndex(++iEndingIndex);
    }

    public void setEndingIndex(int endTime) {
        iEndingIndex = endTime;
        iNucleiMgr.setEndingIndex(endTime);
        clearTree();
        Hashtable oldHash = iAncesTree.getCellsByName();
        buildTree(true);
        Hashtable newHash = iAncesTree.getCellsByName();
        String name = null;
        Cell c = null;
        Enumeration newKeys = newHash.keys();
        Vector newNames = new Vector();
        while(newKeys.hasMoreElements()) {
            name = (String)newKeys.nextElement();
            if(oldHash.containsKey(name)) continue;
            c = (Cell)newHash.get(name);
            newNames.add(name);
        }
        Collections.sort(newNames);
        for (int i=0; i < newNames.size(); i++) {
            c = (Cell)newHash.get((String)newNames.elementAt(i));
            //System.out.println(c.toString(0));
        }
    }

    public void undo() {
        iEditLog.append("UNDO");
        iNucleiMgr.restoreNucleiRecord();
        iNucleiMgr.clearAllHashkeys();
        clearTree();
        buildTree(true);
        setStartingCell((Cell)iRoot.getFirstChild(), 1);
        iEditLog.setModified(true);

    }

    ///////////////////// editing end ///////////////////////////////////

    public void exit() {

       //JOptionPane pane = new JOptionPane(message, JOptionPane.QUESTION_MESSAGE),
        //         JOptionPane.OK_CANCEL_OPTION);
        //pane.set.Xxxx(...); // Configure
        if (iEditLog.getModified()) {
            Object[] options = { "OK", "CANCEL" };
            String msg = "Warning, you have unsaved edits.\nClick OK to continue.";
            int choice = JOptionPane.showOptionDialog(null, msg, "Warning",
                     JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE,
                             null, options, options[0]);
            if (choice == 0) iWinEvtHandler.windowClosing(null);
        } else {
            iWinEvtHandler.windowClosing(null);
        }
    }

    public void editTraverse() {
        iEditTraverse = new EditTraverse(iCurrentCell);
    }

    public void setEditTraverseNull() {
        iEditTraverse = null;
    }

    public void setFocusHome() {
        iHome.grabFocus();
    }

    public void ancestral() {

        //iCurrentCell.setLateTime(250);

        // jtb: hacking this so that the initial tree displays
        if (iCurrentCell.toString().equals("P")) {
          new SulstonTree(this, "Ancestral Tree", new Cell("AB"), true);
        }
        else {
          new SulstonTree(this, "Ancestral Tree", iCurrentCell, true);
        }
        //iShowTriangle = !iShowTriangle;
    }

    public void canonical() {
        //System.out.println("AceTree.test");
        if (iCanonicalTree == null) iCanonicalTree = CanonicalTree.getCanonicalTree();
        new AuxFrame(this, "Sulston Tree", iCanonicalTree);
    }

    public void vtree() {
	    new VTree();
    }

    public void test() {
        //System.out.println("AceTree.test");
        new AceTreeHelp("/org/rhwlab/help/messages/J3Derror.html", 600, 300);

    }

    public CanonicalTree getCanonicalTree() {
        return iCanonicalTree;
    }

    public void threeDview() {
        String s = getImageTitle();
        try {
            iImage3D = new Image3D(this, s);
            iImage3D.insertContent(s);
        } catch(NoClassDefFoundError ee) {
            System.out.println("you need to install Java3D");
            new AceTreeHelp("/org/rhwlab/help/messages/J3Derror.html", 600, 300);
        }
    }

    public void threeDview2() {
        String s = getImageTitle();
        try {
            iImage3D2 = new Image3D2(this, s);
            //iImage3D2.insertContent(s);
        } catch(NoClassDefFoundError ee) {
            System.out.println("you need to install Java3D");
            new AceTreeHelp("/org/rhwlab/help/messages/J3Derror.html", 600, 300);
        }

    }

    public void threeDview2Z() {
        String s = getImageTitle();
        try {
            iImage3D2Z = new Image3D2Z(this, s);
            //iImage3D2.insertContent(s);
        } catch(NoClassDefFoundError ee) {
            System.out.println("you need to install Java3D");
            new AceTreeHelp("/org/rhwlab/help/messages/J3Derror.html", 600, 300);
        }

    }


    public void allCentroidsView() {
        String s = getImageTitle();
        new ImageAllCentroids(this, s);

    }

    public String getImageTitle() {
        String s = makeImageName();
        int k = s.lastIndexOf("-");
        s = s.substring(0, k);
        String s2 = iTifPrefix;
        int j = s2.lastIndexOf(C.Fileseparator);
        if (j > 0) s2 = s2.substring(j + 1);
        return s2 + s;
    }

    public void image3DOff() {
        iImage3D = null;
    }

    public void image3DSave(boolean saveIt) {
        Image3D.setSaveImageState(saveIt);
        if (!saveIt) return;
        if (iImage3D == null) {
            System.out.println("no active image3d");
            threeDview();
            //return;
        } else {
            iImage3D.saveImage();
        }
        //image3D2Save(saveIt);
    }

    public void image3D2Save(boolean saveIt) {
        Image3D2.setSaveImageState(saveIt);
        if (!saveIt) return;
        if (iImage3D2 == null) {
            System.out.println("no active image3d2");
            threeDview2();
            //return;
        } else {
//            iImage3D2.saveImage();
        	println("image3D2Save, " + iImgWin.getTitle());
        	println("image3D2Save, " + iImgWin.getName());
        	String name = iImgWin.getTitle();

// jtb: not sure of why these are here; replacing them with
// something else which is hopefully equivalent.
//        	name = name.substring(0, name.length() - 8);
//        	name = name.substring(4, name.length());
          int i = name.lastIndexOf(java.io.File.separatorChar);
          name = name.substring(i+1, name.length() - 8);
System.out.println("name is now " + name);
            iImage3D2.offScreenRendering(name);
        }
    }

    public void image3D2ZSave(boolean saveIt) {
        Image3D2Z.setSaveImageState(saveIt);
        if (!saveIt) return;
        if (iImage3D2Z == null) {
            System.out.println("no active image3d2");
            threeDview2Z();
            //return;
        } else {
            iImage3D2Z.saveImage();
        }
    }

    public void image2DSave(boolean saveIt) {
        iImgWin.setSaveImageState(saveIt);
        iImgWin.saveImageIfEnabled();
    }


    private void delay(int n) {
        long t = System.currentTimeMillis();
        long e = t + n;
        while (t < e) t = System.currentTimeMillis();

    }

    public boolean nextImage() {
        iImgWin.setSpecialEffect(null);
        boolean b = nextTime();
        updateDisplay();
        return b;
    }

    public boolean prevImage() {
        boolean b = prevTime();
        updateDisplay();
        return b;
    }

    public void updateUseZip(int useZip) {
    	iUseZip = useZip;
    	println("updateUseZip, " + useZip);
    	ImageWindow.cUseZip = useZip;
    	iNucleiMgr.getConfig().iUseZip = useZip;
    }

    public void exportNewick() {
        try {
            Newick newick = new Newick(iRoot);
        } catch(Exception e) {
            System.out.println("ATVapp unavailable");
        } catch(NoClassDefFoundError ee) {
            System.out.println("you need to get ATVTree.jar");
            new AceTreeHelp("/org/rhwlab/help/messages/ATVerror.html");
        }
    }
/*
    public Parameters getParameters() {
        //System.out.println("getParameters: " + iParameters);
        return iParameters;
    }
*/
    public EditLog getEditLog() {
        return iEditLog;
    }

    public Log getDebugLog() {
        return iDLog;
    }

    public int getTimeInc() {
        return iTimeInc;
    }

    public int getPlaneInc() {
        return iPlaneInc;
    }

    public int getImageTime() {
        return iImageTime;
    }

    public int getImagePlane() {
        return iImagePlane;
    }

    public boolean getShowCentroids() {
        return iShowCentroids;
    }

    public void setShowCentroids(boolean show) {
        iShowCentroids = show;
    }

    public boolean getShowAnnotations() {
        return iShowAnnotations;
    }

    public void setShowAnnotations(boolean show) {
        //println("setShowAnnotations: " + show);
        //new Throwable().printStackTrace();
        iShowAnnotations = show;
        if (iShow != null) {
            if (show) iShow.setText(HIDE);
            else iShow.setText(SHOW);
        }
    }

    public Cell getCurrentCell() {
        return iCurrentCell;
    }

    public JFrame getMainFrame() {
        return iMainFrame;
    }

    public Object getDispProps3D() {
        return iDispProps3D;
    }

    public Object getDispProps3D2() {
        return iDispProps3D2;
    }

    public Object getDispProps3D2Z() {
        return iDispProps3D2Z;
    }

    public void setDispProps3D(Object obj) {
        iDispProps3D = obj;
    }

    public void setDispProps3D2(Object obj) {
        iDispProps3D2 = obj;
    }

    public void setDispProps3D2Z(Object obj) {
        iDispProps3D2Z = obj;
    }

    public Cell getRoot() {
        return iRoot;
    }

    public boolean isTracking() {
        //System.out.println("isTracking: " + iTrackPosition + CS + (iTrackPosition != ImageWindow.NONE));
        return iTrackPosition != ImageWindow.NONE;
    }

    public void setOrientation(String orientation) {
        iOrientation = orientation;
        System.out.println("setOrientation: " + iOrientation);
    }

    public String getOrientation() {
        System.out.println("getOrientation: " + iOrientation);
        return iOrientation;
    }

    public Hashtable getNucleiMgrHash() {
        return iNucleiMgrHash;
    }

    public ImageWindow getImageWindow() {
        return iImgWin;
    }

    private void debugShow(int testTime) {
        System.out.println();
        System.out.println("setup for edit");
        Vector nuclei = null;
        nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(testTime - 1);
        Enumeration e = nuclei.elements();
        System.out.println("time: " + (testTime - 1));
        while (e.hasMoreElements()) {
            System.out.println((Nucleus)e.nextElement());
        }

        nuclei = (Vector)iNucleiMgr.getNucleiRecord().elementAt(testTime);
        e = nuclei.elements();
        System.out.println("time: " + (testTime));
        while (e.hasMoreElements()) {
            System.out.println((Nucleus)e.nextElement());
        }

    }



    public final static int
         LEFTCLICKONTREE = 1
        ,RIGHTCLICKONTREE = 2
        ,RIGHTCLICKONIMAGE = 3
        ,CONTROLCALLBACK = 4
        ,NEXTTIME = 5
        ,PREVTIME = 6
        ,RIGHTCLICKONEDITIMAGE = 7
        ;

    private final static int
     WIDTH = 330
    ,HEIGHT200 = 200
    ,HEIGHT100 = 100
    ,HEIGHT30 = 30
    ,XINC = 8
    ,YINC = 12
   ;


    final public static String
     PARAMETERS = "parameters"
    ,POSITION = "Mouse position: "
    ,SPACES15 = "               "
    ,TITLE = "AceTree"
    ,HELPMSG = "you must provide file: "
    ,SEP = ", "
    ,ROOTNAME = "P"
    ;

    final public static String
    NEXTT = "Next"
   ,PREV = "Prev"
   ,UP   = "Up  "
   ,DOWN = "Down"
   ,HOME = "Home"
   ,SHOW = "ShowA"
   ,SHOWC = "ShowC"
   ,HIDE = "HideA"
   ,HIDEC = "HideC"
   ,CLEAR = "Clear"
   ,COPY = "Copy"
   ,EDIT = "Edit"
   ,NOTRACK = "NoTrk"
   ,TRACK = "Track"
   ,SISTER = "Sister"
   ,COLORTOGGLE = "All/G/R/N"
   ,CS = ", "
   ;

    /*
    private static final String [] configParams = {
            "zipFileName"
           ,"tif directory"
           ,"tifPrefix"
           ,"nuclei directory"
           ,"root name"
           ,"starting index"
           ,"ending index"
           ,"use zip"
           ,"namingMethod"
    };
    */

    private static final int
         ZIPFILENAME = 0
        ,TIFDIRECTORY = 1
        ,TIFPREFIX = 2
        ,NUCLEIDIRECTORY = 3
        ,ROOTNAMEI = 4
        ,STARTINGINDEX = 5
        ,ENDINGINDEX = 6
        ,USEZIP = 7
        ,NAMINGMETHOD = 8
        ;

    private static final float
         HALFROUND = 0.5f
        ;

    public void debugTest(boolean b) {
        iDebugTest = b;
    }

    private void createAndShowGUI() {
        JFrame.setDefaultLookAndFeelDecorated(true);
        //iMainFrame = new JFrame(TITLE);
        iMainFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        //JPanel newContentPane = this;
        iMainFrame.getContentPane().setLayout(new BorderLayout());
        //int hInput = iInputCtrl.getY();
        //System.out.println("hInput="+hInput);
        int height = HEIGHT200 + HEIGHT100 + HEIGHT100 + HEIGHT100 + 2*HEIGHT30;
        this.setMinimumSize(new Dimension(WIDTH, height));
        this.setOpaque(true); //content panes must be opaque

        // make this AceTree instance the content pane
        //iAceMenuBar = new AceMenuBar(this);
        iMainFrame.setJMenuBar(iAceMenuBar);
        iMainFrame.getContentPane().add(this, BorderLayout.NORTH);
        iMainFrame.pack();
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
        Dimension windowSize = iMainFrame.getSize();
        //System.out.println("windowSize: " + windowSize);
        iMainFrame.setLocation(Math.max(0,(screenSize.width -windowSize.width)/2),
                Math.max(0,(screenSize.height-windowSize.height)/2));
        iWinEvtHandler = new WindowEventHandler();
        iMainFrame.addWindowListener(iWinEvtHandler);

        iMainFrame.setFocusTraversalPolicy(
                new FocusControl(iHome
                        , iInputCtrl.getTimeField()
                        , iInputCtrl.getNameField()
                        ));
        iMainFrame.setVisible(true);
    }

    public void run(String arg0) {
        createAndShowGUI();
    }

    Orientation		iOrientationPanel;
    void showOrientation() {
    	iOrientationPanel = new Orientation();
    }

    Zafer1			iZafer1Panel;
    void showZafer1() {
    	iZafer1Panel = new Zafer1();
    }

    Juvenesence		iJuvenesencePanel;
    void showJuvenesence() {
    	iJuvenesencePanel = new Juvenesence();
    }

    Lazarus			iLazarusPanel;
    void showLazarus() {
    	iLazarusPanel = new Lazarus();
    }

    Siamese			iSiamesePanel;
    void showSiamese() {
    	iSiamesePanel = new Siamese();
    }

    DeathsAdjacencies 	iDeathsAdjacenciesPanel;
    void showDeathsAdjacencies() {
    	iDeathsAdjacenciesPanel = new DeathsAdjacencies();
    }

///////////////////////////////////////////////////////////////////////////////
///// stuff that gets modified for the AceTree object in EmbryoDB ///////////////////////
/////////////////////////////////////////////////////////////////////////////
    private boolean				iCmdLineRun;

    private class WindowEventHandler extends WindowAdapter {
        public void windowActivated(WindowEvent e) {
            iHome.requestFocusInWindow();
        }
        public void windowClosing(WindowEvent e) {
            System.out.println("AceTree shutdown " + new GregorianCalendar().getTime() + CS + iCmdLineRun);
            iMainFrame.dispose();
            if (iImgWin != null) iImgWin.dispose();
            if (iCmdLineRun) System.exit(0);
        }

    }

///////////////////////////////////////////////////////////////////////////////
///// end of stuff that gets modified for the AceTree object in EmbryoDB ///////////////////////
/////////////////////////////////////////////////////////////////////////////
    public static void main(String[] args) throws IOException {
    	
    	System.out.println("Working Directory = " +
            System.getProperty("user.dir"));
    
        System.out.println("AceTree launched: " + new GregorianCalendar().getTime());
        ManifestX.reportAndUpdateManifest();
        String config = null;
        AceTree ot;
        //System.exit(1);
        if (args.length > 0) {
            System.out.println("AceTree args[0]: " + args[0]);
            ot = getAceTree(args); //new AceTree(args[0]);

        }
        else ot = getAceTree(null);
        ot.run("");
        ot.debugTest(false);
        ot.iCmdLineRun = true;
        System.out.println("main exiting");
    }

    private static void println(String s) {System.out.println(s);}

    public void medianSmooth() {
    	// jtb: trying to add smoothing
    	MedianSmoother smoother
    		= new MedianSmoother(iNucleiMgr.getRoot(), 30);
    	System.out.println("about to smooth");
    	smoother.computeMedians();
    	System.out.println("finished smoothing");
    }

}
