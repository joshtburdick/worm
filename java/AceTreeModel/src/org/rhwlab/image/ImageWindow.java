/*
 * Copyright 2005 University of Washington Genome Sciences
 * All rights reserved
 */
package org.rhwlab.image;
import java.awt.AWTException;
import java.awt.BorderLayout;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.HeadlessException;
import java.awt.Image;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Robot;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.awt.event.WindowListener;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.Vector;
import java.util.zip.ZipEntry;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.OvalRoi;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.io.TiffDecoder;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.stream.FileImageOutputStream;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.border.Border;
import javax.swing.event.MouseInputAdapter;

import net.sf.ij.jaiio.BufferedImageCreator;

import org.rhwlab.acetree.AceTree;
import org.rhwlab.acetree.AnnotInfo;
import org.rhwlab.acetree.NucUtils;
import org.rhwlab.help.AceTreeHelp;
import org.rhwlab.image.Image3D.SublineageDisplayProperty;
import org.rhwlab.image.Image3D.PropertiesTab.SublineageUI;
import org.rhwlab.nucedit.AddOneDialog;
import org.rhwlab.nucedit.EIDialog1;
//import org.rhwlab.nucedit.EIDialog2;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.Cell;
import org.rhwlab.utils.C;
import org.rhwlab.utils.EUtils;

/**
 * Provides a JFrame window to contain the ImageJ ImagePlus object
 *
 * @author biowolp
 * @version 1.0 January 25, 2005
 */
public class ImageWindow extends JFrame
        implements ActionListener, KeyListener, Runnable {
    public ImageCanvas             iImgCanvas;
    static ImagePlus        iImgPlus;
    String                  iTitle;
    static Object []        iSpecialEffect;
    AceTree                 iAceTree;
    Vector                  iAnnotsShown;
    MouseHandler            iMouseHandler;
    boolean                 iMouseEventHandled;
    int                     iImageTime;
    int                     iTimeInc;
    int                     iImagePlane;
    int                     iPlaneInc;
    boolean                 iIsMainImgWindow;
    boolean                 iIsRightMouseButton;
    boolean                 iSaveImage;
    boolean                 iSaveInProcess;
    String                  iSaveImageDirectory;
    boolean                 iUseRobot;
    boolean                 iNewConstruction;
    //private JTabbedPane     iTabbedPane;
    public static ColorSchemeDisplayProperty []     iDispProps;
    private JPanel          iControlPanel;
    protected JMenuBar      iMenuBar;
    protected JToolBar      iToolBar;
    protected JButton	    iHelp;
    protected JButton       iProperties;

    static boolean         cAcbTree = false;

    static byte []          iRpix;
    static byte []          iGpix;
    static byte []          iBpix;



    // static variables and functions

    public static String        cZipTifFilePath;
    public static String        cTifPrefix;
    public static String        cTifPrefixR;
    public static int           cUseZip;
    static ZipImage             cZipImage;
    public static NucleiMgr     cNucleiMgr;
    public static int           cImageWidth;
    public static int           cImageHeight;
    public static int           cLineWidth;
    public static String        cCurrentImageFile;
    public static String        cCurrentImagePart;

    //public static EditImage3    cEditImage3;

    /**
     * this is the constructor that is actually used
     * note that there are many static functions and class variables
     */
    public ImageWindow(String title, ImagePlus imgPlus) {
        super(title);
        iTitle = title;
        iImgPlus = imgPlus;
        ImageCanvas ic = new ImageCanvas(imgPlus);
        iImgCanvas = ic;
        iDispProps = getDisplayProps();

        iToolBar = new JToolBar();
        iToolBar.setLayout(new GridLayout(1,0));
        iToolBar.setMinimumSize(new Dimension(ImageWindow.cImageWidth, 30));
        iToolBar.setPreferredSize(new Dimension(ImageWindow.cImageWidth, 30));
        iHelp = new JButton("Help");
        iHelp.setFocusable(false);
        iHelp.addActionListener(this);
        iToolBar.add(iHelp);
        iProperties = new JButton("Properties");
        iProperties.addActionListener(this);
        iProperties.setFocusable(false);
        iToolBar.add(iProperties);

        fillToolBar();

        Container c = getContentPane();
        JPanel jp = new JPanel();
        jp.setLayout(new BorderLayout());
        jp.add(iToolBar, BorderLayout.NORTH);
        jp.add(ic, BorderLayout.SOUTH);
        c.add(jp);

        pack();
        setVisible(true);
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        WinEventMgr wem = new WinEventMgr();
        addWindowFocusListener(wem);
        addWindowListener(wem);
        iMouseHandler = new MouseHandler(this);
        iImgCanvas.addMouseMotionListener(iMouseHandler);
        iImgCanvas.addMouseListener(iMouseHandler);
        setImageTimeAndPlaneFromTitle();
        iAnnotsShown = new Vector();
        iIsRightMouseButton = false;
        iSaveImage = false;
        iSaveImageDirectory = null;
        iUseRobot = false;

        iImgCanvas.addKeyListener(this);

        setDialogsEnabled(false);
    }


    public static void setStaticParameters(String zipTifFilePath, String tifPrefix, int useZip) {
        //System.out.println("ImageWindow.setStaticParameters entered");
        cZipTifFilePath = zipTifFilePath;
        cTifPrefix = tifPrefix;
        cUseZip = useZip;
        if (cUseZip == 1) cZipImage = new ZipImage(cZipTifFilePath);
        cLineWidth = 2;//LINEWIDTH;
        String [] sa = cTifPrefix.split("/");
        if(sa.length > 1) cTifPrefixR = sa[0] + "R" + C.Fileseparator + sa[1];
        //System.out.println("cZipTifFilePath, cTifPrefix, cTifPrefixR: " +
        //        cZipTifFilePath + CS + cTifPrefix + CS + cTifPrefixR);
    }
    public static void setNucleiMgr(NucleiMgr nucleiMgr) {
        cNucleiMgr = nucleiMgr;
    }

    public static ImagePlus makeImage(String s) {
        //System.out.println("makeImage: " + s + CS + cUseZip);
        cCurrentImageFile = s;
        //new Throwable().printStackTrace();
        ImagePlus ip = null;
        //iSpecialEffect = null;
        //if (iSpecialEffect != null) iSpecialEffect = null;

        //println("makeImage, cUseZip=" + cUseZip);
        switch(cUseZip) {
        case 0:
        case 3:
            ip = doMakeImageFromTif(s);
            break;
        case 1:
            ip = doMakeImageFromZip(s);
            break;
        default:
            ip = doMakeImageFromZip2(s);
            break;

        }
        /*
        if (cUseZip == 0) {
            ip = doMakeImageFromTif(s);
        } else {
            ip = doMakeImageFromZip(s);
        }
        */
        //if (cImageWidth == 0 && ip != null) {
        if (ip != null) {
            cImageWidth = ip.getWidth();
            cImageHeight = ip.getHeight();
            //System.out.println("***ImageWindow: " + cImageWidth + CS + cImageHeight);
        }
        if (ip == null) return iImgPlus;
        else return ip;
    }

    public static ImagePlus doMakeImageFromZip(String s) {
        //System.out.println("ImageWindow.doMakeImageFromZip entered: " + s);
        if (cZipImage == null) cZipImage = new ZipImage(cZipTifFilePath);
        ZipEntry ze = cZipImage.getZipEntry(s);
        ImagePlus ip;
        if (ze == null) {
            ip = new ImagePlus();
            ImageProcessor iproc = new ColorProcessor(cImageWidth, cImageHeight);
            ip.setProcessor(s, iproc);
        }
        else ip = cZipImage.readData(ze);
        //System.out.println("ImageWindow.makeImage exiting");
        return ip;
    }

    public static ImagePlus doMakeImageFromZip2(String s) {
        //System.out.println("ImageWindow.doMakeImageFromZip2 entered: " + s);
        //System.out.println("ImageWindow.doMakeImageFromZip2: " + cZipTifFilePath);
        //if (cZipImage == null) cZipImage = new ZipImage(cZipTifFilePath + File.separator + s);
        cZipImage = new ZipImage(cZipTifFilePath + "/" + s);
        int k1 = s.indexOf("/") + 1;
        String ss = s.substring(k1);
        int k2 = ss.indexOf(".");
        ss = ss.substring(0, k2);
        //System.out.println("using: " + ss);
        ZipEntry ze = null;
        if (cZipImage != null) ze = cZipImage.getZipEntry(ss + ".tif");
        //System.out.println("ZipEntry: " + ze);
        //if (cZipImage == null) cZipImage = new ZipImage(cZipTifFilePath);
        //ZipEntry ze = cZipImage.getZipEntry(s);
        ImagePlus ip;
        if (ze == null) {
            ip = new ImagePlus();
            ImageProcessor iproc = new ColorProcessor(cImageWidth, cImageHeight);
            ip.setProcessor(s, iproc);
        }
        else ip = cZipImage.readData(ze);
        //System.out.println("ImageWindow.makeImage exiting");
        //ip = convertToRGB(ip);
        ColorProcessor iprocColor = (ColorProcessor)ip.getProcessor();
        int [] all = (int [])iprocColor.getPixels();
        byte [] R = new byte[all.length];
        byte [] G = new byte[all.length];
        byte [] B = new byte[all.length];
        //ColorProcessor iproc3 = new ColorProcessor(iproc.getWidth(), iproc.getHeight());
        iprocColor.getRGB(R, G, B);
        //G = bpix;
        //R = getRedChannel(R);
        iRpix = R;
        iGpix = G;
        iBpix = B;
        return ip;
    }

    private static void showError(String fileName) {
	new Throwable().printStackTrace();
        String message = "Exiting: cannot find image\n";
        message += fileName;
        JOptionPane pane = new JOptionPane(message);
        JDialog dialog = pane.createDialog(null, "Error");
        dialog.show();
    }

    public static ImagePlus doMakeImageFromTif(String s) {
		if (cUseZip == 3) s = s.replaceAll("tif", "jpg");
        //println("ImageWindow.doMakeImageFromTif entered: " + s);

        cCurrentImagePart = s;
        //FileInputStream fis;
        ImagePlus ip = null;
        String ss = cZipTifFilePath + C.Fileseparator + s;
        //println("ImageWindow.makeImage entered: " + ss);
        ip = new Opener().openImage(ss);
        if (ip != null) {
            cImageWidth = ip.getWidth();
            cImageHeight = ip.getHeight();
            ip = convertToRGB(ip);
        } else {
            ip = new ImagePlus();
            ImageProcessor iproc = new ColorProcessor(cImageWidth, cImageHeight);
            ip.setProcessor(s, iproc);
        }

        return ip;
    }

    public static ImagePlus readData(FileInputStream fis, boolean bogus) {
        if (fis == null) return null;
        int byteCount;
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        byte[] buf = new byte[4096];
        try {
            InputStream is = (InputStream)fis;
            byte data[] = new byte[DATA_BLOCK_SIZE];

            //  4. read source zipped data and write to uncompressed stream
            while ( (byteCount = is.read(data, 0, DATA_BLOCK_SIZE)) != -1) {
                out.write(data, 0, byteCount);
            }
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }
        return openTiff(new ByteArrayInputStream(out.toByteArray()), true);

    }

    public static ImagePlus readData(FileInputStream fis) {
        if (fis == null) return null;
        byte [] ba = readByteArray(fis);
        return openTiff(new ByteArrayInputStream(ba), true);
    }

    public static byte[] readByteArray(FileInputStream fis) {
        if (fis == null) return null;
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        int byteCount;
        byte[] buf = new byte[4096];
        try {
            InputStream is = (InputStream)fis;
            byte data[] = new byte[DATA_BLOCK_SIZE];

            //  4. read source zipped data and write to uncompressed stream
            while ( (byteCount = is.read(data, 0, DATA_BLOCK_SIZE)) != -1) {
                out.write(data, 0, byteCount);
            }
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }
        return out.toByteArray();

    }

    /** Attempts to open the specified inputStream as a
    TIFF, returning an ImagePlus object if successful. */
    public static ImagePlus openTiff(InputStream in, boolean convertToRGB) {
        //System.out.println("openTiff entered");
        if (in == null) return null;
        FileInfo[] info = null;
        try {
            TiffDecoder td = new TiffDecoder(in, null);
            info = td.getTiffInfo();
        } catch (FileNotFoundException e) {
            IJ.error("TiffDecoder", "File not found: "+e.getMessage());
            return null;
        } catch (Exception e) {
            IJ.error("TiffDecoder", ""+e);
            return null;
        }
        ImagePlus imp = null;
        if (IJ.debugMode) // dump tiff tags
            IJ.log(info[0].info);
        FileOpener fo = new FileOpener(info[0]);
        imp = fo.open(false);
        // detect 8 bit or RGB from the FileInfo object info[0]
        if (info[0].getBytesPerPixel() == 1 && convertToRGB) {
            imp = convertToRGB(imp);
        }
        //IJ.showStatus("");
        return imp;
    }

    /**
     * If the images in the zip archive are 8 bit tiffs,
     * we use that as the green plane of an RGB image processor
     * so the program is always showing RGB images
     *
     * @param ip an Image processor obtained from the image file
     * @return
     */
    private static ImagePlus convertToRGB(ImagePlus ip) {
        //System.out.println("convertToRGB entered");
        ImageProcessor iproc = ip.getProcessor();
        byte [] bpix = (byte [])iproc.getPixels();
        byte [] R = new byte[bpix.length];
        byte [] G = new byte[bpix.length];
        byte [] B = new byte[bpix.length];
        ColorProcessor iproc3 = new ColorProcessor(iproc.getWidth(), iproc.getHeight());
        iproc3.getRGB(R, G, B);
        // special test removal
        G = bpix;
        R = getRedChannel(R);
        // special test insert
        //byte [] ba = getRedChannel(R);
        //System.arraycopy(ba, 0, R, 0, ba.length);
        //System.arraycopy(ba, 0, G, 0, ba.length);
        //System.arraycopy(ba, 0, B, 0, ba.length);
        //R = ba;
        //G = ba;
        //B = ba;

        // end special
        iRpix = R;
        iGpix = G;
        iBpix = B;
        return buildImagePlus(ip);
        //iproc3.setRGB(R, G, B);
        //ip.setProcessor("test", iproc3);
        //return ip;
    }

    private static ImagePlus buildImagePlus(ImagePlus ip) {
        ImageProcessor iproc = ip.getProcessor();
        ColorProcessor iproc3 = new ColorProcessor(iproc.getWidth(), iproc.getHeight());
        iproc3.setRGB(iRpix, iGpix, iBpix);
        ip.setProcessor("test", iproc3);
        return ip;

    }

    protected static ImagePlus makeRedImagePlus(ImagePlus ip) {
        ImageProcessor iproc = ip.getProcessor();
        ColorProcessor iproc3 = new ColorProcessor(iproc.getWidth(), iproc.getHeight());
        iproc3.setRGB(iRpix, new byte[iRpix.length], new byte[iRpix.length]);
        ip.setProcessor("test", iproc3);
        return ip;
    }

    protected static ImagePlus makeGreenImagePlus(ImagePlus ip) {
        ImageProcessor iproc = ip.getProcessor();
        //System.out.println("makeGreenImagePlus: " + iproc);
        ColorProcessor iproc3 = new ColorProcessor(iproc.getWidth(), iproc.getHeight());
        //System.out.println("makeGreenImagePlus2: " + iproc + CS + iGpix  + CS + iRpix);
        iproc3.setRGB(new byte[iRpix.length], iGpix, new byte[iRpix.length]);
        ip.setProcessor("test", iproc3);
        return ip;
    }

    protected static ImagePlus makePlainImagePlus(ImagePlus ip) {
        ImageProcessor iproc = ip.getProcessor();
        ColorProcessor iproc3 = new ColorProcessor(iproc.getWidth(), iproc.getHeight());
        if (cAcbTree) {
            //byte [] added = new byte[iRpix.length];
            //for (int i=0; i < iRpix.length; i++) {
             //   added[i] = (byte)(iRpix[i] + iGpix[i]);
            //}
            iproc3.setRGB(iRpix, iRpix, iRpix);
        } else {
            iproc3.setRGB(new byte[iRpix.length], new byte[iRpix.length], new byte[iRpix.length]);

        }
        ip.setProcessor("test", iproc3);
        return ip;
    }



    private static byte[] getRedChannel(byte [] R) {
        String fileName = makeRedChannelName();
        //System.out.println("getRedChannel: " + fileName);
        File f = new File(fileName);
        if (f.exists()) {
            FileInputStream fis;
            ImagePlus ip = null;
            ip = new Opener().openImage(fileName);
            if (ip != null) {
                ByteProcessor bproc = (ByteProcessor)ip.getProcessor();
                R = (byte [])bproc.getPixels();
            } else {
                System.out.println("getRedChannel, Opener returned null ip");
            }
        } else {
            //System.out.println("getRedChannel, file does not exist");

        }
        return R;

    }

    private static String makeRedChannelName() {
        // 20071108 rehacked this because windows vista was very picky
        // and backslashes were plagueing me
        // the green parsing was working so I created cCurrentImagePart
        // to go from there to red by substituting "tifR" for "tif"
        String s = cCurrentImageFile;
        //int k = s.indexOf(cTifPrefix) + cTifPrefix.length();
        String ss = cCurrentImagePart;
        //System.out.println("getRedChannelName, " + ss);
        ss = ss.substring(3);
        //System.out.println("getRedChannelName, " + ss);
        //s = cZipTifFilePath + C.Fileseparator + cTifPrefixR + s.substring(k);
        s = cZipTifFilePath + C.Fileseparator + "/tifR/" + ss;
        if (cUseZip == 3) s = cZipTifFilePath + C.Fileseparator + "/jpgR/" + ss;

        return s;
    }

    // end of static stuff

    public ImageWindow() {

    }

    /**
     * this constructor for test purposes only
     */
    public ImageWindow(String title, ImagePlus imgPlus, boolean test) {
        super(title);
        iTitle = title;
        iImgPlus = imgPlus;
        ImageCanvas ic = new ImageCanvas(imgPlus);
        iImgCanvas = ic;
        getContentPane().add(ic);
        pack();
        setVisible(true);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    }

    public class ColorSchemeDisplayProperty {
        public String iName;
        public int    iLineageNum;

        public ColorSchemeDisplayProperty(String name, int lineageNum) {
            iName = name;
            iLineageNum = lineageNum;
        }
    }

    public class ColorSchemeUI {
        public JPanel       iPanel;
        public JTextField   iTF;
        public JComboBox    iCB;
        public JLabel       iLabel;

        public ColorSchemeUI(int i) {
            iPanel = new JPanel();
            iPanel.setLayout(new GridLayout(1,2));
            iTF = new JTextField(iDispProps[i].iName, WIDTH);
            iLabel = new JLabel(iDispProps[i].iName);
            String [] list;
            list = COLORS;
            if (i == 5) list = SIZES;
            iCB = new JComboBox(list);
            iCB.setSelectedIndex(iDispProps[i].iLineageNum);
            //iPanel.add(iTF);
            iPanel.add(iLabel);
            iPanel.add(iCB);
            iPanel.setMaximumSize(new Dimension(200,10));
        }
        private String [] COLORS = {
                "red"
                ,"blue"
                ,"green"
                ,"yellow"
                ,"cyan"
                ,"magenta   "
                ,"pink"
                ,"gray"
                ,"white"

        };

        private String [] SIZES = {"1", "2", "3"};

    }


    public class PropertiesTab implements ActionListener {
        JPanel                          iPanel;
        //SublineageDisplayProperty []    iDispProps;
        ColorSchemeUI []                 iCSUI;

        public PropertiesTab() {
            Border blackline = BorderFactory.createLineBorder(Color.black);
            iDispProps = getDisplayProps();
            iCSUI = new ColorSchemeUI[iDispProps.length];
            iPanel = new JPanel();
            iPanel.setLayout(new BorderLayout());
            iPanel.setBorder(blackline);
            //iPanel.setLayout(new BoxLayout(iPanel, BoxLayout.PAGE_AXIS));
            JPanel lineagePanel = new JPanel();
            JPanel dummyPanel = new JPanel();
            JPanel topPart = new JPanel();
            topPart.setLayout(new GridLayout(1,2));
            lineagePanel.setLayout(new GridLayout(0,1));
            lineagePanel.setBorder(blackline);
            //lineagePanel.setMaximumSize(new Dimension(300,400));
            topPart.add(lineagePanel);
            topPart.add(dummyPanel);
            JPanel [] testPanel = new JPanel[iDispProps.length];
            JTextField textField;
            JComboBox cb;
            JPanel labelPanel = new JPanel();
            JLabel sublineage = new JLabel("item");
            JLabel color = new JLabel("color");
            labelPanel.setLayout(new GridLayout(1,2));
            labelPanel.add(sublineage);
            labelPanel.add(color);
            lineagePanel.add(labelPanel);

            for (int i=0; i < iDispProps.length; i++) {
                iCSUI[i] = new ColorSchemeUI(i);
                lineagePanel.add(iCSUI[i].iPanel);
            }
            lineagePanel.setMaximumSize(new Dimension(200, 200));
            iPanel.add(topPart, BorderLayout.NORTH);
            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new GridLayout(1,3));
            //buttonPanel.setMinimumSize(new Dimension(400, 100));

            JButton reset = new JButton("Reset");
            JButton apply = new JButton("Apply");
            JButton cancel = new JButton("Cancel");
            buttonPanel.add(apply);
            reset.addActionListener(this);
            apply.addActionListener(this);
            cancel.addActionListener(this);
            buttonPanel.add(reset);
            buttonPanel.add(apply);
            buttonPanel.add(cancel);
            JPanel botPart = new JPanel();
            botPart.setLayout(new GridLayout(5,1));
            botPart.add(new JPanel());
            botPart.add(buttonPanel);
            botPart.add(new JPanel());
            botPart.add(new JPanel());
            botPart.add(new JPanel());
            iPanel.add(botPart, BorderLayout.CENTER);

            //iPanel.add(buttonPanel, BorderLayout.CENTER);
            //iPanel.add(new JPanel(), BorderLayout.CENTER);

        }

        public void actionPerformed(ActionEvent e) {
            String command = e.getActionCommand();
            if (command.equals("Reset")) {
                iDispProps = getDisplayProps();
                for (int i=0; i < iDispProps.length; i++) {
                    iCSUI[i].iLabel.setText(iDispProps[i].iName);
                    iCSUI[i].iCB.setSelectedIndex(iDispProps[i].iLineageNum);
                }


            } else if (command.equals("Apply")) {
                for (int i=0; i < iDispProps.length; i++) {
                    String name = iCSUI[i].iTF.getText();
                    if (name.length() == 0) name = "-";
                    int num = iCSUI[i].iCB.getSelectedIndex();
                    iDispProps[i].iName = name;
                    iDispProps[i].iLineageNum = num;
                }

            }
        }

        public JPanel getPanel() {
            return iPanel;
        }

        private String [] COLORS = {
                "red"
                ,"blue"
                ,"green"
                ,"yellow"
                ,"cyan"
                ,"magenta   "
                ,"pink"
                ,"gray"
                ,"white"

        };

        private String [] SIZES = {"1", "2", "3"};

        private static final int
            WIDTH = 15
           ;

    }

    public ColorSchemeDisplayProperty [] getDisplayProps() {
        ColorSchemeDisplayProperty [] dispProps = {
                new ColorSchemeDisplayProperty("normal centroid", 1)
                ,new ColorSchemeDisplayProperty("selected centroid", 8)
                ,new ColorSchemeDisplayProperty("annotations", 8)
                ,new ColorSchemeDisplayProperty("upper sister", 4)
                ,new ColorSchemeDisplayProperty("lower sister", 5)
                ,new ColorSchemeDisplayProperty("line size" , 1)
        };
        return dispProps;
    }

    private int getLineageNumber(String name) {
        int num = iDispProps.length;
        for (int i=0; i < iDispProps.length; i++) {
            if (name.indexOf(iDispProps[i].iName) >= 0) {
                num = iDispProps[i].iLineageNum;
                break;
            }
        }
        return num;
    }






    protected JMenuBar createMenuBar() {
        JMenuBar menuBar = new JMenuBar();
        //JMenu menu = new JMenu(FILE);
        //menuBar.add(menu);
        //JMenuItem test = new JMenuItem(SAVEAS);
        //menu.add(test);
        //test.addActionListener(this);
        JMenu menu = null;
        JMenuItem test = null;

        menu = new JMenu("dummy");
        menuBar.add(menu);
        test = new JMenuItem("dummy");
        menu.add(test);
        return menuBar;

    }

    public void setAceTree(AceTree aceTree) {
        iAceTree = aceTree;
    }

    public AceTree getAceTree() {
        return iAceTree;
    }

    public ImagePlus refreshDisplay(String imageName) {
    	//println("refreshDisplay, ");
        if (imageName == null) imageName = iTitle;
        else {
            if (imageName.indexOf(cTifPrefix) == -1) {
                imageName = cTifPrefix +imageName;
            }
            iTitle = imageName;
            setTitle(iTitle);
        }
        //System.out.println("ImageWindow.refreshDisplay2: " + imageName);
        if (iIsMainImgWindow) {
            iTimeInc = iAceTree.getTimeInc();
            iPlaneInc = iAceTree.getPlaneInc();
            iImageTime = iAceTree.getImageTime();
            iImagePlane = iAceTree.getImagePlane();
        } else {
            iTimeInc = 0;
            iPlaneInc = 0;
            setImageTimeAndPlaneFromTitle();
        }
        String random = RANDOMT;
        if (cUseZip > 0) random = RANDOMF;
        int k = imageName.indexOf(random);
        if (k > -1) imageName = imageName.substring(0, k + random.length() - 1 );
        ImagePlus ip = null;

        //System.out.println("ImageWindow.refreshDisplay3: " + System.currentTimeMillis());
        //new IOException().printStackTrace();
        ip = makeImage(imageName);
        if (ip == null) {
            iAceTree.pausePlayerControl();
            System.out.println("no ip for: " + iTitle);
        }

        if (iAceTree == null) return null;

        switch (iAceTree.getColor()) {
            case 1:
                ip = makeGreenImagePlus(ip);
                break;
            case 2:
                ip = makeRedImagePlus(ip);
                break;
            case 3:
                ip = makePlainImagePlus(ip);
                break;
            default:
        }
        //ip = makeGreenImagePlus(ip);

        if (ip != null) iImgPlus.setProcessor(imageName, ip.getProcessor());
        if (iIsMainImgWindow && iAceTree.isTracking()) iAceTree.addMainAnnotation();
        if (iAceTree.getShowCentroids()) showCentroids();
        if (iAceTree.getShowAnnotations()) showAnnotations();
        if (iSpecialEffect != null) showSpecialEffect();
        //iSpecialEffect = null;
        iImgCanvas.repaint();
        return iImgPlus;

    }

    /* (non-Javadoc)
     * @see java.awt.event.KeyListener#keyPressed(java.awt.event.KeyEvent)
     */
    public void keyPressed(KeyEvent e) {
        int code = e.getKeyCode();
        int mods = e.getModifiers();
        boolean shift = (mods & InputEvent.SHIFT_MASK) == InputEvent.SHIFT_MASK;
        boolean ctrl = (mods & InputEvent.CTRL_MASK) == InputEvent.CTRL_MASK;
        println("ImageWindow.keyPressed, " + code + CS + shift + CS + ctrl + CS + e);
        if (shift || ctrl) sendToEIDialog2(code, shift, ctrl);
        else {
        switch(code) {
            case KeyEvent.VK_UP:
                iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.UP));
                break;
            case KeyEvent.VK_DOWN:
                iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.DOWN));
                break;
            case KeyEvent.VK_LEFT:
                iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.PREV));
                break;
            case KeyEvent.VK_RIGHT:
                iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.NEXTT));
                break;
            case KeyEvent.VK_F2:
                iAceTree.actionPerformed(new ActionEvent(this, 0, "F2"));
                break;
            default:
                return;

        }
        }
    }

    private void sendToEIDialog2(int keycode, boolean alt, boolean ctrl) {
    	println("sendToEIDialog2, ");
    	ActionEvent a = null;
        switch(keycode) {
        case KeyEvent.VK_UP:
            iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.UP));
            break;
        case KeyEvent.VK_DOWN:
            iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.DOWN));
            break;
        case KeyEvent.VK_LEFT:
        	//if (ctrl) a = new ActionEvent(this, 0, EIDialog2.LEFT);
        	//else a = new ActionEvent(this, 0, EIDialog2.BIG);
        	AddOneDialog addOne = iAceTree.iAddOneDialog;
            if (addOne != null) addOne.actionPerformed(new ActionEvent(this, 0, AddOneDialog.LEFT));
            break;
        case KeyEvent.VK_RIGHT:
            iAceTree.actionPerformed(new ActionEvent(this, 0, AceTree.NEXTT));
            break;
        case KeyEvent.VK_F2:
            iAceTree.actionPerformed(new ActionEvent(this, 0, "F2"));
            break;
        default:
            return;

    }

    }

    public void quickRefresh() {
        iImgCanvas.repaint();
    }

    public void setSpecialEffect(Object [] specialEffect) {
        iSpecialEffect = specialEffect;
    }

    protected void showSpecialEffect() {
        if (!iAceTree.isTracking()) return;
        int x1 = ((Integer)iSpecialEffect[0]).intValue();
        int y1 = ((Integer)iSpecialEffect[1]).intValue();
        int z1 = ((Integer)iSpecialEffect[2]).intValue();
        int x2 = ((Integer)iSpecialEffect[3]).intValue();
        int y2 = ((Integer)iSpecialEffect[4]).intValue();
        int r2 = ((Integer)iSpecialEffect[5]).intValue();
        int z2 = ((Integer)iSpecialEffect[6]).intValue();
        String s = (String)iSpecialEffect[7];
        int offset = r2 + 4;
        if (y2 < y1) offset = -offset;


        ImageProcessor iproc = getImagePlus().getProcessor();
        //iproc.setColor(Color.magenta);
        iproc.setColor(COLOR[iDispProps[LOWERSIS].iLineageNum]);
        if (z2 <= z1) iproc.setColor(COLOR[iDispProps[UPPERSIS].iLineageNum]);
        //if (z2 <= z1) iproc.setColor(Color.cyan);

        iproc.setLineWidth(cLineWidth);
        //iproc.drawLine((int)x1, (int)y1, (int)x2, (int)y2);
        iproc.drawLine(x1, y1, x2, y2);
        iproc.drawPolygon(EUtils.pCircle((int)x2, (int)y2, (int)r2));
        iproc.drawString("    " + s + "(" + z2 + ")", x2, y2 + offset);
    }

    private void redrawMe() {
        iImgCanvas.repaint();
    }

    protected void setImageTimeAndPlaneFromTitle() {
        //System.out.println("setImage..: " + iTitle);
        int k = iTitle.lastIndexOf(DASHT) + DASHT.length();
        if (k <= 1) {
            iImageTime = 1;
            iImagePlane = 15;
            iTimeInc = 0;
            iPlaneInc = 0;
            String random = RANDOMT;
            if (cUseZip > 0) random = RANDOMF;
            iIsMainImgWindow = iTitle.indexOf(random) == -1;
            return;
        }
        System.out.println("setImage..: " + k);
        String time = iTitle.substring(k, k + 3);
        //System.out.println("setImage..: " + time);
        iImageTime = Integer.parseInt(time);
        String s = iTitle.substring(k);
        //System.out.println("setImage..: " + s);
        k = s.indexOf(DASHP) + DASHP.length();
        String plane = s.substring(k, k + 2);
        //System.out.println("setImage..: " + plane);
        iImagePlane = Integer.parseInt(plane);
        iTimeInc = 0;
        iPlaneInc = 0;
        String random = RANDOMT;
        if (cUseZip > 0) random = RANDOMF;
        iIsMainImgWindow = iTitle.indexOf(random) == -1;
    }


    public ImageCanvas getCanvas() {
        return iImgCanvas;
    }

    public ImagePlus getImagePlus() {
        return iImgPlus;
    }

    ////////////////////////////////////////

    public void addAnnotation(int mx, int my, boolean dontRemove) {
        if (iIsMainImgWindow) {
            iTimeInc = iAceTree.getTimeInc();
            iImageTime = iAceTree.getImageTime();
            iPlaneInc = iAceTree.getPlaneInc();
        } else {
            iTimeInc = 0;
            iPlaneInc = 0;
        }
        double x, y, r;
        boolean g;
        Nucleus n = cNucleiMgr.findClosestNucleus(mx, my, iImagePlane + iPlaneInc, iImageTime + iTimeInc);
            if (cNucleiMgr.hasCircle(n, (double)(iImagePlane + iPlaneInc))) {
                AnnotInfo ai = new AnnotInfo(n.identity, n.x, n.y);
                // now, if this one is not in the vector add it
                // otherwise remove it
                boolean itemRemoved = false;
                boolean itemAlreadyPresent = false;
                String test = n.identity;
                AnnotInfo aiTest = null;
                for (int k=0; k < iAnnotsShown.size(); k++) {
                    aiTest =(AnnotInfo)iAnnotsShown.elementAt(k);
                    if (aiTest.iName.equals(test)) {
                        itemAlreadyPresent = true;
                        if (!dontRemove) {
                            iAnnotsShown.remove(k);
                            itemRemoved = true;
                        }
                        break;
                    }

                }

                if (!itemRemoved && !itemAlreadyPresent) {
                    iAnnotsShown.add(ai);
                }
                // if this was a button 3 mouse click
                // and this is the main window
                // we will make this the current cell and makeDisplayText agree
                if (iIsRightMouseButton && iIsMainImgWindow) {
                    iIsRightMouseButton = false;
                }
            }



    }



    protected void showCentroids() {
        int time = iImageTime + iTimeInc;
        if (time < 0) {
            iImageTime = 1;
            iTimeInc = 0;
        }
        Vector v = (Vector)cNucleiMgr.getNucleiRecord().elementAt(iImageTime + iTimeInc - 1);
        ImageProcessor iproc = getImagePlus().getProcessor();
        iproc.setColor(COLOR[iDispProps[NCENTROID].iLineageNum]);
        iproc.setLineWidth(WIDTHS[iDispProps[LINEWIDTH].iLineageNum]);
        Polygon p = null;
        Enumeration e = v.elements();
        Cell currentCell = iAceTree.getCurrentCell();

        while(e.hasMoreElements()) {
            Nucleus n = (Nucleus)e.nextElement();
            if (n.status < 0) continue;
            double x = cNucleiMgr.nucDiameter(n,
                    (double)(iImagePlane + iPlaneInc));
            if (x > 0) {
                if (currentCell != null && n.hashKey != null && n.hashKey.equals(currentCell.getHashKey()) && iAceTree.isTracking()) {
                    iproc.setColor(COLOR[iDispProps[SCENTROID].iLineageNum]);
                }
                iproc.drawPolygon(EUtils.pCircle(n.x, n.y, (int)(x/2.)));
                iproc.setColor(COLOR[iDispProps[NCENTROID].iLineageNum]);
            }

        }
    }

    private void drawRoi(int plane, Nucleus c, ImageProcessor iproc) {
        double d = cNucleiMgr.nucDiameter(c, plane);
        float fxx = c.x;
        float fyy = c.y;
        fxx -= d/2;
        fyy -= d/2;
        int xx = (int)(fxx + 0.5);
        int yy = (int)(fyy + 0.5);
        int dd = (int)(d + 0.5);

        //int d = (int)(c.d + 0.5);
        //System.out.println("processImage, d=" + d + C.CS + c.d);
        //int xx = c.x - d/2;
        //int yy = c.y - d/2;
        OvalRoi oRoi = new OvalRoi(xx, yy, dd, dd);
        //Color csave = iproc.getColor();
        iproc.setColor(new Color(0, 0, 255));
        oRoi.drawPixels(iproc);
        Rectangle r = oRoi.getBounds();
        int width = iproc.getWidth();
        int offset, i;
        for (int y=r.y; y < (r.y + r.height); y++) {
            offset = y * width;
            for (int x = r.x; x <= (r.x + r.width); x++) {
                i = offset + x;
                if (oRoi.contains(x, y)) {
                    //iproc.drawPixel(x,y);
                    int k = iproc.getPixel(x, y);
                    int m = k & -16711936;
                    //System.out.println("drawRoi: " + k + C.CS + m);
                    //redSum += Math.abs(redPixels[i]);
                }
            }
        }



    }

    protected void showAnnotations() {
        //showWhichAnnotations();
        Vector v = (Vector)cNucleiMgr.getNucleiRecord().elementAt(iImageTime  + iTimeInc - 1);
        int size = v.size();
        int [] x = new int[size];
        int [] y = new int[size];
        Vector annots = new Vector();
        Enumeration e = v.elements();
        while(e.hasMoreElements()) {
            AnnotInfo ai = null;
            Nucleus n = (Nucleus)e.nextElement();
            //if (n.identity.length() > 0 && isInList(n.identity)) {
            if (n.status >= 0 && (isInList(n.identity) != null)) {
                ai = new AnnotInfo(n.identity, n.x, n.y);
                if (cNucleiMgr.hasCircle(n, (double)(iImagePlane + iPlaneInc))) {
                    annots.add(ai);
                }
            }
        }
        drawStrings(annots, this);
        //NucUtils.drawStrings(annots, this);
        //iShow.setText(HIDE);
    }

    private void drawStrings(Vector annots, ImageWindow imgWin) {
        ImagePlus imgPlus = imgWin.getImagePlus();
        ImageProcessor imgProc = imgPlus.getProcessor();
        ImageCanvas imgCan = imgWin.getCanvas();
        //imgProc.setColor(Color.yellow);
        //System.out.println("iDispProps: " + iDispProps);
        imgProc.setColor(COLOR[iDispProps[ANNOTATIONS].iLineageNum]);
        imgProc.setFont(new Font("SansSerif", Font.BOLD, 13));
        Enumeration e = annots.elements();
        while (e.hasMoreElements()) {
            AnnotInfo ai = (AnnotInfo)e.nextElement();
            imgProc.moveTo(imgCan.offScreenX(ai.iX),imgCan.offScreenY(ai.iY));
            imgProc.drawString(ai.iName);
        }
        imgPlus.updateAndDraw();
    }

    private void showWhichAnnotations() {
        for (int i=0; i < iAnnotsShown.size(); i++) {
            System.out.println((AnnotInfo)iAnnotsShown.elementAt(i));

        }

    }

    public void updateCurrentCellAnnotation(Cell newCell, Cell old, int time) {
        //new Throwable().printStackTrace();
        //println("updateCurrentCellAnnotation: " + newCell + CS + old + CS + time);
        AnnotInfo ai = null;
        if (old != null) ai = isInList(old.getName());
        if (ai != null) iAnnotsShown.remove(ai);
        if (time == -1) time = newCell.getTime();
        String s = newCell.getHashKey();
        Nucleus n = null;
        //println("updateCurrentCellAnnotation:2 " + s);
        if (s != null) {
            n = cNucleiMgr.getNucleusFromHashkey(newCell.getHashKey(), time);
            //println("updateCurrentCellAnnotation:3 " + n);
        }
        if ((n != null) && (isInList(newCell.getName()) == null)) {
            ai = new AnnotInfo(newCell.getName(), n.x, n.y);
            iAnnotsShown.add(ai);
        }
    }

    public void clearAnnotations() {
        iAnnotsShown.clear();
    }

    public void addAnnotation(String name, int x, int y) {
        AnnotInfo ai = new AnnotInfo(name, x, y);
        iAnnotsShown.add(ai);
    }

    protected AnnotInfo isInList(String name) {
        //System.out.println("isInList: " + name + CS + iAnnotsShown.size());
        AnnotInfo aiFound = null;
        Enumeration e = iAnnotsShown.elements();
        while(e.hasMoreElements()) {
            AnnotInfo ai = (AnnotInfo)e.nextElement();
            boolean is = ((String)ai.iName).equals(name);
            if (is) {
                aiFound = ai;
                break;
            }
        }
        return aiFound;
    }

    public void saveImageIfEnabled() {
        if (iSaveImage) {
            while(iSaveInProcess);
            new Thread(this).start();
        }
    }

    public void run() {
        iSaveInProcess = true;
        int k = 1000;
        if (iNewConstruction) {
            k = 5000; // long delay needed on new open
            iNewConstruction = false;
        }
        try {
            Thread.sleep(k);
        } catch(InterruptedException ie) {
            ie.printStackTrace();
        }
        saveImage();

    }

	void saveJpeg(BufferedImage bi, String outFileName, int quality) {
		Iterator iter = ImageIO.getImageWritersByFormatName("jpeg");
		javax.imageio.ImageWriter writer = (javax.imageio.ImageWriter)iter.next();
		ImageWriteParam iwp = writer.getDefaultWriteParam();
		iwp.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
		iwp.setCompressionQuality((float)quality/100);   // an integer between 0 and 1
		//println("otherSaveJpeg, " + iwp.canWriteCompressed() + CS + iwp.getCompressionQuality());
		//outFileName = "jpg" + quality + ".jpg";
		File file = new File(outFileName);
		if (file.exists()) file.delete();
		file = new File(outFileName);
		try {
			FileImageOutputStream output = new FileImageOutputStream(file);
			writer.setOutput(output);
			IIOImage image = new IIOImage(bi, null, null);
			writer.write(null, image, iwp);
			writer.dispose();
			output.close();
		} catch(Exception e) {
			e.printStackTrace();
		}
	}



    public void saveImage() {
        String title = makeTitle();
        if (title == null) {
            cancelSaveOperations();
            return;
        }
        Rectangle screenRect = this.getBounds();
        int topAdjust = 23;
        int y = screenRect.y;
        screenRect.y += topAdjust;
        int height = screenRect.height;
        screenRect.height -= topAdjust;
        // create screen shot
        //File f = new File("ij-ImageIO_.jar");
        //if (!f.exists()) {
        //    println("CANNOT SAVE FILES -- MISSING ij-ImageIO_.jar");
        //    return;
        //}

        Robot robot = null;
        BufferedImage image = null;
        if (iUseRobot) {
            try {
                robot = new Robot();
            } catch(AWTException e) {
                println("EXCEPTION -- NO ROBOT -- NOT SAVING");
                iSaveInProcess = false;
                iSaveImage = false;
                iAceTree.iAceMenuBar.resetSaveState();
                return;
            }
            image = robot.createScreenCapture(screenRect);
        } else {
            image = BufferedImageCreator.create((ColorProcessor)iImgPlus.getProcessor());
        }

        saveJpeg(image, title, 20);
        /*
        try {
            //robot = new Robot();
            //BufferedImage image = robot.createScreenCapture(screenRect);
            //BufferedImage image = BufferedImageCreator.create((ColorProcessor)iImgPlus.getProcessor());
            ImageIO.write(image, "jpeg", new File(title));
            //ImageIO.write(image, "png", new File(title));
        //} catch(AWTException awtex) {
        //    awtex.printStackTrace();
        } catch(Exception e) {
            e.printStackTrace();
        }
        */

        System.out.println("file: " + title + " written");
        iSaveInProcess = false;
    }

    public void cancelSaveOperations() {
        println("WARNING: NO IMAGE SAVE PATH -- NOT SAVING!");
        iSaveInProcess = false;
        iSaveImage = false;
        iAceTree.iAceMenuBar.resetSaveState();
        return;

    }

    public String getSaveImageDirectory() {
        if (iSaveImageDirectory != null) return iSaveImageDirectory;
        try {
            Class.forName("net.sf.ij.jaiio.BufferedImageCreator");
        } catch(ClassNotFoundException e) {
            iUseRobot = true;
            println("USING ROBOT FOR IMAGE2D SAVING");
        }

        {
            JFileChooser fc = new JFileChooser("");
            fc.setDialogTitle("Save images to: ");
            fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            int returnVal = fc.showOpenDialog(null);
            String path = null;
            if(returnVal == JFileChooser.APPROVE_OPTION) {
                path = fc.getSelectedFile().getPath();
                iSaveImageDirectory = path;
                System.out.println("Saving images to: " + path);
            }
            return path;
        }

    }

    private String makeTitle() {
        if (iSaveImageDirectory == null) {
            String dir = getSaveImageDirectory();
            iSaveImageDirectory = dir;
            if (dir == null) return null;
        }
        String s = iTitle;
        int j = s.lastIndexOf(C.Fileseparator) + 1;
        int k = s.lastIndexOf(".");
        s = s.substring(j, k) + ".jpeg";
        //s = s.substring(j, k) + ".png";
        s = iSaveImageDirectory + "/" + s;
        return s;
    }

    public void setSaveImageState(boolean saveIt) {
        iSaveImage = saveIt;
    }

    public Object [] getSpecialEffect() {
        return iSpecialEffect;
    }

    private class WinEventMgr extends WindowAdapter {
        public void windowGainedFocus(WindowEvent e) {
            //System.out.println("windowGainedFocus, ");
            //refreshDisplay(null);
        	iAceTree.requestFocus();

        }
        public void windowClosing(WindowEvent e) {
            System.out.println("windowClosing: " + iIsMainImgWindow);
            if (!iIsMainImgWindow) dispose();

        }
    }

    class MouseHandler extends MouseInputAdapter {

        public MouseHandler(ImageWindow iw) {
            super();
        }

        public void mouseMoved(MouseEvent e) {
            iAceTree.mouseMoved(e);
        }

        public void mouseClicked(MouseEvent e) {
        	//println("ImageWindow.mouseClicked, " + cEditImage3 + CS + e);
            int button = e.getButton();
            if (button == 3) {
                iIsRightMouseButton = true;
            } else {
                iIsRightMouseButton = false;
            }
            if (button == 3) {
                Nucleus n = cNucleiMgr.findClosestNucleus(e.getX(), e.getY(), iImagePlane + iPlaneInc, iImageTime + iTimeInc);
                if (n == null) return;
                Cell c = iAceTree.getCellByName(n.identity);
                iAceTree.setCurrentCell(c, iImageTime + iTimeInc, AceTree.RIGHTCLICKONIMAGE);
            } else if (button == 1){
                //System.out.println("mouseClicked " + e.getX());
                addAnnotation(e.getX(), e.getY(), false);
                refreshDisplay(null);
            }
            iAceTree.cellAnnotated(getClickedCellName(e.getX(), e.getY()));
            //if (cEditImage3 != null) cEditImage3.processEditMouseEvent(e);
            processEditMouseEvent(e);
        }

    }

    private String getClickedCellName(int x, int y) {
        int timeInc = 0;
        int planeInc = 0;
        if (iIsMainImgWindow) {
            timeInc = iAceTree.getTimeInc();
            planeInc = iAceTree.getPlaneInc();
        }
        String name = "";
        Nucleus n = cNucleiMgr.findClosestNucleus(x, y, iImageTime + iTimeInc);
        if (cNucleiMgr.hasCircle(n, (double)(iImagePlane + iPlaneInc))) {
            name = n.identity;
        }
        return name;
    }

    protected static final String
    RANDOMF = ".zip0"
   ,RANDOMT = ".tif0"
   ,DASHT = "-t"
   ,DASHP = "-p"
   ;

    public static void main(String[] args) {
        System.out.println("ImageWindow main");
        test3();
    }


    public static void test3() {
        FileInputStream fis;
        ImagePlus ip = null;
        //String ss = "/home/biowolp/AncesTree/temp2/images/050405-t050-p15.tif";
        String ss = "/nfs/waterston1/images/bao/081505/tif/081505_L1-t050-p15.tif";
        try {
            fis = new FileInputStream(ss);
            byte [] ba = readByteArray(fis);
            ip = openTiff(new ByteArrayInputStream(ba), false);
            //ip = readData(fis);
        } catch(IOException ioe) {
            System.out.println("ImageWindow.test3 exception ");
            System.out.println(ioe);
        }
        int width = ip.getWidth();
        int height = ip.getHeight();
        System.out.println("width, height: " + width + CS + height);
        ByteProcessor bproc = new ByteProcessor(width, height);
        //test(bproc);
        ImagePlus ip2 = new ImagePlus("newtest3", bproc);
        new ImageWindow("test", ip2);

        //ij.gui.ImageWindow iImgWin = new ij.gui.ImageWindow(ip2);

        //FileSaver fs = new FileSaver(ip2);
        //fs.saveAsTiff();

    }

    public static void test2() {
        ImageWindow.setStaticParameters("", "", 0);
        String s = "/home/biowolp/AncesTree/temp2/images/050405-t050-p15.tif";
        //String s = "/home/biowolp/0tmp/newtest2.tif";
        ImagePlus ip = ImageWindow.makeImage(s);
        System.out.println("handleImage: " + ip + CS + s);
        //ImageWindow iImgWin = new ImageWindow(s, ip, true);
        ij.gui.ImageWindow iImgWin = new ij.gui.ImageWindow(ip);
        ColorProcessor cproc = (ColorProcessor)ip.getProcessor();
        //byte [] bpix = (byte [])cproc.getPixels();
        int [] pix = (int [])cproc.getPixels();
        byte [] R = new byte[pix.length];
        byte [] G = new byte[pix.length];;
        byte [] B = new byte[pix.length];;
        cproc.getRGB(R, G, B);
        ByteProcessor bproc = new ByteProcessor(cproc.getWidth(), cproc.getHeight());
        bproc.setPixels(R);
        test(bproc);
        R = (byte [])bproc.getPixels();
        cproc.setRGB(R, G, B);

        ImagePlus ip2 = new ImagePlus("newtest2", cproc);

        FileSaver fs = new FileSaver(ip2);
        fs.saveAsTiff();


    }

    public static void test(ImageProcessor ip) {
        byte [] pixels = (byte [])ip.getPixels();
        int width = ip.getWidth();
        OvalRoi oRoi = new OvalRoi(560, 130, 55, 55);
        ip.setRoi(oRoi);
        Rectangle r = ip.getRoi();
        int offset, i;
        for (int y=r.y; y < (r.y + r.height); y++) {
            offset = y * width;
            for (int x = r.x; x < (r.x + r.width); x++) {
                i = offset + x;
                if (oRoi.contains(x, y)) {
                    pixels[i] = (byte)(128.*Math.random());
                }
            }
        }
    }

    public static final Integer
         ANTERIOR = new Integer(1)
        ,POSTERIOR = new Integer(2)
        ,NONE = new Integer(0)
        ;

    private static final String
         CS = ", "
        ;

    private static final int
    DATA_BLOCK_SIZE  = 2048
   //,LINEWIDTH = 1
   ;

    public static final int
         NCENTROID = 0
        ,SCENTROID = 1
        ,ANNOTATIONS = 2
        ,UPPERSIS = 3
        ,LOWERSIS = 4
        ,LINEWIDTH = 5
        ;

    public static final Color [] COLOR = {
            Color.RED
            ,new Color(140,70,255)
            ,Color.GREEN
            ,Color.YELLOW
            ,Color.CYAN
            ,Color.MAGENTA
            ,Color.PINK
            ,Color.LIGHT_GRAY
            ,Color.WHITE
    };

    public static final int [] WIDTHS = {1,2,3};

    private static void println(String s) {System.out.println(s);}


/*
    public void actionPerformed(ActionEvent e) {
        Object o = e.getSource();
        if (o == iProperties) {
            new ImageParmsDialog(this);
        } else if (o == iHelp) {
            String item = "/org/rhwlab/help/html/ImageWindowToolbarHelp.html";
            new AceTreeHelp(item);
        }


    }
*/

    /* (non-Javadoc)
     * @see java.awt.event.KeyListener#keyReleased(java.awt.event.KeyEvent)
     */
    public void keyReleased(KeyEvent e) {
        // TODO Auto-generated method stub

    }


    /* (non-Javadoc)
     * @see java.awt.event.KeyListener#keyTyped(java.awt.event.KeyEvent)
     */
    public void keyTyped(KeyEvent e) {
        // TODO Auto-generated method stub

    }

    private JButton     iAddSeries;
    private JButton     iAddOne;
    public JDialog 		iDialog;

    protected void fillToolBar() {
        iAddOne = new JButton("AddOne");
        iAddOne.addActionListener(this);
        iToolBar.add(iAddOne);
        iAddSeries = new JButton("AddSeries");
        iAddSeries.addActionListener(this);
        iToolBar.add(iAddSeries);

        iAddOne.setFocusable(false);
        iAddSeries.setFocusable(false);
        iAddOne.setFocusable(false);
    }

    public void actionPerformed(ActionEvent e) {
        Object o = e.getSource();
        String s = e.getActionCommand();
        if (s.equals("AddSeries")) {
            setDialogsEnabled(false);
            Cell c = iAceTree.getCurrentCell();
            int time = iImageTime + iTimeInc;
            iDialog = new EIDialog1(iAceTree, this, false, c, time);
        } else if (s.equals("AddOne")) {
        	if (iAceTree.iAddOneDialog != null) return;
            setDialogsEnabled(false);
            Cell c = iAceTree.getCurrentCell();
            int time = iImageTime + iTimeInc;
            //iDialog = new EIDialog2(iAceTree, this, false, c, time);
            iDialog = new AddOneDialog(iAceTree, this, false, c, time);
            iAceTree.iAddOneDialog = (AddOneDialog)iDialog;
        } else {
            if (o == iProperties) {
                new ImageParmsDialog(this);
            } else if (o == iHelp) {
                String item = "/org/rhwlab/help/html/ImageWindowToolbarHelp.html";
                new AceTreeHelp(item);
            }
        }


    }

    public void parentNotifyDialogClosing(JDialog dialog) {
    	iDialog = null;
    	iDialog = iAceTree.iAddOneDialog;
    	println("parentNotifyingDialogClosing, " + iDialog);
        setDialogsEnabled(true);
    }

    public void setDialogsEnabled(boolean enabled) {
        iAddOne.setEnabled(enabled);
        iAddSeries.setEnabled(enabled);
    }

    public void updateCellAnnotation(Cell newCell, String oldName, int time) {
        AnnotInfo ai = isInList(oldName);
        if (ai != null) {
            iAnnotsShown.remove(ai);
            //if (time == -1) time = newCell.getTime();
            Nucleus n = ImageWindow.cNucleiMgr.getNucleusFromHashkey(newCell.getHashKey(), time);
            if (isInList(newCell.getName()) == null) {
                ai = new AnnotInfo(newCell.getName(), n.x, n.y);
                iAnnotsShown.add(ai);
            }
        }
    }

    public void processEditMouseEvent(MouseEvent e) {
    	//println("ImageWindow.processEditMouseEvent, " + e);
        //if (ImageWindow.cEditImage3 == null) return;
        int button = e.getButton();
        //System.out.println("EditImage3: mc: " + button);
        if (iDialog == null) return;
        int type = 0;
        if (iDialog instanceof EIDialog1) type = 1;
        else if (iDialog instanceof AddOneDialog) type = 2;

        println("ImageWindow.processMouseEvent, dialog type = " + type);

        switch(type) {
            case 1:
            // this is an addSeries dialog
            EIDialog1 eid = (EIDialog1)iDialog;
            if (eid.iNothing.isSelected()) {
            	println("Nothing:");
                //eid.processMouseEvent(e);
            }
            if (button == 3) {
                eid.processMouseEvent(e);
                refreshDisplay(null);


            } else if (button == 1) {
                eid.processMouseEvent(e);
            }
            break;

            case 2:
            	// this is an addOne dialog
                AddOneDialog addOne = (AddOneDialog)iDialog;
                addOne.processMouseEvent(e);
                break;
        }
        //super.mouseClicked(e);
    }



}
