package worm.tree;

import javax.swing.*;

import java.awt.event.*;
import java.awt.print.*;
import java.awt.Color;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;
import java.awt.Component;
import java.awt.Dimension;
import java.io.*;
import java.util.HashMap;

import javax.print.*;
import javax.print.attribute.*;
import javax.print.attribute.standard.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;

import slider.RangeSlider;

/** Preferences for the TreeExprViewer. */
public class TreeExprViewerPrefs extends JPanel {

	private final TreePanel tp;
	private JTextField root;
	private JTextField endTime;
	private JCheckBox chckbxShownames;
	private JButton btnPrint;
	private JTextField txtImageFile;
	private Component glue;
    private JLabel lblImageFile;
    private JTextField txtImageFile1;
    private JButton btnPickRed;
    private JButton btnPickGreen;
    private JCheckBox chckbxShowMaxIntensity;
    private JLabel lblMinIntensity;
    private JTextField txtMinintensity;
    private JLabel lblMaxIntensity;
    private JTextField txtMaxintensity;
    private JTextField txtMinintensity_1;
    private JTextField txtMaxintensity_1;
    private JButton btnPickBlue;
    private JTextField txtMinIntensity_2;
    private JTextField txtMaxIntensity_2;
    private JTextField txtImageFile2;
    private JCheckBox chckbxShowCellFates;

	public TreeExprViewerPrefs(TreePanel tp) {
		this.tp = tp;
		
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[]{0, 0, 0, 0, 0, 0, 38, 0, 0};
		gridBagLayout.columnWeights = new double[]{0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0};
		gridBagLayout.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		setLayout(gridBagLayout);
		
		lblImageFile = new JLabel("image file");
		GridBagConstraints gbc_lblImageFile = new GridBagConstraints();
		gbc_lblImageFile.insets = new Insets(0, 0, 5, 5);
		gbc_lblImageFile.gridx = 4;
		gbc_lblImageFile.gridy = 0;
		add(lblImageFile, gbc_lblImageFile);
	
		JLabel lblRoot = new JLabel("root");
		GridBagConstraints gbc_lblRoot = new GridBagConstraints();
		gbc_lblRoot.anchor = GridBagConstraints.EAST;
		gbc_lblRoot.insets = new Insets(0, 0, 5, 5);
		gbc_lblRoot.gridx = 0;
		gbc_lblRoot.gridy = 0;
		add(lblRoot, gbc_lblRoot);
		
		root = new JTextField();
		root.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateTreePanel();
			}
		});
		root.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				updateTreePanel();
			}
		});
		root.setToolTipText("Root of the tree to display.");
		root.setText("P0");
		GridBagConstraints gbc_root = new GridBagConstraints();
		gbc_root.insets = new Insets(0, 0, 5, 5);
		gbc_root.fill = GridBagConstraints.HORIZONTAL;
		gbc_root.gridx = 1;
		gbc_root.gridy = 0;
		add(root, gbc_root);
		root.setColumns(10);
		
		lblMinIntensity = new JLabel("min. intensity");
		GridBagConstraints gbc_lblMinIntensity = new GridBagConstraints();
		gbc_lblMinIntensity.insets = new Insets(0, 0, 5, 5);
		gbc_lblMinIntensity.gridx = 5;
		gbc_lblMinIntensity.gridy = 0;
		add(lblMinIntensity, gbc_lblMinIntensity);
		
		lblMaxIntensity = new JLabel("max. intensity");
		GridBagConstraints gbc_lblMaxIntensity = new GridBagConstraints();
		gbc_lblMaxIntensity.insets = new Insets(0, 0, 5, 5);
		gbc_lblMaxIntensity.gridx = 6;
		gbc_lblMaxIntensity.gridy = 0;
		add(lblMaxIntensity, gbc_lblMaxIntensity);
		
		glue = Box.createGlue();
		glue.setMinimumSize(new Dimension(10, 10));
		GridBagConstraints gbc_glue = new GridBagConstraints();
		gbc_glue.insets = new Insets(0, 0, 5, 0);
		gbc_glue.gridx = 8;
		gbc_glue.gridy = 0;
		add(glue, gbc_glue);
		
		JLabel lblEndTime = new JLabel("end time");
		GridBagConstraints gbc_lblEndTime = new GridBagConstraints();
		gbc_lblEndTime.anchor = GridBagConstraints.EAST;
		gbc_lblEndTime.insets = new Insets(0, 0, 5, 5);
		gbc_lblEndTime.gridx = 0;
		gbc_lblEndTime.gridy = 1;
		add(lblEndTime, gbc_lblEndTime);
		
		endTime = new JTextField();
		endTime.setMaximumSize(new Dimension(60, 2147483647));
		endTime.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateTreePanel();
			}
		});
		endTime.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				updateTreePanel();
			}
		});
		endTime.setText("100");
		GridBagConstraints gbc_endTime = new GridBagConstraints();
		gbc_endTime.insets = new Insets(0, 0, 5, 5);
		gbc_endTime.fill = GridBagConstraints.HORIZONTAL;
		gbc_endTime.gridx = 1;
		gbc_endTime.gridy = 1;
		add(endTime, gbc_endTime);
		endTime.setColumns(10);
		
		btnPickRed = new JButton("Pick red");
		btnPickRed.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				loadImageData(0);
			}
		});
		GridBagConstraints gbc_btnPickRed = new GridBagConstraints();
		gbc_btnPickRed.insets = new Insets(0, 0, 5, 5);
		gbc_btnPickRed.gridx = 3;
		gbc_btnPickRed.gridy = 1;
		add(btnPickRed, gbc_btnPickRed);
		
		txtImageFile = new JTextField();
		txtImageFile.setEditable(false);
		GridBagConstraints gbc_txtImageFile = new GridBagConstraints();
		gbc_txtImageFile.insets = new Insets(0, 0, 5, 5);
		gbc_txtImageFile.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtImageFile.gridx = 4;
		gbc_txtImageFile.gridy = 1;
		add(txtImageFile, gbc_txtImageFile);
		txtImageFile.setColumns(10);
		
		txtMinintensity = new JTextField();
		txtMinintensity.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateTreePanel();
			}
		});
		txtMinintensity.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				updateTreePanel();
			}
		});
		txtMinintensity.setText("0");
		GridBagConstraints gbc_txtMinintensity = new GridBagConstraints();
		gbc_txtMinintensity.insets = new Insets(0, 0, 5, 5);
		gbc_txtMinintensity.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtMinintensity.gridx = 5;
		gbc_txtMinintensity.gridy = 1;
		add(txtMinintensity, gbc_txtMinintensity);
		txtMinintensity.setColumns(10);
		
		txtMaxintensity = new JTextField();
		txtMaxintensity.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateTreePanel();
			}
		});
		txtMaxintensity.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				updateTreePanel();
			}
		});
		txtMaxintensity.setText("1000");
		GridBagConstraints gbc_txtMaxintensity = new GridBagConstraints();
		gbc_txtMaxintensity.insets = new Insets(0, 0, 5, 5);
		gbc_txtMaxintensity.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtMaxintensity.gridx = 6;
		gbc_txtMaxintensity.gridy = 1;
		add(txtMaxintensity, gbc_txtMaxintensity);
		txtMaxintensity.setColumns(10);
		
		btnPickGreen = new JButton("Pick green");
		btnPickGreen.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				loadImageData(1);
			}
		});
		GridBagConstraints gbc_btnPickGreen = new GridBagConstraints();
		gbc_btnPickGreen.insets = new Insets(0, 0, 5, 5);
		gbc_btnPickGreen.gridx = 3;
		gbc_btnPickGreen.gridy = 2;
		add(btnPickGreen, gbc_btnPickGreen);
		
		txtImageFile1 = new JTextField();
		txtImageFile1.setEditable(false);
		txtImageFile1.setColumns(10);
		GridBagConstraints gbc_txtImageFile1 = new GridBagConstraints();
		gbc_txtImageFile1.insets = new Insets(0, 0, 5, 5);
		gbc_txtImageFile1.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtImageFile1.gridx = 4;
		gbc_txtImageFile1.gridy = 2;
		add(txtImageFile1, gbc_txtImageFile1);
		
		txtMinintensity_1 = new JTextField();
		txtMinintensity_1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateTreePanel();
			}
		});
		txtMinintensity_1.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				updateTreePanel();
			}
		});
		txtMinintensity_1.setText("0");
		GridBagConstraints gbc_txtMinintensity_1 = new GridBagConstraints();
		gbc_txtMinintensity_1.insets = new Insets(0, 0, 5, 5);
		gbc_txtMinintensity_1.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtMinintensity_1.gridx = 5;
		gbc_txtMinintensity_1.gridy = 2;
		add(txtMinintensity_1, gbc_txtMinintensity_1);
		txtMinintensity_1.setColumns(10);
		
		txtMaxintensity_1 = new JTextField();
		txtMaxintensity_1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateTreePanel();
			}
		});
		txtMaxintensity_1.addFocusListener(new FocusAdapter() {
			@Override
			public void focusLost(FocusEvent e) {
				updateTreePanel();
			}
		});
		txtMaxintensity_1.setText("1000");
		GridBagConstraints gbc_txtMaxintensity_1 = new GridBagConstraints();
		gbc_txtMaxintensity_1.insets = new Insets(0, 0, 5, 5);
		gbc_txtMaxintensity_1.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtMaxintensity_1.gridx = 6;
		gbc_txtMaxintensity_1.gridy = 2;
		add(txtMaxintensity_1, gbc_txtMaxintensity_1);
		txtMaxintensity_1.setColumns(10);
		
		chckbxShowMaxIntensity = new JCheckBox("show max. intensity");
		chckbxShowMaxIntensity.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				updateTreePanel();
			}
		});
		
		btnPickBlue = new JButton("Pick blue");
		btnPickBlue.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				loadImageData(2);
			}
		});
		GridBagConstraints gbc_btnPickBlue = new GridBagConstraints();
		gbc_btnPickBlue.insets = new Insets(0, 0, 5, 5);
		gbc_btnPickBlue.gridx = 3;
		gbc_btnPickBlue.gridy = 3;
		add(btnPickBlue, gbc_btnPickBlue);
		
		txtImageFile2 = new JTextField();
		txtImageFile2.setEditable(false);
		GridBagConstraints gbc_txtImageFile2 = new GridBagConstraints();
		gbc_txtImageFile2.insets = new Insets(0, 0, 5, 5);
		gbc_txtImageFile2.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtImageFile2.gridx = 4;
		gbc_txtImageFile2.gridy = 3;
		add(txtImageFile2, gbc_txtImageFile2);
		txtImageFile2.setColumns(10);
		
		txtMinIntensity_2 = new JTextField();
		txtMinIntensity_2.setText("0");
		GridBagConstraints gbc_txtMinIntensity_2 = new GridBagConstraints();
		gbc_txtMinIntensity_2.insets = new Insets(0, 0, 5, 5);
		gbc_txtMinIntensity_2.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtMinIntensity_2.gridx = 5;
		gbc_txtMinIntensity_2.gridy = 3;
		add(txtMinIntensity_2, gbc_txtMinIntensity_2);
		txtMinIntensity_2.setColumns(10);
		
		txtMaxIntensity_2 = new JTextField();
		txtMaxIntensity_2.setText("1000");
		GridBagConstraints gbc_txtMaxIntensity_2 = new GridBagConstraints();
		gbc_txtMaxIntensity_2.insets = new Insets(0, 0, 5, 5);
		gbc_txtMaxIntensity_2.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtMaxIntensity_2.gridx = 6;
		gbc_txtMaxIntensity_2.gridy = 3;
		add(txtMaxIntensity_2, gbc_txtMaxIntensity_2);
		txtMaxIntensity_2.setColumns(10);
		chckbxShowMaxIntensity.setSelected(true);
		GridBagConstraints gbc_chckbxShowMaxIntensity = new GridBagConstraints();
		gbc_chckbxShowMaxIntensity.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxShowMaxIntensity.gridx = 1;
		gbc_chckbxShowMaxIntensity.gridy = 4;
		add(chckbxShowMaxIntensity, gbc_chckbxShowMaxIntensity);
		
		btnPrint = new JButton("Print");
		btnPrint.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				print();
			}
		});
		
		chckbxShownames = new JCheckBox("show lineage names");
		chckbxShownames.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				updateTreePanel();
			}
		});
		chckbxShownames.setSelected(true);
		GridBagConstraints gbc_chckbxShownames = new GridBagConstraints();
		gbc_chckbxShownames.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxShownames.gridx = 2;
		gbc_chckbxShownames.gridy = 4;
		add(chckbxShownames, gbc_chckbxShownames);
		
		chckbxShowCellFates = new JCheckBox("show cell fates");
		chckbxShowCellFates.setEnabled(false);
		GridBagConstraints gbc_chckbxShowCellFates = new GridBagConstraints();
		gbc_chckbxShowCellFates.insets = new Insets(0, 0, 0, 5);
		gbc_chckbxShowCellFates.gridx = 2;
		gbc_chckbxShowCellFates.gridy = 5;
		add(chckbxShowCellFates, gbc_chckbxShowCellFates);
		GridBagConstraints gbc_btnPrint = new GridBagConstraints();
		gbc_btnPrint.insets = new Insets(0, 0, 0, 5);
		gbc_btnPrint.gridx = 6;
		gbc_btnPrint.gridy = 5;
		add(btnPrint, gbc_btnPrint);
	}

	private void updateTreePanel() {
		
		// XXX hack
		if (tp.signal.size() == 0)
			return;
		
//		System.out.println("focus lost");
		tp.root = root.getText();
		tp.lastT = Integer.parseInt(endTime.getText());
//		tp.scale = Float.parseFloat(scale.getText());
		tp.showMaxIntensitySoFar = chckbxShowMaxIntensity.isSelected();
		tp.showNames = chckbxShownames.isSelected();

		if (tp.signal.size() >= 1) {
			tp.signal.elementAt(0).loIntensity = Float.parseFloat(txtMinintensity.getText());
			tp.signal.elementAt(0).hiIntensity = Float.parseFloat(txtMaxintensity.getText());
		}
		if (tp.signal.size() >= 2) {
			tp.signal.elementAt(1).loIntensity = Float.parseFloat(txtMinintensity_1.getText());
			tp.signal.elementAt(1).hiIntensity = Float.parseFloat(txtMaxintensity_1.getText());
		}
		if (tp.signal.size() >= 3) {
			tp.signal.elementAt(2).loIntensity = Float.parseFloat(txtMinIntensity_2.getText());
			tp.signal.elementAt(2).hiIntensity = Float.parseFloat(txtMaxIntensity_2.getText());
		}

		tp.repaint();
	}
	
	public JCheckBox getChckbxShownames() {
		return chckbxShownames;
	}

	/** Loads one channel, based on the series name.
	 *  
	 * @param channel  index of the channel to load
	 */
	private void loadFromSeriesName(int channel) {
	
		
	
		
	}
	
	/** Loads one channel's worth of data,
	 * 
	 * @param channel  index of the channel to load
	 */
	private void loadImageData(int channel) {
	
		// show file chooser to choose SCD file 
    	JFileChooser fileChooser = new JFileChooser();
 	
        // ??? what should this be set to?
    	fileChooser.setCurrentDirectory(new File("."));
        fileChooser.setSelectedFile(new File(""));
        
        // check whether we should import the file
        int returnVal = fileChooser.showOpenDialog(this);
        if (returnVal != JFileChooser.APPROVE_OPTION)
        	return;
        
        // make sure it exists
        File file = fileChooser.getSelectedFile();       
        if (!file.exists()) {
            JOptionPane.showMessageDialog(this, "File " + file.toString() + " doesn't exist");
        	return;
        }
		
		// read file (or try to)
		try {
			CDFileParser p = new CDFileParser(new FileReader(file));
			HashMap<String, Cell> h = p.parseFile();
			if (h == null)
				throw new Exception("failed to read file");   // XXX
			TreePanelSignal si = new TreePanelSignal(h);
			switch (channel) {
				case 0: si.color = Color.red; break;
				case 1: si.color = Color.green; break;
				case 2: si.color = Color.blue; break;
			}
			System.out.println("read file " + file);
			
			HashMap<String, Cell> expr = tp.expr;
			// if there isn't an expression dataset already, use this
			if (expr == null)
				expr = h;

			if (tp.signal.size() < channel+1)
				tp.signal.setSize(channel+1);
			tp.signal.setElementAt(si, channel);
			
			// add in (smoothed) signal			
			GaussianSmoother gs =
					new GaussianSmoother(tp.signal.elementAt(channel).cell.get("P0"));
//			gs.smoothMaximizing(expr, channel);  // FIXME enable one of these
			gs.smoothEWMA(expr, channel);
			
			// finally, redraw
			tp.expr = expr;
			updateTreePanel();
		}
		catch (Exception e) {
			e.printStackTrace();
            JOptionPane.showMessageDialog(this, "Failed to read file " + file.toString());
        	return;
		}
		
		// if we reached this point, update the relevant label
		switch (channel) {
		case 0: txtImageFile.setText(file.getName()); break;
		case 1: txtImageFile1.setText(file.getName()); break;
		case 2: txtImageFile2.setText(file.getName()); break;
		}
	}
	
	/** Prints the tree. */
	public void print() {
		
		// the following somewhat cryptic code is from
		// http://docs.oracle.com/javase/7/docs/technotes/guides/jps/spec/appendix_2DtoStream.fm.html#7083
		// FIXME: change background from grey
		
		/* Use the pre-defined flavor for a Printable from an InputStream */
        DocFlavor flavor = DocFlavor.SERVICE_FORMATTED.PRINTABLE;

        /* Specify the type of the output stream */
        String psMimeType = DocFlavor.BYTE_ARRAY.POSTSCRIPT.getMimeType();

        /* Locate factory which can export a GIF image stream as Postscript */
        StreamPrintServiceFactory[] factories =
        StreamPrintServiceFactory.lookupStreamPrintServiceFactories(
                                flavor, psMimeType);
        assert( factories.length > 0 );

        try {
            /* Create a file for the exported postscript */
        	JFileChooser fileChooser = new JFileChooser();
            // ??? what should this be set to?
        	fileChooser.setCurrentDirectory(new File("."));
            fileChooser.setSelectedFile(new File(""));
            
            // select output file
            int returnVal = fileChooser.showSaveDialog(this);
            if (returnVal != JFileChooser.APPROVE_OPTION)
            	return;
            
            // check if it exists
            File file = fileChooser.getSelectedFile();        
            if (file.exists()) {
            	int n = JOptionPane.showConfirmDialog(
            	    this,
            	    "File " + file + " exists. Overwrite?",
            	    "",
            	    JOptionPane.YES_NO_OPTION);
            	if (n == 1)
            		return;
            }
            
            FileOutputStream fos = new FileOutputStream(file);

            /* Create a Stream printer for Postscript */
            StreamPrintService sps = factories[0].getPrintService(fos);

            /* Create and call a Print Job */
            DocPrintJob pj = sps.createPrintJob();
 //           PrintRequestAttributeSet attrSet = new HashPrintRequestAttributeSet();
 //           attrSet.add(MediaSizeName.NA_LETTER);
 //           attrSet.add(OrientationRequested.LANDSCAPE);
//            attrSet.add(new MediaSize(3f, 11f, MediaSize.INCH));
//             attrSet.add(new MediaPrintableArea(0.25f, 0.25f, 3.25f, 7f, MediaPrintableArea.INCH));                     
            DocAttributeSet das = new HashDocAttributeSet();
            das.add(MediaSizeName.NA_LEGAL);
            das.add(OrientationRequested.LANDSCAPE);

            Doc doc = new SimpleDoc(tp, flavor, das);
            pj.print(doc, null);
            fos.close();
            
        } catch (Exception e) {
        	e.printStackTrace();
        }

	}
	protected JTextField getTxtImagefile() {
		return txtImageFile;
	}
	protected JTextField getTxtImageFile1() {
		return txtImageFile1;
	}
	protected JTextField getTxtImageFile2() {
		return txtImageFile2;
	}
}
