/*
 * Created on Apr 28, 2006
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.rhwlab.image;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import org.rhwlab.image.ImageWindow.ColorSchemeUI;

/**
 * @author biowolp
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class ImageParmsDialog extends JDialog implements ActionListener {
    JPanel                          iPanel;
    //SublineageDisplayProperty []    iDispProps;
    ColorSchemeUI []                 iCSUI;
    ImageWindow                     iImgWin;
    JCheckBox                       iAcbTree;
    
    
    public ImageParmsDialog(ImageWindow imgWin) {
        super(imgWin, "Image parameters", false);
        iImgWin = imgWin;
        Border blackline = BorderFactory.createLineBorder(Color.black);
        //iDispProps = imgWin.getDisplayProps();
        iCSUI = new ColorSchemeUI[ImageWindow.iDispProps.length];
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
        JPanel [] testPanel = new JPanel[ImageWindow.iDispProps.length];
        JTextField textField;
        JComboBox cb;
        JPanel labelPanel = new JPanel();
        JLabel sublineage = new JLabel("item");
        JLabel color = new JLabel("color");
        labelPanel.setLayout(new GridLayout(1,2));
        labelPanel.add(sublineage);
        labelPanel.add(color);
        lineagePanel.add(labelPanel);
        
        iAcbTree = new JCheckBox("AcbTree");
        iAcbTree.addActionListener(this);
        ImageWindow.cAcbTree = false;
        dummyPanel.add(iAcbTree);
        dummyPanel.setBorder(blackline);
        
        for (int i=0; i < ImageWindow.iDispProps.length; i++) {
            iCSUI[i] = imgWin.new ColorSchemeUI(i);
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
        setContentPane(iPanel);
        setSize(new Dimension(450, 450));
        setLocationRelativeTo(imgWin);
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        setVisible(true);
        
        
    }

    public void actionPerformed(ActionEvent e) {
        //println("actionPerformed: " + e);
        String command = e.getActionCommand();
        Object o = e.getSource();
        if (o == iAcbTree) {
            boolean b = iAcbTree.isSelected();
            ImageWindow.cAcbTree = b;
            println("ImageParmsDialog.actionPerformed: AcbTree: " + b);
            //println("actionPerformed:2 " + b);
        }
        if (command.equals("Reset")) {
            ImageWindow.iDispProps = iImgWin.getDisplayProps();                
            for (int i=0; i < ImageWindow.iDispProps.length; i++) {
                println("Reset: " + i + CS + ImageWindow.iDispProps[i].iName
                        + CS + ImageWindow.iDispProps[i].iLineageNum);
                iCSUI[i].iLabel.setText(ImageWindow.iDispProps[i].iName);
                iCSUI[i].iCB.setSelectedIndex(ImageWindow.iDispProps[i].iLineageNum);
            }
            
            
        } else if (command.equals("Apply")) {
            println("ImageParmsDialog.actionPerformed: Apply");
            for (int i=0; i < ImageWindow.iDispProps.length; i++) {
                String name = iCSUI[i].iTF.getText();
                if (name.length() == 0) name = "-";
                int num = iCSUI[i].iCB.getSelectedIndex();
                ImageWindow.iDispProps[i].iName = name;
                ImageWindow.iDispProps[i].iLineageNum = num;
            }
            
        }
        //iAceTree.setDispProps3D(iDispProps);
        //updateDisplayedTab();
        
    }
    

    public static void main(String[] args) {
    }
    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ";

}
