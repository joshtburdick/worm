/*
 * Created on Apr 28, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.rhwlab.acetree;


import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JToolBar;

/**
 * @author biowolp
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class PlayerControl extends JPanel implements ActionListener, Runnable {

    AceTree iAceTree;
    JToolBar iToolBar;
    JButton iPlay;
    JButton iPause;
    JButton iReverse;
    JButton iStepForward;
    JButton iStepBack;
    boolean iRunning;
    boolean iForward;
    int iDwell;

    public PlayerControl(AceTree aceTree) {
        super(new BorderLayout());
        iAceTree = aceTree;
        //Create the toolbar.
        iToolBar = new JToolBar("");
        addButtons();
        setPreferredSize(new Dimension(100, 30));
        add(iToolBar, BorderLayout.PAGE_START);
        iRunning = false;
        iForward = true;
        iDwell = 100;
    }

    private void addButtons() {
        iStepBack = makeButton("/images/StepBack16");
        iToolBar.add(iStepBack);
        iReverse = makeButton("/images/PlayBack16");
        iToolBar.add(iReverse);
        iToolBar.add(new JToolBar.Separator());
        iPause = makeButton("/images/Pause16");
        iToolBar.add(iPause);
        iToolBar.add(new JToolBar.Separator());
        iPlay = makeButton("/images/Play16");
        iToolBar.add(iPlay);
        iStepForward = makeButton("/images/StepForward16");
        iToolBar.add(iStepForward);
        setEnabledAll(true);
        iPause.setEnabled(false);
    }

    private void setEnabledAll(boolean enabled) {
        int i = 0;
        Component c;
        //boolean more = true;
        while(true) {
            c = iToolBar.getComponentAtIndex(i++);
            if (c == null) break;
            c.setEnabled(enabled);
        }
    }

    private JButton makeButton(String imageName) {
        JButton b = new JButton();
        String imgLoc = imageName + ".gif";
        URL imageURL = PlayerControl.class.getResource(imgLoc);
        b.setIcon(new ImageIcon(imageURL, "x"));
        b.addActionListener(this);
        return b;
    }

    public void run() {
        boolean b; // enables run to exit when movie hits the wall
        while (iRunning) {
            if (iForward) b = iAceTree.nextImage();
            else b = iAceTree.prevImage();
            if (b) {
                try {
                    Thread.sleep(iDwell);
                } catch(InterruptedException ie) {
                     ie.printStackTrace();
                 }
            } else {
                iRunning = false;
            }
        }
        setEnabledAll(true);
        iPause.setEnabled(false);
    }


    public void stop() {
    	iRunning = false;
    }

    public void pause() {
        iRunning = false;
        setEnabledAll(true);
        iPause.setEnabled(false);

    }

    public void actionPerformed(ActionEvent e) {
        Object o = e.getSource();
        if (o == iPause) {
            //iEventPusher.stop();
            pause();
        }
        if (o == iPlay) {
            //iEventPusher = new EventPusher(iAceTree, 30);
            //iEventPusher.start(true);
            iRunning = true;
            iForward = true;
            setEnabledAll(false);
            iPause.setEnabled(true);
            new Thread(this, "TEST").start();

        } else if (o == iReverse) {
            //iEventPusher = new EventPusher(iAceTree, 30);
            //iEventPusher.start(false);
            if (iRunning) return;
            iRunning = true;
            iForward = false;
            setEnabledAll(false);
            iPause.setEnabled(true);
            new Thread(this, "TEST").start();

        } else if (o == iStepForward) {
            if (iRunning) return;
            iAceTree.nextImage();
        } else if (o == iStepBack) {
            if (iRunning) return;
            iAceTree.prevImage();
        }


    }

    public static void main(String[] args) {
    }
}
