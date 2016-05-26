package worm.image3d;

import javax.swing.JPanel;
import java.awt.*;
import java.awt.event.*;

import javax.swing.JLabel;

import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.JButton;

public class ImagePrefs extends JPanel 
implements ActionListener {
	
	public JTextField redBrightness;
	public JTextField time;
	public JButton updateButton;
	private JTextField blurRadius;
	private JTextField greenBrightness;
	private JLabel lblAlpha;
	private JLabel lblBrightness_1;
	private JLabel lblAlpha_1;
	private JTextField txtRedAlpha;
	private JTextField txtGreenAlpha;
	private JButton btnBack;
	private JButton btnForward;
	private JButton btnResetView;

	/**
	 * Create the panel.
	 */
	public ImagePrefs() {
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[]{0, 0, 0};
		gridBagLayout.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		gridBagLayout.columnWeights = new double[]{0.0, 1.0, Double.MIN_VALUE};
		gridBagLayout.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
		setLayout(gridBagLayout);
		
		JLabel lblBrightness = new JLabel("Red");
		GridBagConstraints gbc_lblBrightness = new GridBagConstraints();
		gbc_lblBrightness.anchor = GridBagConstraints.EAST;
		gbc_lblBrightness.insets = new Insets(0, 0, 5, 5);
		gbc_lblBrightness.gridx = 0;
		gbc_lblBrightness.gridy = 0;
		add(lblBrightness, gbc_lblBrightness);
		
		JLabel lblRedBrightness = new JLabel("brightness");
		GridBagConstraints gbc_lblRedBrightness = new GridBagConstraints();
		gbc_lblRedBrightness.anchor = GridBagConstraints.EAST;
		gbc_lblRedBrightness.insets = new Insets(0, 0, 5, 5);
		gbc_lblRedBrightness.gridx = 0;
		gbc_lblRedBrightness.gridy = 1;
		add(lblRedBrightness, gbc_lblRedBrightness);
		
		redBrightness = new JTextField();
		redBrightness.setText("1");
		GridBagConstraints gbc_redBrightness = new GridBagConstraints();
		gbc_redBrightness.insets = new Insets(0, 0, 5, 0);
		gbc_redBrightness.fill = GridBagConstraints.HORIZONTAL;
		gbc_redBrightness.gridx = 1;
		gbc_redBrightness.gridy = 1;
		add(redBrightness, gbc_redBrightness);
		redBrightness.setColumns(4);
		
		lblAlpha = new JLabel("alpha");
		GridBagConstraints gbc_lblAlpha = new GridBagConstraints();
		gbc_lblAlpha.anchor = GridBagConstraints.EAST;
		gbc_lblAlpha.insets = new Insets(0, 0, 5, 5);
		gbc_lblAlpha.gridx = 0;
		gbc_lblAlpha.gridy = 2;
		add(lblAlpha, gbc_lblAlpha);
		
		txtRedAlpha = new JTextField();
		txtRedAlpha.setText("1");
		GridBagConstraints gbc_txtRedAlpha = new GridBagConstraints();
		gbc_txtRedAlpha.insets = new Insets(0, 0, 5, 0);
		gbc_txtRedAlpha.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtRedAlpha.gridx = 1;
		gbc_txtRedAlpha.gridy = 2;
		add(txtRedAlpha, gbc_txtRedAlpha);
		txtRedAlpha.setColumns(10);
		
		JLabel lblGreen = new JLabel("Green");
		GridBagConstraints gbc_lblGreen = new GridBagConstraints();
		gbc_lblGreen.anchor = GridBagConstraints.EAST;
		gbc_lblGreen.insets = new Insets(0, 0, 5, 5);
		gbc_lblGreen.gridx = 0;
		gbc_lblGreen.gridy = 3;
		add(lblGreen, gbc_lblGreen);
		
		lblBrightness_1 = new JLabel("brightness");
		GridBagConstraints gbc_lblBrightness_1 = new GridBagConstraints();
		gbc_lblBrightness_1.insets = new Insets(0, 0, 5, 5);
		gbc_lblBrightness_1.anchor = GridBagConstraints.EAST;
		gbc_lblBrightness_1.gridx = 0;
		gbc_lblBrightness_1.gridy = 4;
		add(lblBrightness_1, gbc_lblBrightness_1);
		
		greenBrightness = new JTextField();
		greenBrightness.setText("1");
		GridBagConstraints gbc_greenBrightness = new GridBagConstraints();
		gbc_greenBrightness.insets = new Insets(0, 0, 5, 0);
		gbc_greenBrightness.fill = GridBagConstraints.HORIZONTAL;
		gbc_greenBrightness.gridx = 1;
		gbc_greenBrightness.gridy = 4;
		add(greenBrightness, gbc_greenBrightness);
		greenBrightness.setColumns(10);
		
		lblAlpha_1 = new JLabel("alpha");
		GridBagConstraints gbc_lblAlpha_1 = new GridBagConstraints();
		gbc_lblAlpha_1.anchor = GridBagConstraints.EAST;
		gbc_lblAlpha_1.insets = new Insets(0, 0, 5, 5);
		gbc_lblAlpha_1.gridx = 0;
		gbc_lblAlpha_1.gridy = 5;
		add(lblAlpha_1, gbc_lblAlpha_1);
		
		txtGreenAlpha = new JTextField();
		txtGreenAlpha.setText("1");
		GridBagConstraints gbc_txtGreenAlpha = new GridBagConstraints();
		gbc_txtGreenAlpha.insets = new Insets(0, 0, 5, 0);
		gbc_txtGreenAlpha.fill = GridBagConstraints.HORIZONTAL;
		gbc_txtGreenAlpha.gridx = 1;
		gbc_txtGreenAlpha.gridy = 5;
		add(txtGreenAlpha, gbc_txtGreenAlpha);
		txtGreenAlpha.setColumns(10);
		
		JLabel lblBlurRadius = new JLabel("Blur radius");
		GridBagConstraints gbc_lblBlurRadius = new GridBagConstraints();
		gbc_lblBlurRadius.anchor = GridBagConstraints.EAST;
		gbc_lblBlurRadius.insets = new Insets(0, 0, 5, 5);
		gbc_lblBlurRadius.gridx = 0;
		gbc_lblBlurRadius.gridy = 7;
		add(lblBlurRadius, gbc_lblBlurRadius);
		
		blurRadius = new JTextField();
		blurRadius.setText("2");
		GridBagConstraints gbc_blurRadius = new GridBagConstraints();
		gbc_blurRadius.insets = new Insets(0, 0, 5, 0);
		gbc_blurRadius.fill = GridBagConstraints.HORIZONTAL;
		gbc_blurRadius.gridx = 1;
		gbc_blurRadius.gridy = 7;
		add(blurRadius, gbc_blurRadius);
		blurRadius.setColumns(10);
		
		JLabel lblTime = new JLabel("Time");
		GridBagConstraints gbc_lblTime = new GridBagConstraints();
		gbc_lblTime.anchor = GridBagConstraints.EAST;
		gbc_lblTime.insets = new Insets(0, 0, 5, 5);
		gbc_lblTime.gridx = 0;
		gbc_lblTime.gridy = 9;
		add(lblTime, gbc_lblTime);
		
		time = new JTextField();
		time.setText("55");
		GridBagConstraints gbc_time = new GridBagConstraints();
		gbc_time.insets = new Insets(0, 0, 5, 0);
		gbc_time.fill = GridBagConstraints.HORIZONTAL;
		gbc_time.gridx = 1;
		gbc_time.gridy = 9;
		add(time, gbc_time);
		time.setColumns(4);
		
		btnBack = new JButton("Back");
		GridBagConstraints gbc_btnBack = new GridBagConstraints();
		gbc_btnBack.insets = new Insets(0, 0, 5, 5);
		gbc_btnBack.gridx = 0;
		gbc_btnBack.gridy = 10;
		add(btnBack, gbc_btnBack);
		btnBack.addActionListener(this);
		
		btnForward = new JButton("Forward");
		GridBagConstraints gbc_btnForward = new GridBagConstraints();
		gbc_btnForward.insets = new Insets(0, 0, 5, 0);
		gbc_btnForward.gridx = 1;
		gbc_btnForward.gridy = 10;
		add(btnForward, gbc_btnForward);
		gbc_btnBack.insets = new Insets(0, 0, 5, 0);
		gbc_btnBack.gridx = 1;
		gbc_btnBack.gridy = 12;
		btnForward.addActionListener(this);
		
		btnResetView = new JButton("Reset view");
		GridBagConstraints gbc_btnResetView = new GridBagConstraints();
		gbc_btnResetView.insets = new Insets(0, 0, 5, 0);
		gbc_btnResetView.gridx = 1;
		gbc_btnResetView.gridy = 12;
		add(btnResetView, gbc_btnResetView);
		
		updateButton = new JButton("Update");
		GridBagConstraints gbc_btnUpdate = new GridBagConstraints();
		gbc_btnUpdate.gridx = 1;
		gbc_btnUpdate.gridy = 14;
		add(updateButton, gbc_btnUpdate);
	}

	/** Utility to return the current actual numbers. */
	public class Prefs {
	
		public float redBrightness, redAlpha, greenBrightness, greenAlpha;
		
		int time;
		
		public Prefs() {
		}
	}
	
	public ImagePrefs.Prefs getPrefs() {
		ImagePrefs.Prefs p = new ImagePrefs.Prefs();
		
		p.redBrightness = new Float(redBrightness.getText()).floatValue();
		p.redAlpha = new Float(txtRedAlpha.getText()).floatValue();
		p.greenBrightness = new Float(greenBrightness.getText()).floatValue();
		p.greenAlpha = new Float(txtGreenAlpha.getText()).floatValue();
		
		p.time = new Integer(time.getText()).intValue();
		
		return p;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		
		if (e.getActionCommand().equals("Back")) {
			int t = new Integer(time.getText()).intValue();
			time.setText(new Integer(t-1).toString());
			updateButton.doClick();
		}
		
		if (e.getActionCommand().equals("Forward")) {
			int t = new Integer(time.getText()).intValue();
			time.setText(new Integer(t+1).toString());
			updateButton.doClick();
		}
	}
}
