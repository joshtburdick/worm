package edu.murrlab.image;

import java.awt.BorderLayout;
import java.awt.FlowLayout;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import com.jgoodies.forms.layout.FormLayout;
import com.jgoodies.forms.layout.ColumnSpec;
import com.jgoodies.forms.layout.RowSpec;
import java.awt.GridLayout;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import java.awt.event.InputMethodListener;
import java.awt.event.InputMethodEvent;

import org.rhwlab.tree.Cell;
import javax.swing.JRadioButton;
import javax.swing.Box;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.Component;
import javax.swing.ButtonGroup;

public class SmoothDialog extends JDialog {

	private Cell root;
	
	private final JPanel contentPanel = new JPanel();
	private JTextField numPoints;
	private JButton okButton;
	private final ButtonGroup method = new ButtonGroup();
	private JRadioButton gaussian;

	/**
	 * Create the dialog.
	 */
	public SmoothDialog(Cell root) {
		this.root = root;
		setModal(true);
		setTitle("Median smooth");
		setBounds(100, 100, 450, 300);
		getContentPane().setLayout(new BorderLayout());
		contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		getContentPane().add(contentPanel, BorderLayout.CENTER);
		GridBagLayout gbl_contentPanel = new GridBagLayout();
		gbl_contentPanel.columnWidths = new int[]{368, 0};
		gbl_contentPanel.rowHeights = new int[]{95, 42, 0};
		gbl_contentPanel.columnWeights = new double[]{0.0, Double.MIN_VALUE};
		gbl_contentPanel.rowWeights = new double[]{0.0, 0.0, Double.MIN_VALUE};
		contentPanel.setLayout(gbl_contentPanel);
		{
			Box verticalBox = Box.createVerticalBox();
			GridBagConstraints gbc_verticalBox = new GridBagConstraints();
			gbc_verticalBox.anchor = GridBagConstraints.WEST;
			gbc_verticalBox.gridx = 0;
			gbc_verticalBox.gridy = 1;
			contentPanel.add(verticalBox, gbc_verticalBox);
			{
				Box horizontalBox = Box.createHorizontalBox();
				verticalBox.add(horizontalBox);
			}
			{
				Box horizontalBox_1 = Box.createHorizontalBox();
				verticalBox.add(horizontalBox_1);
				{
					JLabel lblNewLabel = new JLabel("Smoothing method");
					horizontalBox_1.add(lblNewLabel);
				}
				{
					Box verticalBox_1 = Box.createVerticalBox();
					horizontalBox_1.add(verticalBox_1);
					{
						gaussian = new JRadioButton("Gaussian convolve");
						method.add(gaussian);
						gaussian.setSelected(true);
						verticalBox_1.add(gaussian);
					}
					{
						JRadioButton median = new JRadioButton("median of window");
						method.add(median);
						verticalBox_1.add(median);
					}
				}
			}
			{
				Box horizontalBox = Box.createHorizontalBox();
				verticalBox.add(horizontalBox);
				{
					JLabel lblNumberOfTime = new JLabel("Window size (in time points)");
					horizontalBox.add(lblNumberOfTime);
					lblNumberOfTime.setHorizontalAlignment(SwingConstants.LEFT);
				}
				{
					numPoints = new JTextField();
					horizontalBox.add(numPoints);

					numPoints.addActionListener(new ActionListener() {
						public void actionPerformed(ActionEvent e) {
							okButton.doClick();
						}
					});
					numPoints.setHorizontalAlignment(SwingConstants.RIGHT);
					numPoints.setText("4");
					numPoints.setColumns(5);
				}
			}
			{
				Component verticalGlue = Box.createVerticalGlue();
				verticalBox.add(verticalGlue);
			}
			{
				Box horizontalBox = Box.createHorizontalBox();
				verticalBox.add(horizontalBox);
			}
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						smooth();
					}
				});
				okButton.setActionCommand("OK");
				buttonPane.add(okButton);
				getRootPane().setDefaultButton(okButton);
			}
			{
				JButton cancelButton = new JButton("Cancel");
				cancelButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						dispose();
					}
				});
				cancelButton.setActionCommand("Cancel");
				buttonPane.add(cancelButton);
			}
		}
	}

	public void smooth() {
		try {
			int k = Integer.parseInt(numPoints.getText());
			if (k >= 1 && k < 500) {
				// MedianSmoother s = new MedianSmoother(root, k);
				
				boolean m = getGaussian().isSelected();
				System.out.println(m);
				Smoother s = new Smoother(root, getGaussian().isSelected(), k);
				System.out.println("computing medians...");
				s.computeMedians();
				System.out.println("done computing medians...");
				dispose();
			}
			else {
				numPoints.setText("4");
			}
		}
		catch (NumberFormatException e) {
			numPoints.setText("4");
		}
	}
	protected JRadioButton getGaussian() {
		return gaussian;
	}
}
