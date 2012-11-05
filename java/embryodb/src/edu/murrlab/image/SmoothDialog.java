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

public class SmoothDialog extends JDialog {

	private Cell root;
	
	private final JPanel contentPanel = new JPanel();
	private JTextField numPoints;
	private JButton okButton;

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
		contentPanel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 5));
		{
			JLabel lblNumberOfTime = new JLabel("Number of time points to smooth");
			lblNumberOfTime.setHorizontalAlignment(SwingConstants.LEFT);
			contentPanel.add(lblNumberOfTime);
		}
		{
			numPoints = new JTextField();

			numPoints.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					okButton.doClick();
				}
			});
			numPoints.setHorizontalAlignment(SwingConstants.RIGHT);
			numPoints.setText("30");
			contentPanel.add(numPoints);
			numPoints.setColumns(5);
		}
		{
			JPanel buttonPane = new JPanel();
			buttonPane.setLayout(new FlowLayout(FlowLayout.RIGHT));
			getContentPane().add(buttonPane, BorderLayout.SOUTH);
			{
				okButton = new JButton("OK");
				okButton.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						medianSmooth();
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

	public void medianSmooth() {
		try {
			int k = Integer.parseInt(numPoints.getText());
			if (k >= 1 && k < 500) {
				MedianSmoother s = new MedianSmoother(root, k);
				System.out.println("computing medians...");
				s.computeMedians();
				System.out.println("done computing medians...");
				dispose();
			}
			else {
				numPoints.setText("30");
			}
		}
		catch (NumberFormatException e) {
			numPoints.setText("30");
		}
	}
}
