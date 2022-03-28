package org.hepforge.alohep.gui;

import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.Map.Entry;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.Border;

import org.hepforge.alohep.AloHEP;
import org.hepforge.alohep.database.AcceleratorData;
import org.hepforge.alohep.database.ParticleData;

public class ParticlePanel extends JPanel{

	
	private static final long serialVersionUID = 1L;

	private AloHEP alohep;
	
	private HashMap<String, VariablePanel> particleVarPanels;
	private HashMap<String, JRadioButton> accRadioButtons;
	private JPanel comboboxPanel;
	private JComboBox<String> parTypeCbox;
	private JPanel parVarMainPanel;
	private JPanel accTypeMainPanel;

	public ParticlePanel(AloHEP alohep)
	{
		this.alohep = alohep;
		
		
		Border blackline = BorderFactory.createMatteBorder(1,1,2,1,Color.BLACK);
		BoxLayout layout = new BoxLayout(this, BoxLayout.Y_AXIS);
		
		setLayout(layout);
		setBorder(blackline);
		
		comboboxPanel = new JPanel();
		comboboxPanel.setLayout(new FlowLayout());
		JLabel parTypeLabel = new JLabel("Lepton/Hadron/Nucleus:");
		parTypeLabel.setFont(MainPanel.STANDARD_FONT);
		comboboxPanel.add(parTypeLabel);
		parTypeCbox = new JComboBox<String>();
		parTypeCbox.setFont(MainPanel.STANDARD_FONT);
		loadComboBox();
		accRadioButtons = new HashMap<String, JRadioButton>();
		parTypeCbox.addActionListener (new ActionListener () {
		    public void actionPerformed(ActionEvent e) {
				updateAccPanel();
				updateParPanel();
		    }
		});
		comboboxPanel.add(parTypeCbox);
		add(comboboxPanel);
		
		particleVarPanels = new HashMap<String, VariablePanel>();
		parVarMainPanel = new JPanel();
		accTypeMainPanel = new JPanel();
		accTypeMainPanel.setLayout(new FlowLayout());
		updateAccPanel();
		updateParPanel();
		add(parVarMainPanel);
		add(accTypeMainPanel);
		
	}
	public void loadComboBox()
	{
		alohep.getDataManager().getParticleDataMap().forEach((K, V) -> {
			parTypeCbox.addItem(K);
		});
	/*	for(Entry<String, ParticleData> en: alohep.getDataManager().getParticleDataMap().entrySet())
		{
			parTypeCbox.addItem(en.getKey());
		}*/
		parTypeCbox.setSelectedIndex(0);
	}
	public void updateAccPanel()
	{
		
		accTypeMainPanel.setLayout(new FlowLayout(5));
		for(Entry<String, JRadioButton> en: accRadioButtons.entrySet())
		{
			en.getValue().setSelected(false);
			en.getValue().setVisible(false);
		}
		ParticleData pd = alohep.getDataManager().getParticleDataMap().get(parTypeCbox.getSelectedItem());
		boolean selected = false;
		
		addRadioButton("Custom", true);
		for(Entry<String, AcceleratorData> en: pd.getAccelerators().entrySet())
		{
			selected = addRadioButton(en.getKey(), selected);
		}
	}
	public boolean addRadioButton(String key, boolean selected)
	{
		JRadioButton rb;
		if(!accRadioButtons.containsKey(key))
			rb = new JRadioButton(key);
		else
			rb = accRadioButtons.get(key);
		rb.setVisible(true);
		rb.setFont(MainPanel.STANDARD_FONT);
		if(!selected)
		{
			rb.setSelected(true);
			selected = true;
		}
		else
		{
			rb.setSelected(false);
		}
		rb.addActionListener(new ActionListener() 
		{
			
			@Override
			public void actionPerformed(ActionEvent e) 
			{
				
					for(Entry<String, JRadioButton> rb: accRadioButtons.entrySet())
					{
						if(e.getActionCommand().equals(rb.getValue().getActionCommand()))
						{
							rb.getValue().setSelected(true);
						}
						else
						{
							rb.getValue().setSelected(false);
						}
					}
					updateParPanel();
				}
		});
		accRadioButtons.put(key, rb);
		accTypeMainPanel.add(rb);
		
		return selected;
	}
	public void updatePanelEditable(boolean isEditable)
	{
		for(Entry<String, VariablePanel> en: particleVarPanels.entrySet())
		{
			en.getValue().getVarText().setEditable(isEditable);
		}
	}
	public void updateParPanel()
	{
		parVarMainPanel.setLayout(new BoxLayout(parVarMainPanel, BoxLayout.Y_AXIS));
		
		updatePanelEditable(true);
		
		if(getSelectedAcc().equals("Custom"))
		{
			updatePanelEditable(true);
			return;
		}
		
		
		for(Entry<String, VariablePanel> en: particleVarPanels.entrySet())
		{
			en.getValue().setVisible(false);
		}
		ParticleData pd = getParticleData();

		AcceleratorData ad = pd.getAccelerators().get(getSelectedAcc());
		
		for(Entry<String, Double> en: ad.getVariables().entrySet())
		{
			VariablePanel panel;
			if(!particleVarPanels.containsKey(en.getKey()))
			{
				panel = new VariablePanel(alohep, en.getKey());
				particleVarPanels.put(en.getKey(), panel);
				parVarMainPanel.add(panel);
			}
			else
			{
				panel = particleVarPanels.get(en.getKey());
				panel.setVisible(true);
			}
			panel.setValue(en.getValue());
		}
		if(!getSelectedAcc().equals("Custom"))
		{
			updatePanelEditable(false);
		}
	}
	public ParticleData getParticleData()
	{
		return alohep.getDataManager().getParticleDataMap().get(parTypeCbox.getSelectedItem().toString());
	}
	public String getSelectedAcc()
	{
		for(Entry<String, JRadioButton> rb: accRadioButtons.entrySet())
		{
			if(rb.getValue().isSelected())
			{
				return rb.getValue().getActionCommand();
			}
		}
		return "";
	}
	public void setRadioButton(String name)
	{
		for(Entry<String, JRadioButton> en: accRadioButtons.entrySet())
		{
			if(en.getKey().equals(name))
				en.getValue().setSelected(true);
			else
				en.getValue().setSelected(false);
		}
		if(getSelectedAcc().equals("Custom"))
		{
			updatePanelEditable(true);
		}
	}
	public String getSelectedParticle() {
		return parTypeCbox.getSelectedItem().toString();
	}
	public String getSelectedParticleType() {
		String str =  parTypeCbox.getSelectedItem().toString();
		return str.contains("-") ? str.substring(0,str.indexOf("-")) : str;
		
	}
	public AcceleratorData getAccData()
	{
		if(getSelectedAcc().equals("Custom"))
			return null;
		return alohep.getDataManager().getParticleDataMap().get(getSelectedParticle()).getAccelerators().get(getSelectedAcc());
	}
	public HashMap<String, VariablePanel> getParticleVarPanels() {
		return particleVarPanels;
	}
	public boolean hasParticleVarPanel(String key) {
		return particleVarPanels.containsKey(key) && particleVarPanels.get(key).isVisible();
	}
	
	public void setParticleVarPanels(HashMap<String, VariablePanel> particleVarPanels) {
		this.particleVarPanels = particleVarPanels;
	}

	
}
