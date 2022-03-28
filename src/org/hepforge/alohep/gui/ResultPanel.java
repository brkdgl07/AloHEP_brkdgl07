package org.hepforge.alohep.gui;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.text.DecimalFormat;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.hepforge.alohep.AloHEP;
import org.hepforge.alohep.calc.MathHelper;
import org.hepforge.alohep.database.ParticleData;

public class ResultPanel extends JPanel{

	private AloHEP alohep;

	public ResultPanel(AloHEP alohep)
	{
		this.alohep = alohep;
	}
	public void finalResults() {
	    

        DecimalFormat format = new DecimalFormat("0.###");
        // Luminosity Panel
		BorderLayout resultPanelLayout = new BorderLayout();

		setLayout(resultPanelLayout);

		
        JPanel luminosityPanel = new JPanel();  
		BoxLayout luminosityLayout = new BoxLayout(luminosityPanel, BoxLayout.Y_AXIS);
		luminosityPanel.setLayout(luminosityLayout);

    	VariablePanel noLum = new VariablePanel(alohep, "noLum", false);
    	noLum.setValue( MathHelper.parseToScientificNotation(alohep.getLuminosityCalc().getLuminosityRaw()));	
    	luminosityPanel.add(noLum);
    	
    	VariablePanel efLum = new VariablePanel(alohep, "efLum", false);
    	efLum.setValue( MathHelper.parseToScientificNotation(alohep.getLuminosityCalc().getLuminosity()));	
    	luminosityPanel.add(efLum);
    	
    	VariablePanel reFac = new VariablePanel(alohep, "reFac", false);
    	reFac.setValue( format.format(efLum.getValue() / noLum.getValue()));	
    	luminosityPanel.add(reFac);
    	
    	
        add(luminosityPanel,BorderLayout.NORTH);
        
        JPanel particlePanel = new JPanel();
		GridLayout particleLayout = new GridLayout(1,2);
		particlePanel.setLayout(particleLayout);
		
        JPanel leftPanel = new JPanel();
		BoxLayout leftLayout = new BoxLayout(leftPanel, BoxLayout.Y_AXIS);

        leftPanel.setLayout(leftLayout);
        JLabel leftParticleName = new JLabel(alohep.getMainPanel().getLeftPanel().getSelectedParticle());
        leftParticleName.setFont(MainPanel.TITLE_FONT);
        leftPanel.add(leftParticleName);
        
        if(alohep.getMainPanel().getLeftPanel().getParticleData().getType() == ParticleData.Type.LINEAR)
        {
        	VariablePanel disX = new VariablePanel(alohep, "disX", false);
        	VariablePanel disY = new VariablePanel(alohep, "disY", false);
        	
        	double[]dis = alohep.getLuminosityCalc().getDisruption(alohep.getLuminosityCalc().getBeamRight(), alohep.getLuminosityCalc().getBeamLeft());
            
        	disX.setValue(format.format(dis[0]));
        	disY.setValue(format.format(dis[1]));
        	
        	leftPanel.add(disX);
        	leftPanel.add(disY);
        	
        	VariablePanel divX = new VariablePanel(alohep, "divX", false);
        	VariablePanel divY = new VariablePanel(alohep, "divY", false);
        	
            double[]div = alohep.getLuminosityCalc().getDivergence(alohep.getLuminosityCalc().getBeamRight());
            
        	divX.setValue(format.format(div[0]));
        	divY.setValue(format.format(div[1]));
        	
        	leftPanel.add(divX);
        	leftPanel.add(divY);
        }
        else
        {
        	VariablePanel bbTunX = new VariablePanel(alohep, "bbTunX", false);
        	VariablePanel bbTunY = new VariablePanel(alohep, "bbTunY", false);
        	
            double[]bbTun = alohep.getLuminosityCalc().getBeamBeam(alohep.getLuminosityCalc().getBeamRight(), alohep.getLuminosityCalc().getBeamLeft());
            
            bbTunX.setValue(format.format(bbTun[0]));
            bbTunY.setValue(format.format(bbTun[1]));
        	
        	leftPanel.add(bbTunX);
        	leftPanel.add(bbTunY);
        }
        particlePanel.add(leftPanel);
        

        JPanel rightPanel = new JPanel();
		BoxLayout rightLayout = new BoxLayout(rightPanel, BoxLayout.Y_AXIS);
        rightPanel.setLayout(rightLayout);
        
        JLabel rightParticleName = new JLabel(alohep.getMainPanel().getRightPanel().getSelectedParticle());
        rightParticleName.setFont(MainPanel.TITLE_FONT);
        rightPanel.add(rightParticleName);        
        if(alohep.getMainPanel().getRightPanel().getParticleData().getType() == ParticleData.Type.LINEAR)
        {
        	VariablePanel disX = new VariablePanel(alohep, "disX", false);
        	VariablePanel disY = new VariablePanel(alohep, "disY", false);
        	
        	double[]dis = alohep.getLuminosityCalc().getDisruption(alohep.getLuminosityCalc().getBeamLeft(), alohep.getLuminosityCalc().getBeamRight());
            
        	disX.setValue(format.format(dis[0]));
        	disY.setValue(format.format(dis[1]));
        	
        	rightPanel.add(disX);
        	rightPanel.add(disY);
        	
        	VariablePanel divX = new VariablePanel(alohep, "divX", false);
        	VariablePanel divY = new VariablePanel(alohep, "divY", false);
        	
            double[]div = alohep.getLuminosityCalc().getDivergence(alohep.getLuminosityCalc().getBeamLeft());
            
        	divX.setValue(format.format(div[0]));
        	divY.setValue(format.format(div[1]));
        	
        	rightPanel.add(divX);
        	rightPanel.add(divY);
        }
        else
        {
        	VariablePanel bbTunX = new VariablePanel(alohep, "bbTunX", false);
        	VariablePanel bbTunY = new VariablePanel(alohep, "bbTunY", false);
        	
            double[]bbTun = alohep.getLuminosityCalc().getBeamBeam(alohep.getLuminosityCalc().getBeamLeft(), alohep.getLuminosityCalc().getBeamRight());
            
            bbTunX.setValue(format.format(bbTun[0]));
            bbTunY.setValue(format.format(bbTun[1]));
        	
        	rightPanel.add(bbTunX);
        	rightPanel.add(bbTunY);
        }
        particlePanel.add(rightPanel);
        add(particlePanel, BorderLayout.SOUTH);
        
        JFrame frame = new JFrame();
        frame.setContentPane(this);
        frame.setTitle("FINAL RESULTS");
        frame.setSize(520,340);
        frame.setVisible(true);
        frame.setResizable(true);
        frame.setLocationRelativeTo(null);

		
	}
}
