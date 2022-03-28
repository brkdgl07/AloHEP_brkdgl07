package org.hepforge.alohep.gfx;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.text.DecimalFormat;

import org.hepforge.alohep.AloHEP;
import org.hepforge.alohep.calc.Bunch;
import org.hepforge.alohep.calc.LuminosityCalc;
import org.hepforge.alohep.calc.Vector3d;

public class MainState extends State {

	private AloHEP alohep;
	private GaussModel gaussModel;
	Bunch bunchLeft;
	Bunch bunchRight;
	Vector3d leftpos[][];
	Vector3d rightpos[][];
	GaussModel.Gauss leftGauss[];
	GaussModel.Gauss rightGauss[];
	private int mode = 1;
	private int step = 0;
	public MainState(AloHEP alohep) {
		super(alohep.getAnimationPanel());
		this.alohep = alohep;
	}

	@Override
	public void init(AnimationPanel animationPanel) {
		gaussModel = new GaussModel();
		leftpos = alohep.getLuminosityCalc().getBeamLeft().getBunch().getAllPositions();
		leftGauss = new GaussModel.Gauss[leftpos[0].length];
		rightpos = alohep.getLuminosityCalc().getBeamRight().getBunch().getAllPositions();
		rightGauss = new GaussModel.Gauss[rightpos[0].length];
		for(int j = 0; j < leftpos[0].length; j++)
		{
			leftGauss[j] = gaussModel.new Gauss(new Point(0,0), 100);
		}		
		for(int j = 0; j < leftpos[0].length; j++)
		{
			rightGauss[j] = gaussModel.new Gauss(new Point(0,0), 100);
		}
		
	}

	@Override
	public void update() {
		if(getinput().isKeyPressing(KeyEvent.VK_RIGHT))
		{
			if(step<LuminosityCalc.NumSteps-1)
				step++;
		}
		if(getinput().isKeyPressing(KeyEvent.VK_LEFT))
		{
			if(step>0)
				step--;
		}
		if(getinput().isKeyClicked(KeyEvent.VK_X))
		{
			getCamera().setAxis(Camera.Axis.X);
		}
		
		if(getinput().isKeyClicked(KeyEvent.VK_Y))
		{
			getCamera().setAxis(Camera.Axis.Y);
		}
		
		if(getinput().isKeyClicked(KeyEvent.VK_Z))
		{
			getCamera().setAxis(Camera.Axis.Z);
		}
		
		if(getinput().isKeyClicked(KeyEvent.VK_M))
		{
			mode = (mode+1) % 2;
		}
		if(getinput().getWheel()>0)
		{
			getCamera().scaleChange(1.1);
			getinput().setWheel(0);
		}
		if(getinput().getWheel()<0)
		{
			getCamera().scaleChange(0.9);
			getinput().setWheel(0);

		}
	}

	
	@Override
	public void render(BufferedImage image) {
		 Graphics2D g = image.createGraphics();
	        g.setComposite(AlphaComposite.Src);
	        g.setColor(new Color(0,0,0,255));
	        g.fillRect(0,0,image.getWidth(),image.getHeight());
		for(int j = 0; j < leftpos[0].length; j++)
		{
			Point p = board.getCamera().getCameraPos(leftpos[step][j]);

			leftGauss[j].setP(p);
			switch(mode)
			{
			case 0:
				gaussModel.drawGauss(g, leftGauss[j]);
				break;
			case 1:
				g.setComposite(AlphaComposite.Src);
				g.setColor(Color.RED);
				g.fillOval(p.x, p.y, 10, 10);
				break;
				default:
					break;
			}
			

		}
		
		for(int j = 0; j < rightpos[0].length; j++)
		{
			Point p = board.getCamera().getCameraPos(rightpos[step][j]);

			rightGauss[j].setP(p);

			switch(mode)
			{
			case 0:
				gaussModel.drawGauss(g, rightGauss[j]);
				break;
			case 1:
				g.setComposite(AlphaComposite.Src);
				g.setColor(Color.BLUE);
				g.fillOval(p.x, p.y, 10, 10);
				break;
				default:
					break;
			}
		}
		g.setComposite(AlphaComposite.Src);
		g.drawString(getCamera().getAxis() + " Axis", Display.width-80, 20);
		g.drawString("Frame: "+(step+1)+"/"+leftpos.length, Display.width-80, 40);
		g.drawString("Scale x"+new DecimalFormat("0.0").format(getCamera().getScaleFactor()), Display.width-80, 60);

        g.dispose();
	}

}
