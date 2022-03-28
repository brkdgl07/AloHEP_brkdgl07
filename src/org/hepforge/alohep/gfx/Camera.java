package org.hepforge.alohep.gfx;

import java.awt.Point;

import org.hepforge.alohep.AloHEP;
import org.hepforge.alohep.calc.LuminosityCalc;
import org.hepforge.alohep.calc.Vector3d;

public class Camera {

	private double x;
	private double y;
	private Vector3d scale;
	private double scaleFactor = 1;
	private AloHEP alohep;
	private AnimationPanel board;
	private Axis axis = Axis.X;
	
	public static enum Axis{X, Y, Z}
	
	public Camera(AloHEP alohep)
	{
		this.alohep = alohep;
		this.board = alohep.getAnimationPanel();
		scale = Vector3d.unitVec.copy().divide(LuminosityCalc.size).multiply(0.5);
		
	}
	public Point getCameraPos(Vector3d pos)
	{
		int x = 0, y = 0;
		switch(axis)
		{
		case X:
			x = (int)((pos.z*scale.z*scaleFactor+0.5)*Display.height+(Display.width-Display.height)/2);
			x = x < 0 ? 1: x;
			y = (int)((pos.y*scale.y*scaleFactor+0.5)*Display.height);
			y = y < 0 ? 1: y;
			break;
		case Y:
			x = (int)((pos.z*scale.z*scaleFactor+0.5)*Display.height+(Display.width-Display.height)/2);
			x = x < 0 ? 1: x;
			y = (int)((pos.x*scale.x*scaleFactor+0.5)*Display.height);
			y = y < 0 ? 1: y;

			break;
		case Z:
			x = (int)((pos.x*scale.x*scaleFactor+0.5)*Display.height+(Display.width-Display.height)/2);
			x = x < 0 ? 1: x;
			y = (int)((pos.y*scale.y*scaleFactor+0.5)*Display.height);
			y = y < 0 ? 1: y;
			break;
		}
		
		return new Point(x, y);
	}

	public Axis getAxis() {
		return axis;
	}
	public void setAxis(Axis axis) {
		this.axis = axis;
	}
	public Vector3d getScale() {
		return scale;
	}
	public void scaleChange(double k)
	{
		scaleFactor*=k;
	}
	public void setScale(Vector3d scale) {
		this.scale = scale;
	}
	public double getScaleFactor() {
		return scaleFactor;
	}
	public void setScaleFactor(double scaleFactor) {
		this.scaleFactor = scaleFactor;
	}
	
}
