package org.hepforge.alohep.gfx;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.MultipleGradientPaint.CycleMethod;
import java.awt.Point;
import java.awt.RadialGradientPaint;
import java.awt.geom.Point2D;

public class GaussModel {

	public static final int SHADOW_SIZE = 1;


	public void drawGauss(Graphics2D g, Gauss l)
	{  
	    g.setComposite(AlphaComposite.DstOut);
	    Point2D center = new Point2D.Float(l.p.x, l.p.y);
	    float[] dist = {0f, 1.0f};
	    Color[] colors = {new Color(255,255,255,5), new Color(0,0,0,0)};
	    RadialGradientPaint p = new RadialGradientPaint(center, l.radius, dist, colors, CycleMethod.NO_CYCLE);
	    g.setPaint(p);
	    g.fillOval(l.p.x-(int)l.radius,l.p.y-(int)l.radius,(int)l.radius*2,(int)l.radius*2);
	}
	public class Gauss
	{
		private Point p;
		private float radius;
		public Gauss(Point p,float radius)
		{
			this.p = new Point(p.x/GaussModel.SHADOW_SIZE,p.y/GaussModel.SHADOW_SIZE);
			this.radius = radius/GaussModel.SHADOW_SIZE;
		}
		public Point getP() {
			return p;
		}
		public void setP(Point p) {
			this.p = p;
		}
	}
}
