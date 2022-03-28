package org.hepforge.alohep.calc;

public class Vector2d {

	public double x,y;
	
	public Vector2d(double x, double y)
	{
		this.x = x;
		this.y = y; 
	}
	public double diagonal()
	{
		return Math.pow(x*x+y*y,0.5);
	}
	public void set(Vector2d vec) {
		x = vec.x;
		y = vec.y;
	}
}
