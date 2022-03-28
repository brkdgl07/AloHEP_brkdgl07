package org.hepforge.alohep.calc;

public class Vector4d {

	public double t,x,y,z;
	
	public static Vector4d MaxValue = new Vector4d(100000, 100000, 100000, 100000);
	
	public Vector4d(double t, double x, double y, double z)
	{
		this.t = t;
		this.x = x;
		this.y = y; 
		this.z = z;
	}
	public Vector4d(Vector4d vec)
	{
		this.t = vec.t;
		this.x = vec.x;
		this.y = vec.y; 
		this.z = vec.z;
	}
	public Vector4d copy()
	{
		return new Vector4d(this);
	}
	public double diagonal()
	{
		return Math.pow(t*t+x*x+y*y+z*z,0.5);
	}
	public double volume()
	{
		return t*x*y*z;
	}
	public Vector4d set(Vector4d vec) 
	{
		t = vec.t;
		x = vec.x;
		y = vec.y;
		z = vec.z;
		return this;
	}
	public Vector4d add(Vector4d vec)
	{
		t += vec.t;
		x += vec.x;
		y += vec.y;
		z += vec.z;
		return this;
	}
	public Vector4d subtract(Vector4d vec)
	{
		t -= vec.t;
		x -= vec.x;
		y -= vec.y;
		z -= vec.z;
		return this;
	}

	public Vector4d multiply(Vector4d vec)
	{
		t = t*vec.t;
		x = x*vec.x;
		y = y*vec.y;
		z = z*vec.z;
		return this;
	}
	public Vector4d divide(Vector4d vec)
	{
		t = t/vec.t;
		x = x/vec.x;
		y = y/vec.y;
		z = z/vec.z;
		return this;
	}
	public Vector4d multiply(double k)
	{
		t*=k;
		x*=k;
		y*=k;
		z*=k;
		return this;
	}
	public Vector4d pow(double k)
	{
		t = Math.pow(t, k);
		x = Math.pow(x, k);
		y = Math.pow(y, k);
		z = Math.pow(z, k);
		return this;
	}
	public Vector4d exp()
	{
		t = Math.pow(Math.E, t);
		x = Math.pow(Math.E, x);
		y = Math.pow(Math.E, y);
		z = Math.pow(Math.E, z);
		
		return this;
	}
	public Vector4d abs()
	{
		t = Math.abs(t);
		x = Math.abs(x);
		y = Math.abs(y);
		z = Math.abs(z);
		return this;
	}
	public Vector3d pos()
	{
		return new Vector3d(x,y,z);
	}
	
	public double scalarMultiply(Vector4d vec)
	{
		return t*vec.t+x * vec.x + y * vec.y + z * vec.z;	
	}
	public boolean isNaN()
	{
		return Double.isNaN(t) || Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(z);
	}
	public boolean isInfinity()
	{
		return Double.isInfinite(t) || Double.isInfinite(x) || Double.isInfinite(y) || Double.isInfinite(z);
	}

}
