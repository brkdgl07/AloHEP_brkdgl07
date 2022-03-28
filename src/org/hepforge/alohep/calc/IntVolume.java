package org.hepforge.alohep.calc;

public class IntVolume {

	private Vector3d center;
	private Vector3d size;
	private Vector3d dq;
	
	public IntVolume(Vector3d center, Vector3d size)
	{
		this.center = center;
		this.size = size;
		dq = new Vector3d(0, 0, 0);
	}
	public void updateCenter(double x, double y, double z)
	{
		this.center.x = x;
		this.center.y = y;
		this.center.z = z;
	}
	public void updateSize(double x, double y, double z)
	{
		this.size.x = x;
		this.size.y = y;
		this.size.z = z;
	}
	public void updateDq(double res)
	{
		dq.x = this.size.x/res;
		dq.y = this.size.y/res;
		dq.z = this.size.z/res;
	}
	public Vector3d getCenter() {
		return center;
	}
	public void setCenter(Vector3d center) {
		this.center = center;
	}
	public Vector3d getSize() {
		return size;
	}
	public void setSize(Vector3d size) {
		this.size = size;
	}
	public Vector3d getDq() {
		return dq;
	}
	public void setDq(Vector3d dq) {
		this.dq = dq;
	}
}
