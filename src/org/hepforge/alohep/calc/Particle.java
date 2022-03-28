package org.hepforge.alohep.calc;

public class Particle {

	private String name;
	private double charge;
	private double mass;
	private double lifetime;
	private double rad;
	
	public Particle(String name, double charge, double mass, double rad)
	{
		this.name = name;
		this.charge = charge;
		this.mass = mass;
		this.rad = rad;
		lifetime = -1;
	}

	public Particle(String name, double charge, double mass, double rad, double lifetime)
	{
		this.name = name;
		this.charge = charge;
		this.mass = mass;
		this.rad = rad;
		this.lifetime = lifetime;
	}
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public double getCharge() {
		return charge;
	}

	public void setCharge(double charge) {
		this.charge = charge;
	}

	public double getMass() {
		return mass;
	}

	public void setMass(double mass) {
		this.mass = mass;
	}
	
	public double getRad() {
		return rad;
	}

	public void setRad(double rad) {
		this.rad = rad;
	}

	public double getLifetime() {
		return lifetime;
	}

	public void setLifetime(double lifetime) {
		this.lifetime = lifetime;
	}
}
