package org.hepforge.alohep.calc;

import org.hepforge.alohep.AloHEP;

public class Bunch {

	public static double c = 2.99e8;
	public static double Coulomb = 1.602e-19;
	public static double eVtoKg = 1.782661921e-36;
	public static double e0 = 8.854187817e-12;
	private AloHEP alohep;
	private Beam beam;
	public static int SigmaResolution = 10;
	public static double SigmaDivider =  1.6;   //8 -> 1.417259;
	public static int Size = (int)Math.pow(SigmaResolution,3);
	public static Vector3d Scale = new Vector3d(3,3,3);
	private SubBunch subBunches[];
	private Slice slices[];
	private Vector3d allPositions[][];
	private Vector3d position;
	private Vector3d velocity;
	private Vector3d sigma;
	private Vector3d subSigma;
	private Vector3d transformedSubSigma;
	private Vector2d emmitance;
	private Vector2d beta;
	private double gamma;
	private double numberOfParticles;
	private double chargeOfParticle;
	private double chargeOfBunch;
	private double massOfParticle;
	private Particle particle;
	private double numberOfCollision = 0;
	private double rParticle;
	public Bunch(AloHEP alohep, Beam beam)
	{
		this.alohep = alohep;
		this.beam = beam;
	}
	
	public void init()
	{
		
		emmitance = new Vector2d(0, 0);
		numberOfParticles = getData("NumParInBun");
		
		beta = new Vector2d(getData("BetaHor"), getData("BetaVer"));
		gamma = getData("EnBeam") * 1.0e9 / (beam.getParticle().getMass());
		emmitance.x = getData("EmNorHor") / gamma;
		emmitance.y = getData("EmNorVer") / gamma;
		double sigmaX = Math.sqrt(emmitance.x * beta.x);
		double sigmaY = Math.sqrt(emmitance.y * beta.y);
		sigma = new Vector3d(sigmaX, sigmaY, getData("BunLen") / 2);
		particle = beam.getParticle();
		
		chargeOfParticle = particle.getCharge() * Coulomb;//  to Coulomb
					
		chargeOfBunch = numberOfParticles * chargeOfParticle;
		massOfParticle = particle.getMass() * eVtoKg;  // eV to kg
		
	}

	public void initsubBunches()
	{
		subBunches = new SubBunch[Size];
		subSigma = sigma.copy().multiply(1/SigmaDivider);
		int count = 0;
		double xPos[] = getLinePosition(sigma.x);
		double yPos[] = getLinePosition(sigma.y);
		double zPos[] = getLinePosition(sigma.z);

		for(int z = 0; z < zPos.length; z++)
		{
			for(int y = 0; y < yPos.length; y++)
			{
				for(int x = 0; x < xPos.length; x++)
				{
					Vector3d position = new Vector3d(this.position.x+xPos[x], this.position.y+yPos[y],this.position.z+zPos[z]);
					Vector3d velocity = this.velocity.copy();
					subBunches[count] = new SubBunch(this, position, velocity);
					count++;
				}
			}
		}
		allPositions = new Vector3d[LuminosityCalc.NumSteps][Size];
		for(count = 0; count < Size; count++)
		{
			subBunches[count].init();
		}
		transformedSubSigma = SubBunch.lengthContraction(subSigma, velocity);
		rParticle = Math.pow(chargeOfParticle,2)/(4*Math.PI*e0*massOfParticle*c*c);

	}
	
	public void initCollision(Bunch bunch)
	{
		for(int count = 0; count < Size; count++)
		{
			subBunches[count].initCollision(bunch);
		}
	}
	
	public void update(Bunch bunch, double dt)
	{

		for(int count = 0; count < Size; count++)
		{
			subBunches[count].update(bunch, dt);
			
			numberOfCollision += dt/(Size*Size)*subBunches[count].getCollisionDensity();
			subBunches[count].setCollisionDensity(0);
		}
	}
	public void updatePosition(int step, double dt, Bunch otherBunch)
	{
		for(int count = 0; count < Size; count++)
		{
			allPositions[step][count] = subBunches[count].position.copy();
			subBunches[count].updatePosition(dt, otherBunch);
		}
	}
	public double[] getLinePosition(double sigma)
	{
		int resolution = 1000;
		int gaussScale = 5;
		double[] pos = new double[SigmaResolution];
		double subArea = 1.0 / SigmaResolution;
		double diffArea = subArea;
		double xstep = sigma*gaussScale/resolution;
		double x = 0;
		int count = SigmaResolution/2;
		for(int i = 0; i < resolution/2; i++)
		{
			if(diffArea >= subArea)
			{
				pos[count] = x;
				pos[SigmaResolution-count-1] = -x;
				count++;
				diffArea = 0;
			}
			diffArea += xstep/(Math.pow(2*Math.PI, 0.5)*sigma)* Math.exp(-x*x/(2*sigma*sigma));
			x+=xstep;
		}
		
		return pos;
		
	}
	
	
	public double getSettings(String key)
	{
		return alohep.getDataManager().getSettingsData().get(key);
	}
	
	public double getData(String key)
	{
		return beam.getData(key);
	}

	public Vector3d getSigma()
	{
		return sigma;
	}
	
	public Vector3d getVelocity() {
		return velocity;
	}
	
	public void setVelocity(Vector3d velocity) {
		this.velocity = velocity;
	}
	
	public Vector3d getPosition() {
		return position;
	}
	
	public void setPosition(Vector3d position) {
		this.position = position;
	}

	public double getGamma() {
		return gamma;
	}

	public void setGamma(double gamma) {
		this.gamma = gamma;
	}

	public Vector2d getEmmitance() {
		return emmitance;
	}

	public void setEmmitance(Vector2d emmitance) {
		this.emmitance = emmitance;
	}

	public Vector2d getBeta() {
		return beta;
	}

	public void setBeta(Vector2d beta) {
		this.beta = beta;
	}

	public double getNumberOfParticles() {
		return numberOfParticles;
	}

	public void setNumberOfParticles(double numberOfParticles) {
		this.numberOfParticles = numberOfParticles;
	}

	public double getChargeOfParticle() {
		return chargeOfParticle;
	}

	public void setChargeOfParticle(double chargeOfParticle) {
		this.chargeOfParticle = chargeOfParticle;
	}

	public double getChargeOfBunch() {
		return chargeOfBunch;
	}

	public void setChargeOfBunch(double chargeOfBunch) {
		this.chargeOfBunch = chargeOfBunch;
	}

	public double getMassOfParticle() {
		return massOfParticle;
	}

	public void setMassOfParticle(double massOfParticle) {
		this.massOfParticle = massOfParticle;
	}

	public SubBunch[] getsubBunches() {
		return subBunches;
	}

	public void setsubBunches(SubBunch[] subBunches) {
		this.subBunches = subBunches;
	}

	public Vector3d getSubSigma() {
		return subSigma;
	}

	public void setSubSigma(Vector3d subSigma) {
		this.subSigma = subSigma;
	}

	public double getNumberOfCollision() {
		return numberOfCollision;
	}

	public void setNumberOfCollision(double numberOfCollision) {
		this.numberOfCollision = numberOfCollision;
	}

	public Vector3d[][] getAllPositions() {
		return allPositions;
	}

	public void setAllPositions(Vector3d[][] allPositions) {
		this.allPositions = allPositions;
	}

	public Vector3d getTransformedSubSigma() {
		return transformedSubSigma;
	}

	public void setTransformedSubSigma(Vector3d transformedSubSigma) {
		this.transformedSubSigma = transformedSubSigma;
	}

	public double getRParticle() {
		return rParticle;
	}

	public void setRParticle(double rParticle) {
		this.rParticle = rParticle;
	}

	public Beam getBeam() {
		return beam;
	}

	public void setBeam(Beam beam) {
		this.beam = beam;
	}
}
