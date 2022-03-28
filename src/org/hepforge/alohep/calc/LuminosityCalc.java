package org.hepforge.alohep.calc;

import java.util.HashMap;


import org.hepforge.alohep.AloHEP;

public class LuminosityCalc {

	private HashMap<String, Particle> particles;  
	
	private AloHEP alohep;
	
	
	public static final double c = 2.99e8;
	public static final double eps0 = 8.85419e-12;
	public static final double mu0 = 4 * Math.PI * 1e-7;
	public static final double eQ = 1.602e-19;
	private static final double k = 1 / (4 * Math.PI * eps0);
	
	
	public static final double Scale = 4;
	public static final double IntScale = 5;
	public static final double gaussScale = Scale; //farkli bir deger olursa calismiyor. Scale ne kadar buyukse %99.9999 o kadar yakin sonuc!	
	public static final int RES_g = 16;
	public static final int RES_g3 = RES_g*RES_g*RES_g; //gaussian dagitilan parcacik cozunurlugu
	public static int NumSteps;
	public double []particlesLeftX = new double[RES_g3];
	public double []particlesLeftY = new double[RES_g3];
	public double []particlesLeftZ = new double[RES_g3];
	public double []particlesRightX = new double[RES_g3];
	public double []particlesRightY = new double[RES_g3];
	public double []particlesRightZ = new double[RES_g3];
	
	public double []particlesLeftVX = new double[RES_g3];
	public double []particlesLeftVY = new double[RES_g3];
	public double []particlesLeftVZ = new double[RES_g3];
	public double []particlesRightVX = new double[RES_g3];
	public double []particlesRightVY = new double[RES_g3];
	public double []particlesRightVZ = new double[RES_g3];
	
	
	public double []particlesLeftPX = new double[RES_g3];
	public double []particlesLeftPY = new double[RES_g3];
	public double []particlesLeftPZ = new double[RES_g3];
	public double []particlesRightPX = new double[RES_g3];
	public double []particlesRightPY = new double[RES_g3];
	public double []particlesRightPZ = new double[RES_g3];
	

	
	public double []gaussMassL       = new double[RES_g3];
	public double []gaussMassR       = new double[RES_g3];
	
	public double [][]ULeft  = new double[RES_g3][4];
	public double [][]URight = new double[RES_g3][4];
	
	public double [][]Lambda_inv_L       = new double[4][4];
	public double [][]Lambda_inv_R       = new double[4][4];
	public double [][]Fields_L_lab         = new double[4][4];
	public double [][]Fields_R_lab         = new double[4][4];
	public double [][]Fields_L       = new double[4][4];
	public double [][]Fields_R       = new double[4][4];
	
	private Vector3d EFieldLeft = new Vector3d(0,0,0);
	private Vector3d BFieldLeft = new Vector3d(0,0,0);
	private Vector3d EFieldRight = new Vector3d(0,0,0);
	private Vector3d BFieldRight = new Vector3d(0,0,0);
	
	private Vector3d FtoLeft= new Vector3d(0,0,0);
	private Vector3d FtoRight= new Vector3d(0,0,0);

	
	private static final int RES_q = 10;
	public static Vector3d size;
	private IntVolume intVolume;
	
	private double totalTime;
	private double presentTime;
	private double dTime;
	private double luminosity;
	double dLxyz;
	
	private Beam beamLeft;
	private Beam beamRight;
	
	public LuminosityCalc(AloHEP alohep)
	{
		this.alohep = alohep;
		initParticleInterface();
	}
	
	public void calculateLuminosity()
	{
		
	}
	
	public void initBeams()
	{
		beamLeft = new Beam(alohep, AloHEP.PanelType.LEFT);
		beamRight = new Beam(alohep, AloHEP.PanelType.RIGHT);
		beamLeft.initBeam();
		beamRight.initBeam();
	}
	public void initCollision()
	{
		double lengthX = Math.max(beamLeft.getBunch().getSigma().x, beamRight.getBunch().getSigma().x) * Scale;
		double lengthY = Math.max(beamLeft.getBunch().getSigma().y, beamRight.getBunch().getSigma().y) * Scale;
		double lengthZ = Math.max(beamLeft.getData("BunLen"), beamRight.getData("BunLen")) * Scale;
		intVolume = new IntVolume(new Vector3d(lengthX, lengthY, lengthZ), new Vector3d(0, 0, 0));   
		size = new Vector3d(lengthX, lengthY, lengthZ);
		totalTime = size.diagonal() / c;
		presentTime = 0;
		luminosity = 0;
		NumSteps = alohep.getDataManager().getSettingsData().get("NumSteps").intValue();
		
		dTime = totalTime / NumSteps;
		double phi = alohep.getDataManager().getSettingsData().get("Angle");
		
		double sin = Math.sin(phi/2);
		double cos = Math.cos(phi/2);

		
		beamLeft.getBunch().setPosition(new Vector3d(0,lengthZ/2*sin, -lengthZ/2*cos));
		beamRight.getBunch().setPosition(new Vector3d(0,lengthZ/2*sin, lengthZ/2*cos+lengthZ*dTime/2));
		
		double gammaL = beamLeft.getBunch().getGamma();
		double gammaR = beamRight.getBunch().getGamma();
		double speedL = c*Math.sqrt(gammaL*gammaL-1)/gammaL;
		double speedR = c*Math.sqrt(gammaR*gammaR-1)/gammaR;

		beamLeft.getBunch().setVelocity(new Vector3d(0,speedL*sin , speedL*cos));
		beamRight.getBunch().setVelocity(new Vector3d(0,speedR*sin , -speedR*cos));
		
		beamLeft.getBunch().initsubBunches();
		beamRight.getBunch().initsubBunches();
		beamLeft.getBunch().initCollision(beamRight.getBunch());
		beamRight.getBunch().initCollision(beamLeft.getBunch());
		
	}
	
	public double updateCollision(int step)
	{
		//locateParticle3d();
		//spaceIntegration2();
		luminosity += dLxyz* dTime;
		presentTime += dTime;

		beamLeft.getBunch().update(beamRight.getBunch(), dTime);
		beamRight.getBunch().update(beamLeft.getBunch(), dTime);

		beamRight.getBunch().updatePosition(step, dTime, beamLeft.getBunch());
		beamLeft.getBunch().updatePosition(step, dTime, beamRight.getBunch());
		
		
		return presentTime / totalTime;
	}
	
	public double getRo(double sigma, double q)
	{
		return Math.pow(Math.E, -(q * q) / (2 * sigma * sigma)) / sigma;
	}

	public double getLuminosity()
	{
		double NL = beamLeft.getData("NumParInBun");
		double NR = beamRight.getData("NumParInBun");
		double ltL = beamLeft.getParticle().getLifetime();
		double ltR = beamRight.getParticle().getLifetime();
		double circumL = beamLeft.hasData("Circum") ? beamLeft.getData("Circum") : 1;
		double circumR = beamRight.hasData("Circum") ? beamRight.getData("Circum") : 1;
		double circleL = beamLeft.hasData("Circle") ? beamLeft.getData("Circle") : Double.MAX_VALUE;
		double circleR = beamRight.hasData("Circle") ? beamRight.getData("Circle") : Double.MAX_VALUE;
		
		double dltL = beamLeft.getBunch().getGamma() * beamLeft.getParticle().getLifetime();
		double dltR = beamRight.getBunch().getGamma() * beamRight.getParticle().getLifetime();

		double tTimeL = circumL * 1000 * circleL / c;
		double tTimeR = circumR * 1000 * circleR / c;
		
		double tTime = Math.min(tTimeL, tTimeR);
		
		double Nsquare = NL * NR;

		if(ltR != -1 && ltL != -1)
		{

			Nsquare = NL * NR * dltL * dltR / (dltL + dltR) * (-Math.pow(Math.E, tTime * -(1 / dltL + 1 / dltR)) + 1) /tTime;

		}
		else if(ltL != -1)
		{
			NL = NL * dltL * (-Math.pow(Math.E, -tTimeL / dltL) + 1) / tTimeL; 
			
			Nsquare = NL * NR;
		}
		else if(ltR != -1)
		{
			NR = NR * dltR * (-Math.pow(Math.E, -tTimeR / dltR) + 1) / tTimeR; 
			Nsquare = NL * NR;
		}
		double fc = Math.min(getFCollision(beamLeft), getFCollision(beamRight));

		
		double phi = alohep.getDataManager().getSettingsData().get("Angle");
		double cos = Math.cos(phi/2);
		System.out.println(beamLeft.getBunch().getNumberOfCollision());
		System.out.println("SubBunch Luminosity: " + (2 * c * Nsquare * fc / (8 * Math.pow(Math.PI,3))*beamLeft.getBunch().getNumberOfCollision() * 1e-4));
		System.out.println("Subsigma: "+beamLeft.getBunch().getSubSigma().print());
		return 2 * c * cos * cos * Nsquare * fc / (8 * Math.pow(Math.PI,3)) * luminosity * 1e-4; // (cm^-2)*(s^-1) 
	}
	public double getFCollision(Beam beam)
	{
		double fc = 0;
		if(beam.hasData("ColFrq"))
		{
			fc = beam.getData("ColFrq");
			return fc;
		}
		else if(beam.hasData("RevFrq"))
		{
			fc = beam.getData("RevFrq");
		}
		else if(beam.hasData("PulFrq"))
		{
			fc = beam.getData("PulFrq");
		}

		else if(beam.hasData("RepFrq"))
		{
			fc = beam.getData("RepFrq") * beam.getData("Circle");
		}
		fc *= beam.getData("NumBunInBeam");
		
		return fc;
	}
	public double getLuminosityRaw()
	{

		double fc = Math.min(getFCollision(beamLeft), getFCollision(beamRight));
		
		double NL = beamLeft.getBunch().getData("NumParInBun");
		double NR = beamRight.getBunch().getData("NumParInBun");
		Vector3d sigmaL = beamLeft.getBunch().getSigma();
		Vector3d sigmaR = beamRight.getBunch().getSigma();
		return fc * NL * NR / (4 * Math.PI * Math.max(sigmaL.x, sigmaR.x) * Math.max(sigmaL.y, sigmaR.y)) * 1e-4;
	}
	
	public double[] getDivergence(Beam B) {
		
		double [] A = new double [2];
		A[0]=0.000;
		A[1]=0.000;

		A[0] = B.getBunch().getSigma().z / B.getBunch().getData("BetaHor") ;
		A[1] = B.getBunch().getSigma().z / B.getBunch().getData("BetaVer") ;
		return A;
	
	}
	
	public double [] getDisruption(Beam A, Beam B) {
		double []  D = new double[2];
		D[0]=0;
		D[1]=0;
		
			double gamma = B.getBunch().getGamma();
			double N     = A.getBunch().getData("NumParInBun");
			double rad   = B.getParticle().getRad(); 
			Vector3d sigma = A.getBunch().getSigma();
			
			D[0] = (2 * N * rad * sigma.z) / ( gamma * sigma.x * (sigma.x + sigma.y) ); 
			D[1] = (2 * N * rad * sigma.z) / ( gamma * sigma.y * (sigma.x + sigma.y) ); 
			

		
		return D;
	
	}
	
	
	public double [] getBeamBeam(Beam A, Beam B) {
		double []  BB = new double[2];
		BB[0]=0;
		BB[1]=0;


			double gamma = B.getBunch().getGamma();
			double N     = A.getBunch().getData("NumParInBun");
			double rad   = B.getParticle().getRad(); 
			double betaStarVer = B.getBunch().getData("BetaVer");
			double betaStarHor = B.getBunch().getData("BetaHor");
			double betaStar = Math.sqrt(betaStarHor * betaStarVer);
			Vector3d sigma = A.getBunch().getSigma();
			
			
			
			BB[0] = (N * rad *  betaStar) / (2 * Math.PI * gamma * sigma.x * (sigma.x+sigma.y));
			BB[1] = (N * rad *  betaStar) / (2 * Math.PI * gamma * sigma.y * (sigma.x+sigma.y));
			
		
		return BB;
	
	}
	
	public void initParticleInterface()
	{
		particles = new HashMap<String, Particle>();
		particles.put("positron", new Particle("positron", 1.0, 0.510998950e6, 2.82e-15));
		particles.put("electron", new Particle("electron", -1.0, 0.510998950e6, 2.82e-15));
		particles.put("proton", new Particle("proton", 1.0, 938.28e6, 1.54e-18));
		particles.put("muon", new Particle("muon", -1.0, 105.658e6, 1.37e-17 ,2.197e-6));
		particles.put("muon+", new Particle("muon+", 1.0, 105.658e6, 1.37e-17 ,2.197e-6));
		particles.put("lead", new Particle("lead", +82, 1.93e11,(87*87/207.2) * 2.82e-15  ));
	}
	public HashMap<String, Particle> getParticles() {
		return particles;
	}
	
	public void setParticles(HashMap<String, Particle> particles) {
		this.particles = particles;
	}

	public  void locateParticle3d()
	{
		Vector3d sigmaL = beamLeft.getBunch().getSigma();
		Vector3d sigmaR = beamRight.getBunch().getSigma();
		
		
		Vector3d sigmaLs = new Vector3d(0,0,0);
		Vector3d sigmaRs = new Vector3d(0,0,0);
			
		Vector2d betaL = beamLeft.getBunch().getBeta();
		Vector2d betaR = beamRight.getBunch().getBeta();
		
		Vector3d posL = beamLeft.getBunch().getPosition();
		Vector3d posR = beamRight.getBunch().getPosition();
		
	    
	
		
		Vector3d speedL = beamLeft.getBunch().getVelocity();
		
		Vector3d speedR = beamRight.getBunch().getVelocity();
		
		//Hourglass effect
		//https://accelconf.web.cern.ch/p01/PAPERS/RPPH141.PDF
		//https://cds.cern.ch/record/941318/files/p361.pdf  370 syf
		
		sigmaLs.x = Math.pow(sigmaL.x*sigmaL.x*(1+(Math.pow(posL.z, 2))/(betaL.x*betaL.x)), 0.5);
		double sigmaLix = Math.pow(sigmaL.x*sigmaL.x*(1+(Math.pow(posL.z-speedL.z*dTime, 2))/(betaL.x*betaL.x)), 0.5);

		sigmaLs.y = Math.pow(sigmaL.y*sigmaL.y*(1+(Math.pow(posL.z, 2))/(betaL.y*betaL.y)), 0.5);
		double sigmaLiy = Math.pow(sigmaL.y*sigmaL.y*(1+(Math.pow(posL.z-speedL.z*dTime, 2))/(betaL.y*betaL.y)), 0.5);
		sigmaLs.z = sigmaL.z;
		
		double sigmaRix = Math.pow(sigmaR.x*sigmaR.x*(1+(Math.pow(posR.z+speedR.z*dTime, 2))/(betaR.x*betaR.x)), 0.5);
		sigmaRs.x = Math.pow(sigmaR.x*sigmaR.x*(1+(Math.pow(posR.z, 2))/(betaR.x*betaR.x)), 0.5);

		double sigmaRiy = Math.pow(sigmaR.y*sigmaR.y*(1+(Math.pow(posR.z+speedR.z*dTime, 2))/(betaR.y*betaR.y)), 0.5);
		sigmaRs.y = Math.pow(sigmaR.y*sigmaR.y*(1+(Math.pow(posR.z, 2))/(betaR.y*betaR.y)), 0.5);
		sigmaRs.z = sigmaR.z;
		
		
		double NL = beamLeft.getData("NumParInBun");
		double NR = beamRight.getData("NumParInBun");
		
		double qL = beamLeft.getParticle().getCharge();
		double qR = beamRight.getParticle().getCharge();
		double qL_SI = qL *1.602e-19;//  to Coulomb
		double qR_SI = qR *1.602e-19;//  to Coulomb
					
		double QL = NL * qL;
		double QR = NR * qR;
		double QL_SI = QL *1.602e-19;//  to Coulomb
		double QR_SI = QR *1.602e-19;//  to Coulomb
		
		
		double gammaL = beamLeft.getBunch().getGamma();
		double gammaR = beamRight.getBunch().getGamma();
		
		double mL = beamLeft.getParticle().getMass() * 1.782661921e-36;  // eV to kg
		double mR = beamRight.getParticle().getMass() * 1.782661921e-36;  // eV to kg
		
		double toplamMass=0;
		
		int i = 0;
		for (int z = 0; z < RES_g; z++) 
		{
			for (int y = 0; y < RES_g; y++) 
			{
				for (int x = 0; x < RES_g; x++) 
				{
					
					if (presentTime == 0)
					{
				
						double particlesLeftiX = posL.x - (Scale *  sigmaLix) + (x * (2 * Scale * sigmaLix ) / RES_g);
						double particlesLeftiY = posL.y - (Scale *  sigmaLiy) + (y * (2 * Scale * sigmaLiy ) / RES_g);
						double particlesRightiX = posR.x - (Scale *  sigmaRix) + (x * (2 * Scale * sigmaRix ) / RES_g);
						double particlesRightiY = posR.y - (Scale *  sigmaRiy) + (y * (2 * Scale * sigmaRiy ) / RES_g);
						
						boolean hourglassEffect = false;
						if(hourglassEffect)
						{
							
							particlesLeftX[i] = posL.x - (Scale *  sigmaLs.x) + (x * (2 * Scale * sigmaLs.x ) / RES_g);
							particlesLeftY[i] = posL.y - (Scale *  sigmaLs.y) + (y * (2 * Scale * sigmaLs.y ) / RES_g);
							particlesLeftZ[i] = posL.z - (Scale *  sigmaL.z) + (z * (2 * Scale * sigmaL.z ) / RES_g);
							particlesRightX[i] = posR.x - (Scale *  sigmaRs.x) + (x * (2 * Scale * sigmaRs.x ) / RES_g);
							particlesRightY[i] = posR.y - (Scale *  sigmaRs.y) + (y * (2 * Scale * sigmaRs.y ) / RES_g);
							particlesRightZ[i] = posR.z - (Scale *  sigmaR.z) + (z * (2 * Scale * sigmaR.z ) / RES_g);
						}
						else
						{
							particlesLeftX[i] = posL.x - (Scale *  sigmaL.x) + (x * (2 * Scale * sigmaL.x ) / RES_g);
							particlesLeftY[i] = posL.y - (Scale *  sigmaL.y) + (y * (2 * Scale * sigmaL.y ) / RES_g);
							particlesLeftZ[i] = posL.z - (Scale *  sigmaL.z) + (z * (2 * Scale * sigmaL.z ) / RES_g);
							particlesRightX[i] = posR.x - (Scale *  sigmaR.x) + (x * (2 * Scale * sigmaR.x ) / RES_g);
							particlesRightY[i] = posR.y - (Scale *  sigmaR.y) + (y * (2 * Scale * sigmaR.y ) / RES_g);
							particlesRightZ[i] = posR.z - (Scale *  sigmaR.z) + (z * (2 * Scale * sigmaR.z ) / RES_g);
							
						
							
						}
						
						if(hourglassEffect)
						{
							particlesLeftVX[i] = (particlesLeftX[i] - particlesLeftiX)/dTime;
							particlesLeftVY[i] = (particlesLeftY[i] - particlesLeftiY)/dTime;
							particlesRightVX[i] = (particlesRightX[i] - particlesRightiX)/dTime;
							particlesRightVY[i] = (particlesRightY[i] - particlesRightiY)/dTime;
							
						}
						else
						{
							particlesLeftVX[i] = 0;
							particlesLeftVY[i] = 0;
							particlesRightVX[i] = 0;
							particlesRightVY[i] = 0;
						}
						particlesLeftVZ[i] = speedL.z;
						particlesRightVZ[i] = speedR.z;
						
						Vector3d positionL = new Vector3d(particlesLeftX[i]-posL.x,particlesLeftY[i]-posL.y,particlesLeftZ[i]-posL.z);
						Vector3d positionR = new Vector3d(particlesRightX[i]-posR.x,particlesRightY[i]-posR.y,particlesRightZ[i]-posR.z);
						
						if(hourglassEffect)
						{

							gaussMassL[i] = calculateMass(NL, sigmaLs, positionL, mL );
							gaussMassR[i] = calculateMass(NR, sigmaRs, positionR, mR );
								
						}
						else
						{
							gaussMassL[i] = calculateMass(NL, sigmaL, positionL, mL );											
							gaussMassR[i] = calculateMass(NR, sigmaR, positionR, mR );
								
						}
						
						
						particlesLeftPX[i] = gammaL * gaussMassL[i] * particlesLeftVX[i];
						particlesLeftPY[i] = gammaL * gaussMassL[i] * particlesLeftVY[i];
						particlesLeftPZ[i] = gammaL * gaussMassL[i] * particlesLeftVZ[i];

						particlesRightPX[i] = gammaR * gaussMassR[i] * particlesRightVX[i];
						particlesRightPY[i] = gammaR * gaussMassR[i] * particlesRightVY[i];
						particlesRightPZ[i] = gammaR * gaussMassR[i] * particlesRightVZ[i];
					    
						
		
						
					}
					else
					{
						// https://uspas.fnal.gov/materials/19Knoxville/lec%202.pdf
						// Siegmann, Magnetisma
						// https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node153.html
				
						
						
						double xLR = (particlesRightX[i] - posL.x);							
						double yLR = (particlesRightY[i] - posL.y);
						double zLR = (particlesRightZ[i] - posL.z);
						double dLR = Math.sqrt((xLR * xLR) + (yLR * yLR) + (zLR * zLR));
						
						
						double xRL = (particlesLeftX[i] - posR.x);							
						double yRL = (particlesLeftY[i] - posR.y);
						double zRL = (particlesLeftZ[i] - posR.z);
						double dRL = Math.sqrt((xRL * xRL) + (yRL * yRL) + (zRL * zRL));
						
						// SOLDAKI BUNCH STATIK OLDUGU REFERANS SİSTEMİ 

						// soldaki bunch merkezinin sagdan gelen gaussMass'lere olan uzakliklari (;// uzunluk buzulmeleri ile)

					    
						
						double xL = xLR  * Math.sqrt(gammaL*gammaL*Math.abs(xLR)/dLR + Math.sqrt(1 - Math.pow(Math.abs(xLR)/dLR,2) ) ); 									
						double yL = yLR  * Math.sqrt(gammaL*gammaL*Math.abs(yLR)/dLR + Math.sqrt(1 - Math.pow(Math.abs(yLR)/dLR,2) ) );
						double zL = zLR  * Math.sqrt(gammaL*gammaL*Math.abs(zLR)/dLR + Math.sqrt(1 - Math.pow(Math.abs(zLR)/dLR,2) ) );
////////////////////////////////
						double xR = xRL  * Math.sqrt(gammaR*gammaR*Math.abs(xRL)/dRL + Math.sqrt(1 - Math.pow(Math.abs(xRL)/dRL,2) ) ); 									
						double yR = yRL  * Math.sqrt(gammaR*gammaR*Math.abs(yRL)/dRL + Math.sqrt(1 - Math.pow(Math.abs(yRL)/dRL,2) ) );
						double zR = zRL  * Math.sqrt(gammaR*gammaR*Math.abs(zRL)/dRL + Math.sqrt(1 - Math.pow(Math.abs(zRL)/dRL,2) ) );

						
						double dL = Math.sqrt((xL * xL) + (yL * yL) + (zL * zL));
////////////////////////////////
						double dR = Math.sqrt((xR * xR) + (yR * yR) + (zR * zR));
						
						
	
						//soldaki bunch'ın durgun oldugu referans sistemindeki elektrik alan hesabi

						double ELeftX_L =0;					
						double ELeftY_L =0;
						double ELeftZ_L =0;
						double BLeftX_L=0;
						double BLeftY_L=0;
						double BLeftZ_L=0;
						
						
						double ERightX_R =0;  
						double ERightY_R =0;  
						double ERightZ_R =0;  
						double BRightX_R=0;
						double BRightY_R=0;
						double BRightZ_R=0;
						
						
						
						for(int p=0; p < RES_g3 ; p++) {
							
							xLR = (particlesRightX[i] - particlesLeftX[p]);							
							yLR = (particlesRightY[i] - particlesLeftY[p]);
							zLR = (particlesRightZ[i] - particlesLeftZ[p]);
							dLR = Math.sqrt((xLR * xLR) + (yLR * yLR) + (zLR * zLR));
							
						
							xL = xLR  * Math.sqrt(gammaL*gammaL*Math.abs(xLR)/dLR + Math.sqrt(1 - Math.pow(Math.abs(xLR)/dLR,2) ) ); 									
							yL = yLR  * Math.sqrt(gammaL*gammaL*Math.abs(yLR)/dLR + Math.sqrt(1 - Math.pow(Math.abs(yLR)/dLR,2) ) );
							zL = zLR  * Math.sqrt(gammaL*gammaL*Math.abs(zLR)/dLR + Math.sqrt(1 - Math.pow(Math.abs(zLR)/dLR,2) ) );

							
							ELeftX_L += Math.pow(Math.E,-1/(dL*dL)) *  (xL/dL) * k*(qL_SI*gaussMassL[p]/mL)/(dL*dL) ;	
							ELeftY_L += Math.pow(Math.E,-1/(dL*dL)) *  (yL/dL) * k*(qL_SI*gaussMassL[p]/mL)/(dL*dL) ;
							ELeftZ_L += Math.pow(Math.E,-1/(dL*dL)) *  (zL/dL) * k*(qL_SI*gaussMassL[p]/mL)/(dL*dL) ;

							
							BLeftX_L=0;
							BLeftY_L=0;
							BLeftZ_L=0;
							
							
							////////////////////////
							
							
							xRL = (particlesLeftX[i] - particlesRightX[p]);							
							yRL = (particlesLeftY[i] - particlesRightY[p]);
							zRL = (particlesLeftZ[i] - particlesRightZ[p]);
							dRL = Math.sqrt((xRL * xRL) + (yRL * yRL) + (zRL * zRL));
							
							xR = xRL  * Math.sqrt(gammaR*gammaR*Math.abs(xRL)/dRL + Math.sqrt(1 - Math.pow(Math.abs(xRL)/dRL,2) ) ); 									
							yR = yRL  * Math.sqrt(gammaR*gammaR*Math.abs(yRL)/dRL + Math.sqrt(1 - Math.pow(Math.abs(yRL)/dRL,2) ) );
							zR = zRL  * Math.sqrt(gammaR*gammaR*Math.abs(zRL)/dRL + Math.sqrt(1 - Math.pow(Math.abs(zRL)/dRL,2) ) );

							ERightX_R += Math.pow(Math.E,-1/(dR*dR)) * (xR/dR) * k*(qR_SI*gaussMassR[p]/mR)/(dR*dR) ;
							ERightY_R += Math.pow(Math.E,-1/(dR*dR)) * (yR/dR) * k*(qR_SI*gaussMassR[p]/mR)/(dR*dR) ;
							ERightZ_R += Math.pow(Math.E,-1/(dR*dR)) * (zR/dR) * k*(qR_SI*gaussMassR[p]/mR)/(dR*dR) ;

							BRightX_R=0;
							BRightY_R=0;
							BRightZ_R=0;
							
						}
					

						Fields_L[0][0]=0;
						Fields_L[0][1]=-ELeftX_L/c;
						Fields_L[0][2]=-ELeftY_L/c;
						Fields_L[0][3]=-ELeftZ_L/c;
						
						Fields_L[1][0]=-Fields_L[0][1];
						Fields_L[1][1]=0;
						Fields_L[1][2]=-BLeftZ_L;
						Fields_L[1][3]=+BLeftY_L;
						
						Fields_L[2][0]=-Fields_L[0][2];
						Fields_L[2][1]=-Fields_L[1][2];
						Fields_L[2][2]=0;
						Fields_L[2][3]=-BLeftX_L;
						
						Fields_L[3][0]=-Fields_L[0][3];
						Fields_L[3][1]=-Fields_L[1][3];
						Fields_L[3][2]=-Fields_L[2][3];;
						Fields_L[3][3]=0;
						
						
////////////////////////////////
						

						Fields_R[0][0]=0;
						Fields_R[0][1]=-ERightX_R/c;
						Fields_R[0][2]=-ERightY_R/c;
						Fields_R[0][3]=-ERightZ_R/c;
						
						Fields_R[1][0]=-Fields_R[0][1];
						Fields_R[1][1]=0;
						Fields_R[1][2]=-BRightZ_R;
						Fields_R[1][3]=+BRightY_R;
						
						Fields_R[2][0]=-Fields_R[0][2];
						Fields_R[2][1]=-Fields_R[1][2];
						Fields_R[2][2]=0;
						Fields_R[2][3]=-BRightX_R;
						
						Fields_R[3][0]=-Fields_R[0][3];
						Fields_R[3][1]=-Fields_R[1][3];
						Fields_R[3][2]=-Fields_R[2][3];;
						Fields_R[3][3]=0;
						

						
						
						
						// LABFRAME'E KOORDINAT DONUSUMU SONRASI ELEKTRIK ALANLAR
						
						
						Lambda_inv_L[0][0]=gammaL;           
						Lambda_inv_L[0][1]=gammaL*speedL.x/c; 
						Lambda_inv_L[0][2]=gammaL*speedL.y/c; 
						Lambda_inv_L[0][3]=gammaL*speedL.z/c;
						
						Lambda_inv_L[1][0]=Lambda_inv_L[0][1];
						Lambda_inv_L[1][1]=1 + (gammaL-1)* Math.pow(speedL.x/speedL.diagonal(),2);
						Lambda_inv_L[1][2]=0 + (gammaL-1)* speedL.x*speedL.y/ Math.pow(speedL.diagonal(),2);
						Lambda_inv_L[1][3]=0 + (gammaL-1)* speedL.x*speedL.z/ Math.pow(speedL.diagonal(),2);
						
						Lambda_inv_L[2][0]=Lambda_inv_L[0][2];
						Lambda_inv_L[2][1]=Lambda_inv_L[1][2];
						Lambda_inv_L[2][2]=1 + (gammaL-1)* Math.pow(speedL.y/speedL.diagonal(),2);
						Lambda_inv_L[2][3]=0 + (gammaL-1)* speedL.y*speedL.z/ Math.pow(speedL.diagonal(),2);
						
						Lambda_inv_L[3][0]=Lambda_inv_L[0][3];
						Lambda_inv_L[3][1]=Lambda_inv_L[1][3];
						Lambda_inv_L[3][2]=Lambda_inv_L[2][3];
						Lambda_inv_L[3][3]=1 + (gammaL-1)* Math.pow(speedL.z/speedL.diagonal(),2);
						
////////////////////////////////
						
						Lambda_inv_R[0][0]=gammaR;           
						Lambda_inv_R[0][1]=gammaR*speedR.x/c; 
						Lambda_inv_R[0][2]=gammaR*speedR.y/c; 
						Lambda_inv_R[0][3]=gammaR*speedR.z/c;
						
						Lambda_inv_R[1][0]=Lambda_inv_R[0][1];
						Lambda_inv_R[1][1]=1 + (gammaR-1)* Math.pow(speedR.x/speedR.diagonal(),2);
						Lambda_inv_R[1][2]=0 + (gammaR-1)* speedR.x*speedR.y/ Math.pow(speedR.diagonal(),2);
						Lambda_inv_R[1][3]=0 + (gammaR-1)* speedR.x*speedR.z/ Math.pow(speedR.diagonal(),2);
						
						Lambda_inv_R[2][0]=Lambda_inv_R[0][2];
						Lambda_inv_R[2][1]=Lambda_inv_R[1][2];
						Lambda_inv_R[2][2]=1 + (gammaR-1)* Math.pow(speedR.y/speedR.diagonal(),2);
						Lambda_inv_R[2][3]=0 + (gammaR-1)* speedR.y*speedR.z/ Math.pow(speedR.diagonal(),2);
						
						Lambda_inv_R[3][0]=Lambda_inv_R[0][3];
						Lambda_inv_R[3][1]=Lambda_inv_R[1][3];
						Lambda_inv_R[3][2]=Lambda_inv_R[2][3];
						Lambda_inv_R[3][3]=1 + (gammaR-1)* Math.pow(speedR.z/speedR.diagonal(),2);
						
						
						double ELeftX = 0;
						double ELeftY = 0;
						double ELeftZ = 0;
						double BLeftX = 0;
						double BLeftY = 0;
						double BLeftZ = 0;
////////////////////////////////
						double ERightX = 0;
						double ERightY = 0;
						double ERightZ = 0;
						double BRightX = 0;
						double BRightY = 0;
						double BRightZ = 0;
						
						for (int mu = 0; mu<=3; mu++) {
							for (int nu = 0; nu<=3; nu++) {
								
								ELeftX +=Lambda_inv_L[1][nu] * Lambda_inv_L[0][mu] * Fields_L[mu][nu];
								ELeftY +=Lambda_inv_L[2][nu] * Lambda_inv_L[0][mu] * Fields_L[mu][nu];
								ELeftZ +=Lambda_inv_L[3][nu] * Lambda_inv_L[0][mu] * Fields_L[mu][nu];
								
								BLeftX +=Lambda_inv_L[2][nu] * Lambda_inv_L[3][mu] * Fields_L[mu][nu];
								BLeftY +=Lambda_inv_L[1][nu] * Lambda_inv_L[3][mu] * Fields_L[mu][nu];
								BLeftZ +=Lambda_inv_L[1][nu] * Lambda_inv_L[2][mu] * Fields_L[mu][nu];
////////////////////////////////							
								
								ERightX +=Lambda_inv_R[1][nu] * Lambda_inv_R[0][mu] * Fields_R[mu][nu];
								ERightY +=Lambda_inv_R[2][nu] * Lambda_inv_R[0][mu] * Fields_R[mu][nu];
								ERightZ +=Lambda_inv_R[3][nu] * Lambda_inv_R[0][mu] * Fields_R[mu][nu];
								
								BRightX +=Lambda_inv_R[2][nu] * Lambda_inv_R[3][mu] * Fields_R[mu][nu];
								BRightY +=Lambda_inv_R[1][nu] * Lambda_inv_R[3][mu] * Fields_R[mu][nu];
								BRightZ +=Lambda_inv_R[1][nu] * Lambda_inv_R[2][mu] * Fields_R[mu][nu];
								
							}	
						}
						
						ELeftX *= -c;
						ELeftY *= -c;
						ELeftZ *= -c;
						BLeftX *= -1;
						BLeftY *= +1;
						BLeftZ *= -1;
////////////////////////////////						
						ERightX *= -c;
						ERightY *= -c;
						ERightZ *= -c;
						BRightX *= -1;
						BRightY *= +1;
						BRightZ *= -1;
						
						
						
						
						
						Fields_L_lab[0][0]=0;
						Fields_L_lab[0][1]=-ELeftX/c;
						Fields_L_lab[0][2]=-ELeftY/c;
						Fields_L_lab[0][3]=-ELeftZ/c;
						
						Fields_L_lab[1][0]=-Fields_L_lab[0][1];
						Fields_L_lab[1][1]=0;
						Fields_L_lab[1][2]=-BLeftZ;
						Fields_L_lab[1][3]=+BLeftY;
						
						Fields_L_lab[2][0]=-Fields_L_lab[0][2];
						Fields_L_lab[2][1]=-Fields_L_lab[1][2];
						Fields_L_lab[2][2]=0;
						Fields_L_lab[2][3]=-BLeftX;
						
						Fields_L_lab[3][0]=-Fields_L_lab[0][3];
						Fields_L_lab[3][1]=-Fields_L_lab[1][3];
						Fields_L_lab[3][2]=-Fields_L_lab[2][3];;
						Fields_L_lab[3][3]=0;
////////////////////////////////	
						Fields_R_lab[0][0]=0;
						Fields_R_lab[0][1]=-ERightX/c;
						Fields_R_lab[0][2]=-ERightY/c;
						Fields_R_lab[0][3]=-ERightZ/c;
						
						Fields_R_lab[1][0]=-Fields_R_lab[0][1];
						Fields_R_lab[1][1]=0;
						Fields_R_lab[1][2]=-BRightZ;
						Fields_R_lab[1][3]=+BRightY;
						
						Fields_R_lab[2][0]=-Fields_R_lab[0][2];
						Fields_R_lab[2][1]=-Fields_R_lab[1][2];
						Fields_R_lab[2][2]=0;
						Fields_R_lab[2][3]=-BRightX;
						
						Fields_R_lab[3][0]=-Fields_R_lab[0][3];
						Fields_R_lab[3][1]=-Fields_R_lab[1][3];
						Fields_R_lab[3][2]=-Fields_R_lab[2][3];;
						Fields_R_lab[3][3]=0;
						
						//ARTIK ELIMIZDE LABFRAME FILEDS VAR SIRA LORENTZ KUVVETINDE
						
						URight[i][0] = c; 
						URight[i][1] = -particlesRightVX[i]; 
						URight[i][2] = -particlesRightVY[i]; 
						URight[i][3] = -particlesRightVZ[i]; 
						
////////////////////////////////						
						ULeft[i][0] = c; 
						ULeft[i][1] = -particlesLeftVX[i]; 
						ULeft[i][2] = -particlesLeftVY[i]; 
						ULeft[i][3] = -particlesLeftVZ[i]; 
						

							for (int beta = 0; beta<=3; beta++) { //Lorentz force in 4 dimensions, no change in zero^th dimension :  c
						
						particlesRightPX[i] +=  (qR_SI*gaussMassR[i]/mR) *    Fields_L_lab[1][beta] * URight[i][beta] / gammaR; // Eger beamstrahlung dikkate alinirsa gammaL ve gammaR sabit olmamali					
						particlesRightPY[i] +=  (qR_SI*gaussMassR[i]/mR) *    Fields_L_lab[2][beta] * URight[i][beta] / gammaR;
						particlesRightPZ[i] +=  (qR_SI*gaussMassR[i]/mR) *    Fields_L_lab[3][beta] * URight[i][beta] / gammaR;
////////////////////////////////						
						particlesLeftPX[i] +=  (qL_SI*gaussMassL[i]/mL) *    Fields_R_lab[1][beta] * ULeft[i][beta] / gammaL; // Eger beamstrahlung dikkate alinirsa gammaL ve gammaR sabit olmamali					
						particlesLeftPY[i] +=  (qL_SI*gaussMassL[i]/mL) *    Fields_R_lab[2][beta] * ULeft[i][beta] / gammaL;
						particlesLeftPZ[i] +=  (qL_SI*gaussMassL[i]/mL) *    Fields_R_lab[3][beta] * ULeft[i][beta] / gammaL;
							
															}
						
					    particlesRightVX[i] = gaussMassR[i] == 0 ? Vector3d.MaxValue.x : particlesRightPX[i] / (gammaR*gaussMassR[i]);
					    particlesRightVY[i] = gaussMassR[i] == 0 ? Vector3d.MaxValue.y : particlesRightPY[i] / (gammaR*gaussMassR[i]);
					    particlesRightVZ[i] = gaussMassR[i] == 0 ? Vector3d.MaxValue.z : particlesRightPZ[i] / (gammaR*gaussMassR[i]);
////////////////////////////////					    
					    particlesLeftVX[i] = gaussMassL[i] == 0 ? Vector3d.MaxValue.x : particlesLeftPX[i] / (gammaL*gaussMassL[i]);
					    particlesLeftVY[i] = gaussMassL[i] == 0 ? Vector3d.MaxValue.y : particlesLeftPY[i] / (gammaL*gaussMassL[i]);
					    particlesLeftVZ[i] = gaussMassL[i] == 0 ? Vector3d.MaxValue.z : particlesLeftPZ[i] / (gammaL*gaussMassL[i]);
							

			
					        particlesLeftX[i] +=particlesLeftVX[i]*dTime;
						    particlesLeftY[i] +=particlesLeftVY[i]*dTime;
						    particlesLeftZ[i] +=c*dTime;//1*Math.abs(particlesLeftVZ[i])
									    
////////////////////////////////


							
					        particlesRightX[i] +=particlesRightVX[i]*dTime;
						    particlesRightY[i] +=particlesRightVY[i]*dTime;
						    particlesRightZ[i] -=c*dTime;//1*Math.abs(particlesRightVZ[i])

												
						    
		if(i==44) {
System.out.println(presentTime / totalTime);
System.out.println("Vx: "+ particlesRightVX[i]);
System.out.println("Vy: "+ particlesRightVY[i]);
System.out.println("X: "+ particlesRightX[i]);
System.out.println("Y: "+ particlesRightY[i]);
System.out.println("ELeftX: "+ ELeftX);
System.out.println("--------------------------------");
System.out.println("VxL: "+ particlesLeftVX[i]);
System.out.println("VyL: "+ particlesLeftVY[i]);
System.out.println("XL: "+ particlesLeftX[i]);
System.out.println("YL: "+ particlesLeftY[i]);
System.out.println("ERightX: "+ ERightX);
System.out.println("--------------------------------");
		}
						
					}
					i++;
					
					
					
				}
			}

	    }

	}




	public void spaceIntegration2()
	{
		
		Vector3d sigmaL = beamLeft.getBunch().getSigma();
		Vector3d sigmaR = beamRight.getBunch().getSigma();
		Vector3d posL = beamLeft.getBunch().getPosition();
		Vector3d posR = beamRight.getBunch().getPosition();

		Bunch bminX = sigmaL.x <= sigmaR.x ? beamLeft.getBunch() : beamRight.getBunch();
		Bunch bminY = sigmaL.y <= sigmaR.y ? beamLeft.getBunch() : beamRight.getBunch();
		Bunch bminZ = sigmaL.z <= sigmaR.z ? beamLeft.getBunch() : beamRight.getBunch();
		
		intVolume.updateCenter(bminX.getPosition().x,bminY.getPosition().y , bminZ.getPosition().z);
		intVolume.updateSize(bminX.getSigma().x*IntScale,bminY.getSigma().y*IntScale , bminZ.getSigma().z*IntScale);
		intVolume.updateDq(RES_q);
		
		double phi = alohep.getDataManager().getSettingsData().get("Angle");
		double sin = Math.sin(phi/2);
		double cos = Math.cos(phi/2);
		
		dLxyz = 0;
		double zi = intVolume.getCenter().z - intVolume.getSize().z;
		double yi = intVolume.getCenter().y - intVolume.getSize().y;
		double xi = intVolume.getCenter().x - intVolume.getSize().x;
		double zf = intVolume.getCenter().z + intVolume.getSize().z;
		double yf = intVolume.getCenter().y + intVolume.getSize().y;
		double xf = intVolume.getCenter().x + intVolume.getSize().x;
		double dz = intVolume.getDq().z;
		double dy = intVolume.getDq().y;
		double dx = intVolume.getDq().x;
		

		double dLxyzp = getRo(particlesLeftZ, sigmaL.z, zi, posL.z)  * getRo(particlesRightZ, sigmaR.z, zi, posR.z);
		
		for(double z = zi; z <= zf; z+=dz)
		{
			double dLxy = 0;
			double dLxyp = getRo(particlesLeftY, sigmaL.y, yi, posL.y )  * getRo(particlesRightY, sigmaR.y, yi, posR.y );

			for(double y = yi; y <= yf; y+=dy)
			{
				double dLx = 0;
				double dLxp = getRo(particlesLeftX, sigmaL.x, xi, posL.x ) * getRo(particlesRightX, sigmaR.x, xi, posR.x );
			
				for(double x = xi+dx; x <= xf; x+=dx)
				{
					double dLxn = getRo(particlesLeftX, sigmaL.x, x, posL.x ) * getRo(particlesRightX, sigmaR.x, x, posR.x );
					dLx += (dLxp+dLxn);
					dLxp = dLxn;
				}
				
				double dLxyn = getRo(particlesLeftY, sigmaL.y, y, posL.y ) * getRo(particlesRightY, sigmaR.y, y , posR.y );
				dLxy += (dLxyp+dLxyn)*dLx;
				dLxyp = dLxyn;
			}
			double dLxyzn = getRo(particlesLeftZ, sigmaL.z, z, posL.z ) * getRo(particlesRightZ, sigmaR.z, z, posR.z  );
			dLxyz += (dLxyzp+dLxyzn)*dLxy;
			dLxyzp = dLxyzn;
		}
		dLxyz *= dx * dy * dz / 8;
	}

	
	public double getRo(double[] particleArray, double sigma, double q, double pos) // burada q = x,y, veya z olacak
	{
		double ro = 0;
		for(int i=0; i<particleArray.length; i++) 
		{
			double q1 = q - pos;
			double q2 = q - particleArray[i];

			double res = (2 * gaussScale * sigma ) / (RES_g);

			ro +=  Math.sqrt(1/(2*Math.PI)) * Math.pow(Math.E, -( q1 * q1) / (2*sigma*sigma))   *  Math.pow(Math.E, -(q2 * q2) / (2*res*res)) / (sigma*RES_g*RES_g) ; // sqrt(2*pi)*sigma nÄ±n dÄ±ÅŸÄ±nda Res_g*Res_g'ye boluyoruz cunku 3 boyuttan dolayi bi o kadar parcacik daha var el aldigimiz
		}

		return ro;
	}
	
	public double calculateMass(double N, Vector3d sigma, Vector3d pos, double massOfParticle) {
		double resx = (2 * gaussScale * sigma.x ) / (RES_g);
		double resy = (2 * gaussScale * sigma.y ) / (RES_g);
		double resz = (2 * gaussScale * sigma.z ) / (RES_g);
		
		double mass = resx * Math.pow(Math.E, -( pos.x * pos.x) / (2*sigma.x * sigma.x))    / (sigma.x);
		       mass*=resy * Math.pow(Math.E, -( pos.y * pos.y) / (2*sigma.y * sigma.y))    / (sigma.y);
		       mass*=resz * Math.pow(Math.E, -( pos.z * pos.z) / (2*sigma.z * sigma.z))    / (sigma.z);
		       mass*=(N*massOfParticle/Math.pow(2*Math.PI,1.5));

		return mass;
	}
	public double sum (double[] array) {
		
		double sonuc =0;
		
		for(int i=0; i<array.length; i++) {
			sonuc+=array[i];
		}
		
		return sonuc;
	}
	
	
	public Vector3d vectorMultiply(Vector3d A, Vector3d B) // (sonuc = A x B)
	{
		Vector3d sonuc = new Vector3d(0,0,0);
		sonuc.x = A.y * B.z - A.z * B.y;
		sonuc.y = A.z * B.x - A.x * B.z;
		sonuc.z = A.x * B.y - A.y * B.x;
		return sonuc;
		
	}
	public double vScalarMultiply(Vector3d A, Vector3d B) // (sonuc = A . B)
	{
		double sonuc = 0;
	
		sonuc = A.x * B.x + A.y * B.y + A.z * B.z;
		
		return sonuc;
		
	}
	public Vector3d sScalarMultiply(double k, Vector3d A) // (sonuc = k . A)
	{
		Vector3d sonuc = new Vector3d(0,0,0);
	
		sonuc.x = k*A.x;
		sonuc.y = k*A.y;
		sonuc.z = k*A.z;
		
		return sonuc;
		
	}
	
	public Vector3d vectorAdd(Vector3d A, Vector3d B) // (sonuc = A + B)
	{
		Vector3d sonuc = new Vector3d(0,0,0);
	
		sonuc.x = A.x + B.x;
		sonuc.y = A.y + B.y;
		sonuc.z = A.z + B.z;
		
		return sonuc;
		
	}
	public Vector3d vectorSubtract(Vector3d A, Vector3d B) // (sonuc = A - B)
	{
		Vector3d sonuc = new Vector3d(0,0,0);
	
		sonuc.x = A.x - B.x;
		sonuc.y = A.y - B.y;
		sonuc.z = A.z - B.z;
		
		return sonuc;
		
	}
	
	public Vector3d vectorDivide(Vector3d A , double k) // (sonuc = A / k)
	{
		Vector3d sonuc = new Vector3d(0,0,0);
	
		sonuc.x = A.x / k;
		sonuc.y = A.y / k;
		sonuc.z = A.z / k;
		
		return sonuc;
		
	}

	public Beam getBeamLeft() {
		return beamLeft;
	}

	public void setBeamLeft(Beam beamLeft) {
		this.beamLeft = beamLeft;
	}

	public Beam getBeamRight() {
		return beamRight;
	}

	public void setBeamRight(Beam beamRight) {
		this.beamRight = beamRight;
	}
	
	

}
