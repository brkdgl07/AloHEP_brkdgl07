package org.hepforge.alohep.calc;

public class SubBunch {

	public Vector3d position;
	private Vector3d velocity;
	private Vector3d momentum;
	
	public static final double c = 2.99e8;
	public static final double eps0 = 8.85419e-12;
	public static final double mu0 = 4 * Math.PI * 1e-7;
	public static final double eQ = 1.602e-19;
	private static final double k = 1 / (4 * Math.PI * eps0);
	private Bunch bunch;
	private double mass;
	private double charge;
	private double gamma;

	public double deltaX;
	public double deltaY;
	
	private double collisionDensity = 0;
	
	
	public SubBunch(Bunch bunch, Vector3d position, Vector3d velocity)
	{
		this.bunch = bunch;
		this.position = position;
		this.velocity = velocity;
	}
	public void init()
	{
		gamma = bunch.getGamma();
		charge = bunch.getChargeOfBunch() / Bunch.Size;
		mass = bunch.getMassOfParticle()*bunch.getNumberOfParticles()/Bunch.Size;
		momentum = velocity.copy().multiply(gamma*mass);
	}
	
	public void initCollision(Bunch bunch)
	{
		deltaX = -getDisruption(bunch)[0]* position.x;
		deltaY = -getDisruption(bunch)[1]* position.y;
	}

	public void update(Bunch otherBunch, double dt)
	{	

		Vector3d EFieldLab = new Vector3d(0,0,0); //DONE
		Vector3d BFieldLab = new Vector3d(0,0,0); //DONE

		
		for(int i = 0;i < Bunch.Size; i++)
		{
			
			SubBunch other = otherBunch.getsubBunches()[i]; //DONE
			/*
			Vector3d distanceVec = position.copy().subtract(other.getPosition()); //DONE

			Vector3d velocity = other.getVelocity(); //DONE
			double speed = velocity.diagonal(); //DONE
			Vector3d transformedDistanceVec = lengthContraction(distanceVec, velocity);

			Vector3d coordScale = new Vector3d(Vector3d.unitVec.copy().divide(bunch.getTransformedSubSigma()));

			Vector3d scaleTransformedDistanceVec = transformedDistanceVec.multiply(coordScale);
			
			double r = scaleTransformedDistanceVec.diagonal(); //DONE
			
			double charge = other.charge*(2*Math.PI*erf(r/Math.sqrt(2))-2*Math.sqrt(2*Math.PI)*Math.exp(-r*r/2)*r);
			//volume integral of gaussian distribution 

			Vector3d EField = transformedDistanceVec.copy().normalize().multiply(k*charge/(r*r));
			Vector3d relativeVelocity = velocity.copy().multiply(-1);
			Vector3d EFieldProjectionVec = relativeVelocity.copy().multiply(EField.scalarMultiply(relativeVelocity)/(speed*speed));
			Vector3d EFieldRejectionVec = EField.copy().subtract(EFieldProjectionVec);
			
			EFieldLab.add(EFieldRejectionVec.multiply(other.getGamma()));
			
			Vector3d normalizedVelocity = relativeVelocity.copy().normalize();
			BFieldLab.add(normalizedVelocity.copy().multiply(-other.getGamma()).vectorMultiply(EField.multiply(1/c))); // eksi isareti burada geliyor yukarida carpmaya gerek var mi?
			
			*/collisionDensity += 1/(bunch.getSubSigma().volume()*otherBunch.getSubSigma().volume())*collisionDensity(otherBunch, other);
		}
/*
		
		Vector3d dP = (EFieldLab.copy().add(velocity.copy().vectorMultiply(BFieldLab)).multiply(1e5*charge*dt));

		momentum.add(dP);*/
	}
	
	public double [] getDisruption(Bunch otherbunch) {
		double []  D = new double[2];
		D[0]=0;
		D[1]=0;
		
			
			double N     = otherbunch.getData("NumParInBun");
			 
			Vector3d sigma = otherbunch.getSigma();
			
			D[0] = (2 * N * bunch.getBeam().getParticle().getRad() * sigma.z) / ( gamma * sigma.x * (sigma.x + sigma.y) );
			D[1] = (2 * N * bunch.getBeam().getParticle().getRad() * sigma.z) / ( gamma * sigma.y * (sigma.x + sigma.y) ); 


		
		return D;
	
	}
	public void updatePosition(double dt, Bunch otherBunch)
	{
		
		/*double P = momentum.diagonal();
		
		double newSpeed = P / Math.sqrt(mass*mass + P*P/(c*c));
		gamma = 1/Math.sqrt(1-newSpeed*newSpeed/c/c);
		this.velocity = momentum.copy().normalize().multiply(newSpeed);*/
		position.add(this.velocity.copy().multiply(dt));
		position.add(new Vector3d(deltaX*c*dt/(2*otherBunch.getSigma().z),deltaY*c*dt/(2*otherBunch.getSigma().z),0));

		System.out.println(position.print());
		
	}

	/*public void update(Bunch otherBunch, double dt)
	{	
		double dpX = 0;
		double dpY = 0;
		for(int i = 0;i < Bunch.Size; i++)
		{
			
			SubBunch other = otherBunch.getsubBunches()[i];
			Vector3d s = otherBunch.getSubSigma();
			double r = position.copy().subtract(other.getPosition()).diagonal();
			dpX += -2*otherBunch.getNumberOfParticles()/Bunch.Size*otherBunch.getRParticle()/(otherBunch.getGamma()*r)*(1-Math.exp(-r*r/(s.x*(s.x+s.y))))
					*(bunch.getGamma()*bunch.getMassOfParticle()*c);
			dpY += -2*otherBunch.getNumberOfParticles()/Bunch.Size*otherBunch.getRParticle()/(otherBunch.getGamma()*r)*(1-Math.exp(-r*r/(s.y*(s.x+s.y))))
					*(bunch.getGamma()*bunch.getMassOfParticle()*c);
			//collisionDensity += 1/(bunch.getSubSigma().volume()*otherBunch.getSubSigma().volume())*collisionDensity(otherBunch, other);
		}
		momentum.add(new Vector3d(dpX, dpY,0));
	}*/

	public static Vector3d lengthContraction(Vector3d lengthVec, Vector3d velocity)
	{
		double speed = velocity.diagonal();
		
		Vector3d projectionVec = velocity.copy().multiply(lengthVec.copy().scalarMultiply(velocity)/(speed*speed)); //DONE
		Vector3d rejectionVec = lengthVec.copy().subtract(projectionVec);

		double projectionVecDiagonal = projectionVec.diagonal();
		double gamma = 1/Math.sqrt(1-speed*speed/c/c);
		
		projectionVec.normalize().multiply(projectionVecDiagonal*gamma);
		Vector3d transformedDistanceVec = projectionVec.copy().add(rejectionVec);
		
		return transformedDistanceVec;
	}
	public double collisionDensity(Bunch otherBunch, SubBunch other)
	{
		double k = 20;
		double t = 20;
		Vector3d distance = position.copy().subtract(other.position);
		if(Math.abs(distance.x) > bunch.getSubSigma().x*t || Math.abs(distance.y) > bunch.getSubSigma().y*t || Math.abs(distance.z) > bunch.getSubSigma().z*t)
			return 0;
		Vector3d center = position.copy().add(distance.multiply(0.5));
		return (integral(center.x+bunch.getSubSigma().x*k, position.x, other.position.x, bunch.getSubSigma().x, otherBunch.getSubSigma().x)
			   -integral(center.x-bunch.getSubSigma().x*k, position.x, other.position.x, bunch.getSubSigma().x, otherBunch.getSubSigma().x))
			  *(integral(center.y+bunch.getSubSigma().y*k, position.y, other.position.y, bunch.getSubSigma().y, otherBunch.getSubSigma().y)
			   -integral(center.y-bunch.getSubSigma().y*k, position.y, other.position.y, bunch.getSubSigma().y, otherBunch.getSubSigma().y))
			  *(integral(center.z+bunch.getSubSigma().z*k, position.z, other.position.z, bunch.getSubSigma().z, otherBunch.getSubSigma().z)
			   -integral(center.z-bunch.getSubSigma().z*k, position.z, other.position.z, bunch.getSubSigma().z, otherBunch.getSubSigma().z));
				
	}
	
	public Vector4d lorentzTransform(double gamma, Vector3d v, Vector4d st)
	{
		double s = v.diagonal();
		double l[][]= {{gamma, -gamma*v.x/c, -gamma*v.y/c, -gamma*v.z/c}
		,{-gamma*v.x/c, 1+(gamma-1)*v.x*v.x/(s*s),(gamma-1)*v.x*v.y/(s*s), (gamma-1)*v.x*v.z/(s*s)}
		,{-gamma*v.y/c, (gamma-1)*v.y*v.x/(s*s), 1+(gamma-1)*v.y*v.y/(s*s), (gamma-1)*v.y*v.z/(s*s)}
		,{-gamma*v.z/c, (gamma-1)*v.z*v.x/(s*s), (gamma-1)*v.z*v.y/(s*s), 1+(gamma-1)*v.z*v.z/(s*s)}};
		
		
		st.t = l[0][0]*st.t + l[0][1]*st.x + l[0][2]*st.y +l[0][3]*st.z;
		st.x = l[1][0]*st.t + l[1][1]*st.x + l[1][2]*st.y +l[1][3]*st.z;
		st.y = l[2][0]*st.t + l[2][1]*st.x + l[2][2]*st.y +l[2][3]*st.z;
		st.z = l[3][0]*st.t + l[3][1]*st.x + l[3][2]*st.y +l[3][3]*st.z;

		return st;
	}
	
	private double errorTerm(double x, double dist1, double dist2, double sigma1, double sigma2)
	{
		return (-dist1*sigma2*sigma2-dist2*sigma1*sigma1+x*(sigma1*sigma1+sigma2*sigma2))/(Math.pow(2*(sigma1*sigma1+sigma2*sigma2), 0.5)*sigma1*sigma2);
	}
	
	private double integral(double x, double dist1, double dist2, double sigma1, double sigma2)
	{
		if(erf(errorTerm(x, dist1, dist2, sigma1, sigma2)) > 1)
			return 0;
		return Math.pow(Math.PI/2, 0.5)*sigma1*sigma2*Math.exp(-Math.pow(dist1-dist2, 2)/(2*(sigma1*sigma1+sigma2*sigma2)))
				*erf(errorTerm(x, dist1, dist2, sigma1, sigma2))/Math.pow(sigma1*sigma1+sigma2*sigma2, 0.5);
		
	}
	private double erf(double x)
	{
		double d = 1-1/Math.pow((1+0.278393*x+0.230389*x*x+0.000972*x*x*x+0.078108*x*x*x*x),4);
		return x > 0 ? d : -d;
	}
/*	private double erf(double x)
	{
		return Math.atan(3*x)/(Math.PI/2);
	}
	private double erf(double x)
	{
		
		double constant = 2 / Math.pow(Math.PI, 0.5);
		double value = 0;
		for(int i = 0; i < 11; i++)
		{
			value += Math.pow(-1, i)*Math.pow(x, 2*i+1)/(2*i+1)/factorial(i);
		}
		return constant*value;
	}
	*/
	
	private int factorial(int i)
	{
		int answer = 1;
		for(int k = 2; k <= i; k++)
			answer*=k;
		return answer;
	}
	public Vector3d getPosition() {
		return position;
	}
	public void setPosition(Vector3d position) {
		this.position = position;
	}
	public Vector3d getVelocity() {
		return velocity;
	}
	public void setVelocity(Vector3d velocity) {
		this.velocity = velocity;
	}
	public Vector3d getMomentum() {
		return momentum;
	}
	public void setMomentum(Vector3d momentum) {
		this.momentum = momentum;
	}
	public double getCollisionDensity() {
		return collisionDensity;
	}
	public void setCollisionDensity(double collisionDensity) {
		this.collisionDensity = collisionDensity;
	}
	public double getMass() {
		return mass;
	}
	public void setMass(double mass) {
		this.mass = mass;
	}
	public double getCharge() {
		return charge;
	}
	public void setCharge(double charge) {
		this.charge = charge;
	}
	public double getGamma() {
		return gamma;
	}
	public void setGamma(double gamma) {
		this.gamma = gamma;
	}
	
}
