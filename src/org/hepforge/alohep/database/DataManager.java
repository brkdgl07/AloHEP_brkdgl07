package org.hepforge.alohep.database;

import java.util.HashMap;
import java.util.LinkedHashMap;

public class DataManager {

	private LinkedHashMap<String, VariableData> variableDataMap;
	private LinkedHashMap <String, Double> settingsData;
	private HashMap<String, ParticleData> particleDataMap;
	
	
	public DataManager()
	{
		variableDataMap = new LinkedHashMap<String, VariableData>();
		settingsData = new LinkedHashMap<String, Double>();
		particleDataMap = new HashMap<String, ParticleData>();
		
		loadVariableData();
		loadSettingsData();
		loadParticlesData();
		
		
	}
	public void loadParticlesData()
	{
		
		AcceleratorData ad;
		ParticleData pd;
		
		pd = new ParticleData();
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2e10);
		ad.put("EnBeam", 250.0);
		ad.put("BetaVer", 4.8e-4);
		ad.put("BetaHor", 1.1e-2);
		ad.put("EmNorVer", 3.5e-8);
		ad.put("EmNorHor", 1.0e-5);
		ad.put("PulFrq", 5.0);
		ad.put("NumBunInBeam", 1312.0);
		ad.put("BunLen", 3e-4);
		ad.put("BunSpace", 554e-9);
		pd.put("ILC", ad);

		ad = new AcceleratorData();
		ad.put("NumParInBun", 1e10);
		ad.put("EnBeam", 5000.0);
		ad.put("BetaVer", 9.9e-5);
		ad.put("BetaHor", 11e-3);
		ad.put("EmNorVer", 3.5e-8);
		ad.put("EmNorHor", 1.0e-5);
		ad.put("PulFrq", 5.0e3);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 2e-5);
		ad.put("BunSpace", 2e-4);
		pd.put("PWFA", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2.3e10);
		ad.put("EnBeam", 1.0);
		ad.put("BetaVer", 0.5e-2);
		ad.put("BetaHor", 0.5e-2);
		ad.put("EmNorVer", 7.57e-6);
		ad.put("EmNorHor", 7.57e-6);
		ad.put("BunLen", 0.5e-2);
		ad.put("ColFrq", 500e6);
		pd.put("TAC", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1e10);
		ad.put("EnBeam", 1000.0);
		ad.put("BetaVer", 1.0);
		ad.put("BetaHor", 1.0);
		ad.put("EmNorVer", 1.0);
		ad.put("EmNorHor", 1.0);
		ad.put("PulFrq", 1.0);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 2.0);
		pd.put("Test", ad);
		
		pd.setLinear();
		particleDataMap.put("electron-linac", pd);
		
		pd = new ParticleData();
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2e10);
		ad.put("EnBeam", 250.0);
		ad.put("BetaVer", 4.8e-4);
		ad.put("BetaHor", 1.1e-2);
		ad.put("EmNorVer", 3.5e-8);
		ad.put("EmNorHor", 1.0e-5);
		ad.put("PulFrq", 5.0);
		ad.put("NumBunInBeam", 1312.0);
		ad.put("BunLen", 3e-4);
		ad.put("BunSpace", 554e-9);
		pd.put("ILC", ad);

		ad = new AcceleratorData();
		ad.put("NumParInBun", 1e10);
		ad.put("EnBeam", 5000.0);
		ad.put("BetaVer", 9.9e-5);
		ad.put("BetaHor", 11e-3);
		ad.put("EmNorVer", 3.5e-8);
		ad.put("EmNorHor", 1.0e-5);
		ad.put("PulFrq", 5.0e3);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 2e-5);
		ad.put("BunSpace", 2e-4);
		pd.put("PWFA", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1e10);
		ad.put("EnBeam", 1000.0);
		ad.put("BetaVer", 1.0);
		ad.put("BetaHor", 1.0);
		ad.put("EmNorVer", 1.0);
		ad.put("EmNorHor", 1.0);
		ad.put("PulFrq", 1.0);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 2.0);
		pd.put("Test", ad);
		
		pd.setLinear();
		particleDataMap.put("positron-linac", pd);
		
		pd = new ParticleData();

		ad = new AcceleratorData();
		ad.put("NumParInBun", 1.8e10);
		ad.put("EnBeam", 3.56);
		ad.put("BetaVer", 0.5e-2);
		ad.put("BetaHor", 0.5e-2);
		ad.put("EmNorVer", 27e-6);
		ad.put("EmNorHor", 27e-6);
		ad.put("BunLen", 0.5e-2);
		ad.put("ColFrq", 500e6);
		pd.put("TAC", ad);
		
		particleDataMap.put("positron-ring", pd);
		
		pd = new ParticleData();
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1.15e11);
		ad.put("EnBeam", 7000.0);
		ad.put("BetaVer", 0.55);
		ad.put("BetaHor", 0.55);
		ad.put("EmNorVer", 3.75e-6);
		ad.put("EmNorHor", 3.75e-6);
		ad.put("RevFrq", 11245.0);
		ad.put("NumBunInBeam", 2808.0);
		ad.put("BunLen", 7.55e-2);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 27.0);
		pd.put("LHC", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1.6e10);
		ad.put("EnBeam", 4000.0);
		ad.put("BetaVer", 0.8);
		ad.put("BetaHor", 0.8);
		ad.put("EmNorVer", 2e-6);
		ad.put("EmNorHor", 2e-6);
		ad.put("RevFrq", 11245.0);
		ad.put("NumBunInBeam", 338.0);
		ad.put("BunLen", 7.55e-2);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 27.0);
		pd.put("LHC-lead", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1e11);
		ad.put("EnBeam", 50000.0);
		ad.put("BetaVer", 1.1);
		ad.put("BetaHor", 1.1);
		ad.put("EmNorVer", 2.2e-6);
		ad.put("EmNorHor", 2.2e-6);
		ad.put("RevFrq", 2998.0);
		ad.put("NumBunInBeam", 10600.0);
		ad.put("BunLen", 0.08);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 100.0);
		pd.put("FCC", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2.2e11);
		ad.put("EnBeam", 50000.0);
		ad.put("BetaVer", 0.1);
		ad.put("BetaHor", 0.1);
		ad.put("EmNorVer", 2.2e-6);
		ad.put("EmNorHor", 2.2e-6);
		ad.put("RevFrq", 2998.0);
		ad.put("NumBunInBeam", 10600.0);
		ad.put("BunLen", 0.08);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 100.0);
		pd.put("FCC upgrade", ad);
		
		particleDataMap.put("proton", pd);
		
		pd = new ParticleData();
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 7e7);
		ad.put("EnBeam", 5.74e5);
		ad.put("BetaVer", 0.5);
		ad.put("BetaHor", 0.5);
		ad.put("EmNorVer", 1.5e-6);
		ad.put("EmNorHor", 1.5e-6);
		ad.put("RevFrq", 11245.0);
		ad.put("NumBunInBeam", 592.0);
		ad.put("BunLen", 7.94e-2);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 27.0);
		pd.put("LHC", ad);
				
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1.4e8);
		ad.put("EnBeam", 4.1e6);
		ad.put("BetaVer", 1.1);
		ad.put("BetaHor", 1.1);
		ad.put("EmNorVer", 1.5e-6);
		ad.put("EmNorHor", 1.5e-6);
		ad.put("RevFrq", 2998.0);
		ad.put("NumBunInBeam", 432.0);
		ad.put("BunLen", 0.08);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 100.0);
		pd.put("FCC", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 1.2e8);
		ad.put("EnBeam", 328000.0);
		ad.put("BetaVer", 0.8);
		ad.put("BetaHor", 0.8);
		ad.put("EmNorVer", 1.5e-6);
		ad.put("EmNorHor", 1.5e-6);
		ad.put("RevFrq", 11245.0);
		ad.put("NumBunInBeam", 338.0);
		ad.put("BunLen", 7.55e-2);
		ad.put("BunSpace", 25e-9);
		ad.put("Circum", 27.0);
		pd.put("LHC-p", ad);
		
		particleDataMap.put("lead", pd);

		
		pd = new ParticleData();
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2e12);
		ad.put("EnBeam", 750.0);
		ad.put("BetaVer", 0.01);
		ad.put("BetaHor", 0.01);
		ad.put("EmNorVer", 2.5e-5);
		ad.put("EmNorHor", 2.5e-5);
		ad.put("RepFrq", 15.0);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 0.01);
		ad.put("Circum", 2.5);
		ad.put("Circle", 1000.0);
		pd.put("MAF", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2e12);
		ad.put("EnBeam", 1500.0);
		ad.put("BetaVer", 0.005);
		ad.put("BetaHor", 0.005);
		ad.put("EmNorVer", 2.5e-5);
		ad.put("EmNorHor", 2.5e-5);
		ad.put("RepFrq", 12.0);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 0.005);
		ad.put("Circum", 4.5);
		ad.put("Circle", 1000.0);
		pd.put("MAF1500", ad);
		
		particleDataMap.put("muon", pd);
		
		pd = new ParticleData();
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2e12);
		ad.put("EnBeam", 750.0);
		ad.put("BetaVer", 0.01);
		ad.put("BetaHor", 0.01);
		ad.put("EmNorVer", 2.5e-5);
		ad.put("EmNorHor", 2.5e-5);
		ad.put("RepFrq", 15.0);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 0.01);
		ad.put("Circum", 2.5);
		ad.put("Circle", 1000.0);
		pd.put("MAF", ad);
		
		ad = new AcceleratorData();
		ad.put("NumParInBun", 2e12);
		ad.put("EnBeam", 1500.0);
		ad.put("BetaVer", 0.005);
		ad.put("BetaHor", 0.005);
		ad.put("EmNorVer", 2.5e-5);
		ad.put("EmNorHor", 2.5e-5);
		ad.put("RepFrq", 12.0);
		ad.put("NumBunInBeam", 1.0);
		ad.put("BunLen", 0.005);
		ad.put("Circum", 4.5);
		ad.put("Circle", 1000.0);
		pd.put("MAF1500", ad);
		
		
		particleDataMap.put("muon+", pd);

	}
	public void loadVariableData()
	{
		variableDataMap.put("NumParInBun", new VariableData("Number of particle per bunch(N):",""));
		variableDataMap.put("EnBeam", new VariableData("Particle beam energy:","GeV"));
		variableDataMap.put("BetaVer", new VariableData("Vertical Beta function of particle beam at IP:","m"));
		variableDataMap.put("BetaHor", new VariableData("Horizontal Beta function of particle beam at IP:","m"));
		variableDataMap.put("Beta", new VariableData("Beta function of particle beam at IP:","m"));
		variableDataMap.put("EmNorVer", new VariableData("Norm. Vertical Emittance of particle beam:","m"));
		variableDataMap.put("EmNorHor", new VariableData("Norm. Horizontal Emittance of particle beam:","m"));
		variableDataMap.put("EmNor", new VariableData("Norm. Emittance of particle beam:","m"));
		variableDataMap.put("PulFrq", new VariableData("Pulse Frequency of beams:","Hz"));
		variableDataMap.put("RevFrq", new VariableData("Revolution Frequency of beam:","Hz"));
		variableDataMap.put("RepFrq", new VariableData("Repetition Rate of beam:","Hz"));
		variableDataMap.put("ColFrq", new VariableData("Collision Frequency of beams:","Hz"));


		variableDataMap.put("NumBunInBeam", new VariableData("Bunches in particle beam:",""));
		variableDataMap.put("BunLen", new VariableData("Particle Beam Bunch Length:","m"));
		variableDataMap.put("BunSpace", new VariableData("Bunch Spacing of Particle Beam:","m"));
		variableDataMap.put("Circum", new VariableData("Circumference:","km"));
		variableDataMap.put("Circle", new VariableData("Number of Turns:","turn"));

		
		variableDataMap.put("PowLim", new VariableData("Power Limit for Particle Beam:","MW"));
		variableDataMap.put("DisrLim", new VariableData("Disruption Limit:",""));
		variableDataMap.put("BeamParLim", new VariableData("Beam-Beam Parameter Limit:",""));
		variableDataMap.put("NumRepPar", new VariableData("Number of Representative Particles:",""));
		variableDataMap.put("Angle", new VariableData("Angle:","rad"));
		variableDataMap.put("NumSteps", new VariableData("Number of Steps:",""));
		
		variableDataMap.put("noLum", new VariableData("Nominal Luminosity: ",""));
		variableDataMap.put("efLum", new VariableData("Effective Luminosity: ",""));
		variableDataMap.put("reFac", new VariableData("Enhancement/Reduction Factor: ",""));
		
		
		variableDataMap.put("disX", new VariableData("<html>  Disruption (D<sub>x</sub>): <html>",""));
		variableDataMap.put("disY", new VariableData("<html>  Disruption (D<sub>y</sub>): <html>",""));
		variableDataMap.put("divX", new VariableData("<html>  Divergence (A<sub>x</sub>): <html>",""));
		variableDataMap.put("divY", new VariableData("<html>  Divergence (A<sub>y</sub>): <html>",""));
		variableDataMap.put("bbTunX", new VariableData("<html>  BB Tuneshift (&xi<sub>x</sub>): <html>",""));
		variableDataMap.put("bbTunY", new VariableData("<html>  BB Tuneshift (&xi<sub>y</sub>): <html>",""));
		
		
	}
	
	public void loadSettingsData()
	{
		settingsData.put("PowLim", 50.0);
		settingsData.put("DisrLim", 25.0);
		settingsData.put("BeamParLim", 0.01);
		settingsData.put("NumRepPar", 1000.0);
		settingsData.put("NumSteps", 16.0);
		settingsData.put("Angle", 0.0);


	}
	public String getTitle(String key)
	{
		return variableDataMap.get(key).getTitle();
	}
	public String getUnit(String key)
	{
		return variableDataMap.get(key).getUnit();
	}
	public LinkedHashMap<String, VariableData> getVariableDataMap() {
		return variableDataMap;
	}
	public void setVariableDataMap(LinkedHashMap<String, VariableData> variableDataMap) {
		this.variableDataMap = variableDataMap;
	}
	public LinkedHashMap<String, Double> getSettingsData() {
		return settingsData;
	}
	public void setSettingsData(LinkedHashMap<String, Double> settingsData) {
		this.settingsData = settingsData;
	}
	public HashMap<String, ParticleData> getParticleDataMap() {
		return particleDataMap;
	}
	public void setParticleDataMap(HashMap<String, ParticleData> particleDataMap) {
		this.particleDataMap = particleDataMap;
	}
	
}
