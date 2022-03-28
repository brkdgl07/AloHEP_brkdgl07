package org.hepforge.alohep.gfx;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.HashMap;

import org.hepforge.alohep.AloHEP;

public class AnimationPanel{

	
	private String title = "AloHEP";
	private int fps = 60;
	private boolean running = false;
	private Thread render_thread;
	private Thread update_thread;
	
	private Display display;
	private Camera camera;
	private Input input;
	private HashMap<String, State> states;
	private State state;
	private AloHEP alohep;
	
	public AnimationPanel(AloHEP alohep)
	{
		this.alohep = alohep;
		camera = new Camera(alohep);
		input = new Input();
		states = new HashMap<String,State>();
		display = new Display(this,title);

	}
	public void start() {
		state.init(this);
		state.setInit(true);
		if(running)
			return;
		
		
		update_thread = new Thread()
		{
			@Override
			public void run()
			{
				running = true;
				while(running)
				{
					
					state.update();
					try {
						Thread.sleep(17);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		};
		render_thread = new Thread()
		{
			public void stopthread()
			{
				if(!running)
					return;
				running = false;
			}
			@Override
			public void run() 
			{
				running = true;
				double timePerTick = 1000000000 / fps;
				double delta = 0;
				long now;
				long lastTime = System.nanoTime();
				long timer = 0;
				int ticks = 0;
				while(running)
				{
					now = System.nanoTime();
					delta += (now - lastTime)/timePerTick;
					timer += now - lastTime;
					lastTime = now;
					if(delta >= 1)
					{
						display.getPanel().repaint();

						ticks++;
						delta--;
					}
					if(timer >= 1000000000)
					{
						//System.out.println("Frame: "+ticks+" ticks");
						ticks = 0;
						timer = 0;
					}
				}
				stopthread();
				
			}
		};
		
		update_thread.start();
		render_thread.start();		
	}
	public void render(Graphics g){

		if(state != null && state.isInit())
		{
			state.render(display.getShadow());
		}
		
	}

	public void setState(String s)
	{
		setState(states.get(s));
		
	}
	public void putState(String name,State state)
	{
		states.put(name, state);
	}
	public void remState(String name)
	{
		states.remove(name);
	}
	public State getState()
	{
		return state;
	}
	public String getTitle() {
		return title;
	}
	public int getFps() {
		return fps;
	}
	public boolean isRunning() {
		return running;
	}

	public Display getDisplay() {
		return display;
	}
	public Input getInput() {
		return input;
	}
	public void setTitle(String title) {
		this.title = title;
	}
	public void setFps(int fps) {
		this.fps = fps;
	}
	public void setRunning(boolean running) {
		this.running = running;
	}
	public Camera getCamera() {
		return camera;
	}
	public HashMap<String, State> getStates() {
		return states;
	}
	public void setStates(HashMap<String, State> states) {
		this.states = states;
	}
	public void setState(State state) {
		this.state = state;
	}

	
}
