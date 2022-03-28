package org.hepforge.alohep.gfx;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.Transparency;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class Display {

	private AnimationPanel board;
	private JFrame frame;
	private JPanel panel;
	private boolean fullscreen = true;
	private String title;
	private BufferedImage image;
	private BufferedImage shadow;
	private GraphicsConfiguration config =
            GraphicsEnvironment.getLocalGraphicsEnvironment()
                .getDefaultScreenDevice()
                .getDefaultConfiguration();
	private Graphics graphics;
	private static int min_width = 800,min_height = 640;
	public static int width, height;
	
	public Display(AnimationPanel board,String title){
		this.board = board;
		this.title = title;
		createDisplay();

	}

	private void createDisplay(){
		frame = new JFrame(title);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		if(fullscreen)
		{
			frame.setExtendedState(JFrame.MAXIMIZED_BOTH);
			frame.setUndecorated(true);
		}
		else
		{
			frame.setSize(min_width,min_height);
		}
		frame.setResizable(false);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		width = frame.getSize().width;
		height = frame.getSize().height;
		image = create(width, height, false);
		shadow = create(width, height, true);
		graphics = image.createGraphics();
		panel = new JPanel()
				{

					private static final long serialVersionUID = 4180554716147604945L;

					@Override
					public void paintComponent(Graphics g)
					{
						super.paintComponent(g);
						board.render(g);
						graphics.drawImage(shadow,0,0, width, height, null);
						g.drawImage(image,0,0,width, height,null);
						graphics.setColor(Color.WHITE);
						graphics.fillRect(0, 0, width, height);
					}
				};
		panel.setPreferredSize(frame.getSize());
		panel.setMaximumSize(frame.getSize());
		panel.setMinimumSize(frame.getSize());
		panel.setFocusable(false);
		frame.addKeyListener(board.getInput());
		frame.add(panel);
		frame.pack();
		panel.addMouseListener(board.getInput());
		panel.addMouseWheelListener(board.getInput());
		
		
	}
    private BufferedImage create(int width, int height, boolean alpha) {
        return config.createCompatibleImage(width, height, alpha
                ? Transparency.TRANSLUCENT : Transparency.OPAQUE);
    }
	public boolean fullscreen()
	{
		if(!fullscreen)
		{
			frame.setExtendedState(JFrame.MAXIMIZED_BOTH);
			frame.setUndecorated(true);
			panel.setPreferredSize(frame.getSize());
			panel.setMaximumSize(frame.getSize());
			panel.setMinimumSize(frame.getSize());
			width = frame.getSize().width;
			height = frame.getSize().height;
			return fullscreen;
		}
		else
		{

			frame.setSize(min_width,min_height);
			frame.setUndecorated(false);
			panel.setPreferredSize(frame.getSize());
			panel.setMaximumSize(frame.getSize());
			panel.setMinimumSize(frame.getSize());
			width = frame.getSize().width;
			height = frame.getSize().height;
			return fullscreen;
		}
	}
	public JPanel getPanel(){
		return panel;
	}
	
	public JFrame getFrame(){
		return frame;
	}

	public BufferedImage getImage() {
		return image;
	}

	public void setImage(BufferedImage image) {
		this.image = image;
	}

	public BufferedImage getShadow() {
		return shadow;
	}

	public void setShadow(BufferedImage shadow) {
		this.shadow = shadow;
	}

	
}
