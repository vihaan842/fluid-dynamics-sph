import java.awt.*;
import java.io.EOFException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import javax.swing.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Set;
public class Main {
	static int width = 600;
	static JButton[] buttons;
	static ButtonHandler[] handlers;
	static final int BUTTON_COUNT = 2;
	static final int SLIDER_COUNT = 2;
	static JSlider[] sliders;
	static SliderHandler[] shandelers;
	static Simulation simulation;
	static double getMin(int i) throws Exception {
		switch (i) {
			case (0):
				return -12;
			case (1):
				return 5;
		}
		throw new Exception();
	}
	static double getMax(int i) throws Exception {
		switch (i) {
			case (0):
				return 12;
			case (1):
				return 100;
		}
		throw new Exception();
	}
	static double getVal(int i) throws Exception {
		return simulation.getProperty(i);
	}
	static void setVal(int i, double d) throws Exception {
		simulation.setProperty(i, d);
		return;
	}
	static String getName(int i) throws Exception {
		switch (i) {
			case (0):
				return "G (m/s)";
			case (1):
				return "Time step (ms)";
		}
		throw new Exception();
	}
	public static void main(String[] args) throws Exception {
		Simulation sim = simulation = new Simulation(1.0 / 15.0);
		Display displ = new Display(width, width, sim);
		JFrame frm = new JFrame();
		JPanel jp = new JPanel();
		jp.setLayout(new GridLayout(2, BUTTON_COUNT + SLIDER_COUNT));
		buttons = new JButton[BUTTON_COUNT];
		handlers = new ButtonHandler[BUTTON_COUNT];
		for (int i = 0; i < buttons.length; i++) {
			buttons[i] = new JButton("Button " + i);
			handlers[i] = new ButtonHandler(i);
			buttons[i].addActionListener(handlers[i]);
			jp.add(buttons[i]);
			jp.add(new JLabel("Button " + i));
		}
		sliders = new JSlider[BUTTON_COUNT];
		shandelers = new SliderHandler[BUTTON_COUNT];
		for (int i = 0; i < buttons.length; i++) {
			sliders[i] = new JSlider((int) getMin(i), (int) getMax(i), (int) getVal(i));
			sliders[i].setMajorTickSpacing((((int) getMax(i)) - ((int) getMin(i))) / 4);
			//sliders[i].setMinorTickSpacing(1);
			sliders[i].setPaintTicks(true);
			sliders[i].setPaintLabels(true);
			shandelers[i] = new SliderHandler(i);
			sliders[i].addChangeListener(shandelers[i]);
			jp.add(sliders[i]);
			jp.add(new JLabel(getName(i)));
		}
		frm.add(displ);
		frm.add(jp, BorderLayout.NORTH);
		frm.setSize(width + 50, width + 150);
		//f.setLayout(null);
		frm.setVisible(true);
		Graphics gr = displ.getGraphics();
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		while (true) {
		    for (int i = 0; i < 5; i++) {
			sim.step();
		    }
		    displ.paint(gr);
		}
	}
}
class ButtonHandler implements ActionListener {
	int id;
	JSlider slider;
	ButtonHandler(int i) {
		id = i;
	}
	public void actionPerformed(ActionEvent event) {
	}
}
class SliderHandler implements ChangeListener {
	int id;
	SliderHandler(int i) {
		id = i;
	}
	public void stateChanged(ChangeEvent event) {
		try {
			Main.setVal(id, Main.sliders[id].getValue());
		}
		catch (Exception e) {
			System.err.println("Exception in setting value for property " + id + ": " + e);
			e.printStackTrace();
		}
	}
}
class Display extends Canvas {
	protected int viewWidth;
	protected int viewHeight;
	double scaling;
	Simulation sim;
	int radius;
	double pressureCoeff = 150.0;
	double mvmtCoeff = 150.0;
	Display(int x, int y, Simulation s) {
		viewWidth = x;
		viewHeight = y;
		scaling = ((double) x) / Simulation.SIZE;
		sim = s;
		radius = (int) (Simulation.PARTICLE_RADIUS * scaling);
	}
	public void paint(Graphics g) {
		//TODO freeze display when re-painting
		//g.clearRect(0, 0, viewWidth, viewHeight);
		double maxPressure = 0.0;
		double maxMvmt = 0.0;
		int pressureSamplingRate = 5;
		
		for (int j = 0; j < viewHeight; j += pressureSamplingRate) {
			for (int i = 0; i < viewWidth; i += pressureSamplingRate) {
				double maxColor;
				double[] prs = pressure(((double) i) / scaling, ((double) j) / scaling);
				double pressure = Math.pow(prs[0], 0.2);
				double dx = 127.9 + (prs[1] * mvmtCoeff);
				double dy = 127.9 + (prs[2] * mvmtCoeff);
				double pres = pressure;
				pressure *= pressureCoeff;
				if (pressure == 0.0) {
					dx = 0;
					dy = 0;
				}
				double mC = Math.max(Math.abs(prs[1]), Math.abs(prs[2])) * 2.0;
				if (pressure > 255) {
					pressure = 255;
				}
				if (dx >= 256) {
					dx = 255;
				}
				else if (dx < 0) {
					dx = 0;
				}
				if (dy > 255) {
					dy = 255;
				}
				else if (dy < 0) {
					dy = 0;
				}
				g.setColor(new Color((int) pressure, (int) pressure, (int) pressure));
				g.fillRect(i, j, pressureSamplingRate, pressureSamplingRate);
				if (pres > maxPressure) {
					maxPressure = pres;
				}
				if (mC > maxMvmt) {
					maxMvmt = mC;
				}
			}
		}
		pressureCoeff = 200.0 / maxPressure;
		mvmtCoeff = 200.0 / maxMvmt;
		g.setColor(new Color(255, 255, 255));
		g.drawString("Fluid Simulator", 40, 550);
		g.drawString("t=" + Double.toString(sim.time) + "s", 40, 570);
		
		/*
		g.clearRect(0, 0, viewWidth, viewHeight);
		for (double[] pos : sim.positions) {
			g.drawOval((int) (pos[0] * scaling), (int) (pos[1] * scaling), radius, radius);
		}
		*/
	}
	double[] pressure(double x, double y) {
		double s = 0.0;
		double dx = 0.00;
		double dy = 0.00;
		double[] r = new double[]{x, y};
		Set<Integer>[] cs = sim.findChunks(r);
		for (Set<Integer> chunk: cs) {
		    if (chunk != null) {
			for (Integer particle: chunk) {
			    s += sim.kernel(Simulation.subtract(r, sim.positions[particle]));
			    dx += s * sim.velocities[particle][0];
			    dy += s * sim.velocities[particle][1];
			}
		    }
		}
		if (dx != 0.0) {
			dx /= s;
		}
		if (dy != 0.0) {
			dy /= s;
		}
		return new double[]{Simulation.PRESSURE_CONSTANT * Math.pow(s, Simulation.PRESSURE_POWER), dx, dy};
	}
}
