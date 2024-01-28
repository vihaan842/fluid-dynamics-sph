import java.awt.*;
import java.io.EOFException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import javax.swing.JFrame;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Set;
public class Main {
	static int width = 1024;
	public static void main(String[] args) throws Exception {
		Simulation sim = new Simulation(1.0 / 30.0);
		Display displ = new Display(width, width, sim);
		JFrame frm = new JFrame();
		frm.add(displ);
		frm.setSize(width + 50, width + 50);
		//f.setLayout(null);
		frm.setVisible(true);
		Graphics gr = displ.getGraphics();
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		while (true) {
			displ.paint(gr);
			//br.readLine();
			sim.step();
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
		int pressureSamplingRate = 10;
		
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
		g.drawString("Fluid Simulator", 40, 990);
		g.drawString("t=" + Double.toString(sim.time) + "s", 40, 1010);
		
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
