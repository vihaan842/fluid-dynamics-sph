import java.awt.*;
import java.io.EOFException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import javax.swing.JFrame;
import java.io.BufferedReader;
import java.io.InputStreamReader;
public class Main {
	static int width = 1024;
	public static void main(String[] args) throws Exception {
		Simulation sim = new Simulation(1.0 / 15.0);
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
		int pressureSamplingRate = 10;
		
		for (int j = 0; j < viewHeight; j += pressureSamplingRate) {
			for (int i = 0; i < viewWidth; i += pressureSamplingRate) {
				double pres;
				int pressure = (int) ((pres = Math.pow(pressure(((double) i) / scaling, ((double) j) / scaling), 0.2)) * pressureCoeff);
				if (pressure > 255) {
					pressure = 255;
				}
				g.setColor(new Color((int) pressure, (int) pressure, (int) pressure));
				g.fillRect(i, j, pressureSamplingRate, pressureSamplingRate);
				if (pres > maxPressure) {
					maxPressure = pres;
				}
			}
		}
		pressureCoeff = 200.0 / maxPressure;
		g.setColor(new Color(255, 255, 255));
		g.drawString("Fluid Simulator", 40, 940);
		g.drawString("t=" + Double.toString(sim.time) + "s", 40, 980);
		

		/*
		for (double[] pos : sim.positions) {
			g.drawOval((int) (pos[0] * scaling), (int) (pos[1] * scaling), radius, radius);
		}
		*/
	}
	double pressure(double x, double y) {
		double dim1min = x - Simulation.diameter;
		double dim1max = x + Simulation.diameter;
		double dim2min = y - Simulation.diameter;
		double dim2max = y + Simulation.diameter;
		double s = 0.0;
		double[] r = new double[]{x, y};
		for (int j = 0; j < sim.PARTICLES; j++) {
		    if (sim.positions[j][0] <= dim1min) {
			continue;
		    }
		    if (sim.positions[j][0] >= dim1max) {
			continue;
		    }
		    if (sim.positions[j][1] <= dim2min) {
			continue;
		    }
		    if (sim.positions[j][1] >= dim2max) {
			continue;
		    }
		    s += sim.kernel(Simulation.subtract(r, sim.positions[j]));
		}
		return Simulation.PRESSURE_CONSTANT * Math.pow(s, Simulation.PRESSURE_POWER);
	}
}
