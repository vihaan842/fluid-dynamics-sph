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
		Simulation sim = new Simulation(0.1);
		Display displ = new Display(width, width, sim);
		JFrame frm = new JFrame();
		frm.add(displ);
		frm.setSize(width, width);
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
	protected int boundLeft;
	protected int boundTop;
	double scaling;
	Simulation sim;
	int radius;
	Display(int x, int y, Simulation s) {
		viewWidth = x;
		viewHeight = y;
		scaling = ((double) x) / Simulation.SIZE;
		sim = s;
		radius = (int) (Simulation.PARTICLE_RADIUS * scaling);
	}
	public void paint(Graphics g) {
		//TODO freeze display when re-painting
		g.clearRect(boundLeft, boundTop, viewWidth, viewHeight);
		g.drawString("Fluid Simulator", boundLeft + 40, boundTop + 40);
		g.drawString("t=" + Double.toString(sim.time) + "s", boundLeft + 40, boundTop + 80);
		for (double[] pos : sim.positions) {
			g.drawOval((int) (pos[0] * scaling), (int) (pos[1] * scaling), radius, radius);
		}
	}
}
