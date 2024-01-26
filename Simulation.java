public class Simulation {
    private static final int DIMENSIONS = 2;
    public static final double SIZE = 10.0;
    private static final double SMOOTHING_LENGTH_SCALE = 0.5;
    public static final int PARTICLES = 500;
    private static final double PRESSURE_POWER = 1.0 + (1.0 / 2.0);
    private static final double DENSITY_POWER = PRESSURE_POWER - 2;
    private static final double PRESSURE_CONSTANT = 0.01;
    private static final double WALL_CONSTANT = 0.5;
    private static final double DAMPING_FACTOR = 0.2;
    private static final double G = 0.001;
    private static final double[] ZERO = {0.0, 0.0};
    public double[][] positions;
    public double[][] velocities;
    public Simulation() {
	positions = new double[PARTICLES][DIMENSIONS];
	for (int i = 0; i < PARTICLES; i++) {
	    for (int j = 0; j < DIMENSIONS; j++) {
		positions[i][j] = Math.random() * SIZE;
	    }
	}
	velocities = new double[PARTICLES][DIMENSIONS];
    }
    
    // find the magnitude of a vector with 2 elements
    public double magnitude(double[] vec) {
	return Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
    }

    // subtract two vectors with 2 elements
    public double[] subtract(double[] a, double[] b) {
	return new double[] {a[0]-b[0], a[1]-b[1]};
    }

    // add two vectors with 2 elements
    public double[] add(double[] a, double[] b) {
	return new double[] {a[0]+b[0], a[1]+b[1]};
    }

    // subtract two vectors with 2 elements
    public double[] multiply(double c, double[] vec) {
	return new double[] {c*vec[0], c*vec[1]};
    }
    
    public double kernel(double[] position) {
	// We will use the cubic spline smoothing kernel
	// We first calculate q = ||r||/h
	double q = magnitude(position) / SMOOTHING_LENGTH_SCALE;
	// If q >= 2, then the spline is equal to 0
	if (q >= 2) {
	    return 0.0;
	}
	// If q is between 1 and 2, it is equal to (2-q)^3
	else if (q >= 1) {
	    double a = 2.0-q;
	    return a*a*a;
	}
	// Otherwise, it is equal to (2-q)^3 - 4(1-q)^3
	else {
	    double a = 2.0-q;
	    double b = 1.0-q;
	    return a*a*a-4*b*b*b;
	}
    }

    public double[] kernel_gradient(double[] position) {
	// To find this gradient, we note that q = sqrt(x^2+y^2+...)
	// We let W = value of the kernel
	// dW/dx = dW/dq * dq/dx
	// (2-q)^3 = 8 - 4q + 2q^2 - q^3
	// d/dq (8 - 4q + 2q^2 - q^3) = -4 + 6q - 3q^2
	// (2-q)^3+4(1-q)^3 = 8 - 4q + 2q^2 - q^3 + 4 - 12q + 12q^2 - 4q^3
	//                  = 12 - 16q + 14q^2 - 5q^3
	// d/dq (12 - 16q + 14q^2 - 5q^3) = -16 + 28q - 15q^2
	// dq/dx = 2x/2/sqrt(x^2+y^2+...) = x / q
	double q = magnitude(position) / SMOOTHING_LENGTH_SCALE;
	if (q >= 2) {
	    return ZERO;
	}
	else if (q >= 1) {
	    // We distributed the 1/q, and then use the x,y,etc. in position
	    return multiply(-4/q+6-3*q, position);
	}
	else {
	    return multiply(-16/q+28-15*q, position);
	}
    }
    
    public double density(int index) {
	// Calculating the exact density is too difficult, so we use an approximation
	// We estimate the density around a single particle as the sum of the masses
	// times the smoothing kernel of the distance of all the particles.
	// Since we assume that the mass is 0, we just sum the smoothing kernel values
	double[] r = positions[index];
	double s = 0.0;
	for (int i = 0; i < PARTICLES; i++) {
	    s += kernel(subtract(r, positions[i]));
	}
	return s;
    }

    public void step() {
	// we first calculate the densities around each particle
	double[] densities = new double[PARTICLES];
	for (int i = 0; i < PARTICLES; i++) {
	    densities[i] = density(i);
	}
	// we then use this to determine the acceleration for all the particles
	double[][] acceleration = new double[PARTICLES][DIMENSIONS];
	for (int i = 0; i < PARTICLES; i++) {
	    // <Insert explanation of derivation>
	    double s = 0.0;
	    for (int j = 0; j < PARTICLES; j++) {
		if (i != j) {
		    acceleration[i] =
			add(acceleration[i],
			    multiply(-PRESSURE_CONSTANT *
				     (
				      Math.pow(densities[i], DENSITY_POWER) + 
				      Math.pow(densities[j], DENSITY_POWER)
				      ),
				     kernel_gradient(subtract(positions[i],
							      positions[j]))));
		    acceleration[i][1] += G;
		}
	    }
	}
	// now, we have to update the velocities
	for (int i = 0; i < PARTICLES; i++) {
	    velocities[i] = add(velocities[i], acceleration[i]);
	    // only for 2d
	    if (positions[i][0] < 0.0 || positions[i][0] >= SIZE) {
		velocities[i][0] *= DAMPING_FACTOR-1.0;
	    }
	    if (positions[i][1] < 0.0 || positions[i][1] >= SIZE) {
		velocities[i][1] *= DAMPING_FACTOR-1.0;
	    }
	}
	// now, we update the positions
	for (int i = 0; i < PARTICLES; i++) {
	    positions[i] = add(positions[i], velocities[i]);
	}
    }
}
