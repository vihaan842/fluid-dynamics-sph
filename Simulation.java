import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.stream.IntStream;

public class Simulation {
    private static final int DIMENSIONS = 2;
    public static final double SIZE = 150.0;// meter
    public static final double PARTICLE_RADIUS = 6.0;// meter
    public static final double diameter = PARTICLE_RADIUS * 2;// meter
    private static final double SMOOTHING_LENGTH_SCALE = PARTICLE_RADIUS / 2.0;
    public static final int PARTICLES = 1500;
    public static final double PRESSURE_POWER = 1.0 + (1.0 / 2.0);
    public static final double DENSITY_POWER = PRESSURE_POWER - 1;
    public static final double PRESSURE_CONSTANT = 10.0;
    private static final double DAMPING_FACTOR = 0.9;
    private static final double NORMALIZATION_CONSTANT = 5.0 / (14.0 * Math.PI * SMOOTHING_LENGTH_SCALE * SMOOTHING_LENGTH_SCALE);// Assuming 2 dimensions
    private static final double G = 9.81;// meter/(second*second)
    private static final double[] ZERO_VECTOR = {0.0, 0.0};
    public double[][] positions;// meter
    public double[][] velocities;// meter/second
    public double[][] accelerations;// meter/(second*second)
    public double[][] new_accelerations;// meter/(second*second)
    public double[] densities;
    public final double chunkSize = diameter / 2;
    private Set<Integer>[][] chunks;
    private int chunk_size;
    public double timeStep;// seconds
    public double time;
    private double mvct = 0.0;
    public Simulation(double step) {
	chunk_size = (int)Math.ceil(SIZE/chunkSize);
	chunks = new Set[chunk_size][chunk_size];
	for (int i = 0; i < chunk_size; i++) {
	    for (int j = 0; j < chunk_size; j++) {
		chunks[i][j] = new CopyOnWriteArraySet();
	    }
	}
	positions = new double[PARTICLES][DIMENSIONS];
	for (int i = 0; i < PARTICLES; i++) {
	    for (int j = 0; j < DIMENSIONS; j++) {
		positions[i][j] = (Math.random() * SIZE * 0.5 + 30.0) / ((double) (2 - j));
	    }
	    chunks[(int)(positions[i][0]/chunkSize)][(int)(positions[i][1]/chunkSize)]
		.add(i);
	}
	velocities = new double[PARTICLES][DIMENSIONS];
	accelerations = new double[PARTICLES][DIMENSIONS];
	new_accelerations = new double[PARTICLES][DIMENSIONS];
	densities = new double[PARTICLES];
	timeStep = step;
	time = 0;
    }
    
    // find the magnitude of a vector with 2 elements
    public static double magnitude(double[] vec) {
	return Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
    }

    // subtract two vectors with 2 elements
    public static double[] subtract(double[] a, double[] b) {
	return new double[] {a[0]-b[0], a[1]-b[1]};
    }

    // add two vectors with 2 elements
    public static double[] add(double[] a, double[] b) {
	return new double[] {a[0]+b[0], a[1]+b[1]};
    }

    // scalar multiple of a vectors with 2 elements
    public static double[] scalar_multiple(double k, double[] vec) {
	return new double[] {k * vec[0], k * vec[1]};
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
	    return a * a * a * NORMALIZATION_CONSTANT;
	}
	// Otherwise, it is equal to (2-q)^3 - 4(1-q)^3
	else {
	    double a = 2.0-q;
	    double b = 1.0-q;
	    return (a * a * a - 4 * b * b * b) * NORMALIZATION_CONSTANT;
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
	    return ZERO_VECTOR;
	}
	else if (q >= 1) {
	    // We distributed the 1/q, and then use the x,y,etc. in position
	    return scalar_multiple((-4/q+6-3*q) * NORMALIZATION_CONSTANT, position);
	}
	else {
	    return scalar_multiple((-16/q+28-15*q) * NORMALIZATION_CONSTANT, position);
	}
    }

    public Set<Integer>[] findChunks(double[] r) {
	int chunkX = (int)(r[0] / chunkSize);
	int chunkY = (int)(r[1] / chunkSize);
	Set<Integer>[] cs = new Set[9];
	cs[0] = chunks[chunkX][chunkY];
	if (chunkX > 0) {
	    if (chunkY > 0) {
		cs[1] = chunks[chunkX-1][chunkY-1];
	    }
	    cs[2] = chunks[chunkX-1][chunkY];
	    if (chunkY < chunk_size-1) {
		cs[3] = chunks[chunkX-1][chunkY+1];
	    }
	}
	if (chunkX < (chunk_size - 1)) {
	    if (chunkY > 0) {
		cs[4] = chunks[chunkX+1][chunkY-1];
	    }
	    cs[5] = chunks[chunkX+1][chunkY];
	    if (chunkY < chunk_size-1) {
		cs[6] = chunks[chunkX+1][chunkY+1];
	    }
	}
	if (chunkY > 0) {
	    cs[7] = chunks[chunkX][chunkY-1];
	}
	if (chunkY < (chunk_size - 1)) {
	    cs[8] = chunks[chunkX][chunkY+1];
	}
	return cs;
    }
    
    public double density(int index) {
	// Calculating the exact density is too difficult, so we use an approximation
	// We estimate the density around a single particle as the sum of the masses
	// times the smoothing kernel of the distance of all the particles.
	// Since the mass is the same for all particles, we just sum the smoothing kernel values
	double[] r = positions[index];
	Set<Integer>[] cs = findChunks(r);

	double s = 0.0;
	for (Set<Integer> chunk: cs) {
	    if (chunk != null) {
		for (Integer particle: chunk) {
		    s += kernel(subtract(r, positions[particle]));
		}
	    }
	}
	return s;
    }

    public void step() {
	/*
	mvct += 1.0 / (Math.random() * 10.0 * timeStep);
	while (mvct >= 1.0) {
		mvct--;
		int i = (int) (Math.random() * PARTICLES);
		positions[i][0] = 20.0 + (Math.random() * 10.0);
		positions[i][1] = SIZE - 40 - (Math.random() * 100.0);
		velocities[i][0] = Math.random();
		velocities[i][1] = Math.random();
	}
	*/
	// we first update the positions using the old velocity and acceleration values,
	// making sure to update the chunk values if necessary
	IntStream.range(0, PARTICLES).parallel().forEach((i) -> {
		int oldChunkX = (int)(positions[i][0] / chunkSize);
		int oldChunkY = (int)(positions[i][1] / chunkSize);
		positions[i] = add(positions[i],
				   add(scalar_multiple(timeStep, velocities[i]),
				       scalar_multiple(timeStep*timeStep/2.0,
						       accelerations[i])));
		if (positions[i][0] < 0.0 || positions[i][0] >= SIZE) {// TODO Prevent particles from becoming stuck in walls when `timeStep' is small (less than 0.03) and when the damping factor is very high (greater than 0.95)
			positions[i][0] -= DAMPING_FACTOR * 1.5 * velocities[i][0] * timeStep;
			velocities[i][0] *= DAMPING_FACTOR - 1.0;
		}
		if (positions[i][1] < 0.0 || positions[i][1] >= SIZE) {
			positions[i][1] -= DAMPING_FACTOR * 1.5 * velocities[i][1] * timeStep;
			velocities[i][1] *= DAMPING_FACTOR - 1.0;
		}
		if (positions[i][0] < 0) {
			positions[i][0] = 0;
		}
		else if (positions[i][0] >= SIZE) {
			positions[i][0] = SIZE - 0.001;
		}
		if (positions[i][1] < 0) {
			positions[i][1] = 0;
		}
		else if (positions[i][1] >= SIZE) {
			positions[i][1] = SIZE - 0.001;
		}
		if ((int)(positions[i][0] / chunkSize) != oldChunkX ||
		    (int)(positions[i][1] / chunkSize) != oldChunkY) {
		    chunks[oldChunkX][oldChunkY].remove(i);
		    chunks[(int)(positions[i][0] / chunkSize)][(int)(positions[i][1] / chunkSize)]
			.add(i);
		}
	    });
	// we then calculate the densities around each particle
	IntStream.range(0, PARTICLES).parallel().forEach((i) -> {
		densities[i] = Math.pow(density(i), DENSITY_POWER);
	    });
	// we then use this to determine the new accelerations for all the particles
	IntStream.range(0, PARTICLES).parallel().forEach((i) -> {
		Set<Integer>[] cs = findChunks(positions[i]);
		new_accelerations[i][0] = 0.0;
		new_accelerations[i][1] = G;
		// <Insert explanation of derivation>
		double s = 0.0;

		for (Set<Integer> chunk: cs) {
		    if (chunk != null) {
			for (Integer j: chunk) {
			    if (i != j) {
				new_accelerations[i] =
				    add(new_accelerations[i],
					scalar_multiple(-PRESSURE_CONSTANT *
							(
							 densities[i] + 
							 densities[j]
							 ),
							kernel_gradient(subtract(positions[i],
										 positions[j]))));
			    }
			}
		    }
		}
		// new_accelerations[i] =
		    // add(new_accelerations[i],
			// scalar_multiple(DAMPING_FACTOR-1.0, velocities[i]));
	    });
	// now, we have to update the velocities
	IntStream.range(0, PARTICLES).parallel().forEach((i) -> {
		velocities[i] = add(velocities[i],
				    scalar_multiple(timeStep/2,
						    add(accelerations[i],
							new_accelerations[i])));
		// only for 2d
		// if (positions[i][0] < 0.0 || positions[i][0] >= SIZE) {
		//     velocities[i][0] *= -1.0;
		// }
		// if (positions[i][1] < 0.0 || positions[i][1] >= SIZE) {
		//     velocities[i][1] *= -1.0;
		// }
	    });
	accelerations = new_accelerations;
	time += timeStep;
    }
}
