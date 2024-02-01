/*
 * Simulation.java
 */
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.stream.IntStream;
public class Simulation {
    private static final int DIMENSIONS = 2;
    public static final double SIZE = 150.0;// meter
    public static final double PARTICLE_RADIUS = 6.0;// meter
    public static final double diameter = PARTICLE_RADIUS * 2;// meter
    private static final double SMOOTHING_LENGTH_SCALE = PARTICLE_RADIUS / 2.0;
    public static final double PRESSURE_POWER = 1.0 + (1.0 / 2.0);
    public static final double DENSITY_POWER = PRESSURE_POWER - 1;
    public static final double PRESSURE_CONSTANT = 10.0;
    private static final double DAMPING_FACTOR = 0.9;
    private static final double NORMALIZATION_CONSTANT = 5.0 / (14.0 * Math.PI * SMOOTHING_LENGTH_SCALE * SMOOTHING_LENGTH_SCALE);// Assuming 2 dimensions
    private static double VISCOSITY = 0.125; // note: this has no real world analog; this is just for simulation purposes
    private static double G = 9.81;// meter/(second*second)
    private static final double[] ZERO_VECTOR = {0.0, 0.0};
    public int particles = 3000;
    public int capacity = particles;
    public double[][] positions;// meter
    public double[][] velocities;// meter/second
    public double[][] viscosities;// meter/second
    public double[][] accelerations;// meter/(second*second)
    public double[][] new_accelerations;// meter/(second*second)
    public double[] densities;
    public final int chunkScale = 2;
    public final double chunkSize = diameter / (double)chunkScale;
    private Set<Integer>[][] chunks;
    private int numChunks;
    private double timeStep;// seconds
    public double time;
    private double mvct = 0.0;
    public Simulation(double step) {// This simply sets up the simulation
	numChunks = (int)Math.ceil(SIZE/chunkSize);
	chunks = new Set[numChunks][numChunks];
	for (int i = 0; i < numChunks; i++) {
	    for (int j = 0; j < numChunks; j++) {
		chunks[i][j] = new CopyOnWriteArraySet();
	    }
	}
	positions = new double[particles][DIMENSIONS];
	for (int i = 0; i < particles; i++) {// Distribute the particles
	    for (int j = 0; j < DIMENSIONS; j++) {
		positions[i][j] = (Math.random() * SIZE * 0.5 + 30.0) / ((double) (2 - j));
	    }
	    chunks[(int)(positions[i][0]/chunkSize)][(int)(positions[i][1]/chunkSize)]
		.add(i);
	}
	velocities = new double[particles][DIMENSIONS];
	viscosities = new double[particles][DIMENSIONS];
	accelerations = new double[particles][DIMENSIONS];
	new_accelerations = new double[particles][DIMENSIONS];
	densities = new double[particles];
	timeStep = step;
	time = 0;// Set the elapsed time to 0 seconds
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
	// Calculate the kernel at a point, considering the
	// particle as being at the origin
	// The kernel is used to model interactions between particles,
	// with particles that are closer getting more influence
	// We ignore particles that are outside of the radius
	
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
	// Find the gradient of the kernel at a point, considering
	//the particle as being at the origin
	
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
	// Finding the chunks that may effect a particle, i.e. the chunks that are actually in the
	// radius of the kernel
	// This allows us to skip a bunch of particles and greatly speeds up the simulation
	// we find a square of chunks that is the appropriate size and is centered on the chunk that
	// the particle is in
	int chunkX = (int)(r[0] / chunkSize);
	int chunkY = (int)(r[1] / chunkSize);
	Set<Integer>[] cs = new Set[(2*chunkScale+1)*(2*chunkScale+1)];
	for (int i = -chunkScale; i < chunkScale+1; i++) {
	    if (chunkX + i >= 0 && chunkX + i < numChunks) {
		for (int j = -chunkScale; j < chunkScale+1; j++) {
		    if (chunkY + j >= 0 && chunkY + j < numChunks) {
			cs[i+chunkScale+(j+chunkScale)*(2*chunkScale+1)] =
			    chunks[chunkX + i][chunkY + j];
		    }
		};
	    }
	};
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

    public synchronized void step() {
	// This procedure performs simulation for 1 step of time
	// This uses "leap frog" integration, which is second order instead of first
	// order like Euler's method, which is good because our equation is for acceleration.
	// Usually, leap frog looks like:
	//  a(i) = A(x(i))
	//  v(i+1/2) = v(i-1/2) + a(i)delta(t)
	//  x(i+1) = x(i) + v(i+1/2)delta(t)
	// However, this is not very useful because it requires us to calculate half steps
	// Instead, we turn it into a form that uses full steps:
	//  x(i+1) = x(i) + v(i)delta(t) + a(t)/2*delta(t)^2
	//  v(i+1) = v(i) + 1/2(a(i) + a(i+1))delta(t)
	
	// we first update the positions using the old velocity and acceleration values,
	// making sure to update the chunk values if necessary
	IntStream.range(0, particles).parallel().forEach((i) -> {
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
		int newChunkX = (int)(positions[i][0] / chunkSize);
		int newChunkY = (int)(positions[i][1] / chunkSize);
		if ((newChunkX != oldChunkX ||
		     newChunkY != oldChunkY) &&
		    newChunkX >= 0 && newChunkX < numChunks &&
		    newChunkY >= 0 && newChunkY < numChunks &&
		    oldChunkX >= 0 && oldChunkX < numChunks &&
		    oldChunkY >= 0 && oldChunkY < numChunks) {
		    chunks[oldChunkX][oldChunkY].remove(i);
		    chunks[newChunkX][newChunkY].add(i);
		}
	    });
	// we then calculate the densities around each particle
	IntStream.range(0, particles).parallel().forEach((i) -> {
		densities[i] = Math.max(0.001, Math.pow(density(i), DENSITY_POWER));
	    });
	// we then use this to determine the new accelerations for all the particles
	IntStream.range(0, particles).parallel().forEach((i) -> {// For each particle
		// Please refer to https://pmocz.github.io/manuscripts/pmocz_sph.pdf for a full explanation
		// The following is a fast summary:
		// Remember, a(i) = A(x(i))
		// We know that A(r) = int(A(r')d(r-r')dr'), where d is the Dirac delta function
		// We then approximate d as the smoothing kernel, and we get
		//  A(r) = int(A(r')W(r-r')dr'), where W is the smoothing kernel
		// We then turn this into a discrete sum
		//  A(r) = sum(A(x(j))W(r-r(j))delta(V(j)))
		// We can now find the gradient and lapacian quite easily
		// We then consider Euler's equation:
		//  ro*dv/dt = -gradient(P) + f
		// We can rewrite this as:
		//  gradient(P)/ro = gradient(P/ro) + P/ro^2 * gradient(ro)
		// Then, we can use this to write:
		//  dv(i)/dt = -sum((P(i)/ro(i)^2 + P(j)/ro(j)^2)gradient(W(r - r(j)))) + b(i)
		
		new_accelerations[i][0] = 0.0;// We keep a running sum of accelerations / forces
		new_accelerations[i][1] = G;// Add the gravitational force to the acceleration forces
		double s = 0.0;
		Set<Integer>[] cs = findChunks(positions[i]);// We find the chunks of space that may contain other particles that may
		// be effected by the particle throught repulsive forces
		for (Set<Integer> chunk: cs) {// For each of those chunks, 
		    if (chunk != null) {
			for (Integer j: chunk) {// For each particle j within the chunk,
			    if (i != j) {// except for particle i,
				new_accelerations[i] =
				    add(new_accelerations[i],// We add to the running acceleration sum
					scalar_multiple(-PRESSURE_CONSTANT *// The product of
							(
							 densities[i] + 
							 densities[j]
							 ),// the calculated pressure
							kernel_gradient(subtract(positions[i],
										 positions[j]))));// and the gradient of the kernel function of particle j at particle i's position
			    }
			}
		    }
		}
	    });
	// Now, we have to update the velocities
	// This is again according to the leap frog technique
	IntStream.range(0, particles).parallel().forEach(i -> {// For each particle
		velocities[i] = add(velocities[i],// v + (timeStep / 2) * (a + a_new) -> v
				    scalar_multiple(timeStep/2,
						    add(accelerations[i],
							new_accelerations[i])));
	    });
	// We now calculate the viscosities, which ensure that our particles behave
	// more like a fluid rather than a bunch of particles
	// We simply add this onto the velocities afterwards
	IntStream.range(0, particles).parallel().forEach(i -> {
	        Set<Integer>[] cs = findChunks(positions[i]);
		viscosities[i] = ZERO_VECTOR;
		for (Set<Integer> chunk: cs) {
		    if (chunk != null) {
			for (Integer j: chunk) {
			    if (i != j) {
				double[] r = subtract(positions[j], positions[i]);
				viscosities[i] = add(viscosities[i],
						     scalar_multiple(VISCOSITY *
								     kernel(r) /
								     densities[j],
								     subtract(velocities[j], velocities[i])));
			    }
			}
		    }
		}
	    });
	// we add the viscosities to the velocities
	IntStream.range(0, particles).parallel().forEach(i -> {
		velocities[i] = add(velocities[i], viscosities[i]);
	    });
	// We update the accelerations
	accelerations = new_accelerations;
	time += timeStep;
    }
    public void addParticles(int n, double x, double y, double radius) {
	// This function simply adds a specified amount of 
	// particles within a specified radius of a specified point
	// This has various problems right now and we recommend that you
	// not use it.
	if (particles + n > capacity) {
	    System.out.println("Increasing capacity");
	    double[][] old_positions = positions;
	    double[][] old_velocities = velocities;
	    double[][] old_accelerations = accelerations;
	    double[][] old_new_accelerations = new_accelerations;
	    capacity *= 2;
	    positions = new double[capacity][DIMENSIONS];
	    velocities = new double[capacity][DIMENSIONS];
	    viscosities = new double[capacity][DIMENSIONS];
	    accelerations = new double[capacity][DIMENSIONS];
	    new_accelerations = new double[capacity][DIMENSIONS];
	    densities = new double[capacity];
	    IntStream.range(0, particles).parallel().forEach(i -> {
		    positions[i] = old_positions[i];
		    velocities[i] = old_velocities[i];
		    accelerations[i] = old_accelerations[i];
		    new_accelerations[i] = old_new_accelerations[i];
		});
	}

	IntStream.range(0, n).parallel().forEach(i -> {
		double r = Math.random() * radius;
		double a = Math.random() * 2 * 3.14;
		positions[particles+i][0] = x + r * Math.cos(a);
		positions[particles+i][1] = y + r * Math.sin(a);
		chunks[(int)((x + r * Math.cos(a)) / chunkSize)][(int)((y + r * Math.sin(a)) / chunkSize)]
		    .add(particles+i);
	    });
	
	particles += n;
    }
    public synchronized double getProperty(int i) throws Exception {
	// For fetching simulation properties to show in the UI
	switch (i) {
	case (0):
	    return G;
	case (1):
	    return timeStep * 1000;
	case (2):
	    return VISCOSITY * 1000;
	}
	throw new Exception();
    }
    public synchronized void setProperty(int i, double d) throws Exception {
	// For changing simulation properties from the UI
	switch (i) {
	case (0):
	    G = d;
	    return;
	case (1):
	    timeStep = d / 1000.0;
	    return;
	case (2):
	    VISCOSITY = d / 1000.0;
	    return;
	}
	throw new Exception();
    }
}
