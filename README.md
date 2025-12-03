# 2D Molecular Dynamics Simulation in MATLAB

Educational code for simulating particle interactions in two dimensions using molecular dynamics.

## Overview

This repository contains MATLAB code for a 2D Molecular Dynamics (MD) simulation. The simulation demonstrates fundamental concepts in computational physics and chemistry by modeling the motion and interactions of particles in a two-dimensional periodic box.

**Source:** Based on public domain molecular dynamics algorithms (original author unknown).

## What is Molecular Dynamics?

Molecular Dynamics is a computational method that simulates the physical movements of atoms and molecules. The particles are allowed to interact for a fixed period of time, giving a view of the dynamic evolution of the system. This educational code simplifies the concept to 2D for easier visualization and understanding.

## Features

- **Lennard-Jones Potential**: Particles interact via the classic Lennard-Jones 12-6 potential, which models attractive and repulsive forces
- **Velocity Verlet Integration**: Uses a time-reversible and symplectic integration scheme for accurate energy conservation
- **Periodic Boundary Conditions**: Particles wrap around the simulation box, effectively simulating an infinite system
- **Real-time Visualization**: Shows particle positions, velocity vectors, energy evolution, and temperature
- **Energy Conservation**: Monitors kinetic, potential, and total energy throughout the simulation

## Requirements

- MATLAB (tested with R2019b and later)
- No additional toolboxes required

## How to Run

1. Clone or download this repository
2. Open MATLAB and navigate to the repository directory
3. Run the main simulation script:
   ```matlab
   md_simulation_2d
   ```

The simulation will run and display four panels:
1. **Particle Positions**: Current positions of all particles in the box
2. **Velocity Vectors**: Direction and magnitude of particle velocities
3. **Energy Evolution**: Kinetic, potential, and total energy over time
4. **Temperature Evolution**: System temperature compared to target temperature

## Parameters

You can modify the following parameters in `md_simulation_2d.m`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 50 | Number of particles in the simulation |
| `L` | 10 | Size of the simulation box (L × L) |
| `dt` | 0.01 | Time step for integration |
| `n_steps` | 1000 | Total number of simulation steps |
| `T_target` | 1.0 | Target temperature (in reduced units) |
| `epsilon` | 1.0 | Lennard-Jones energy parameter |
| `sigma` | 1.0 | Lennard-Jones length parameter |
| `r_cut` | 2.5 | Cutoff radius for force calculations |

## Physics Concepts Demonstrated

1. **Lennard-Jones Potential**: 
   - U(r) = 4ε[(σ/r)^12 - (σ/r)^6]
   - Models repulsion at short range and attraction at medium range

2. **Velocity Verlet Algorithm**:
   - First-order time-reversible integrator
   - Better energy conservation than simple Euler method

3. **Periodic Boundary Conditions**:
   - Minimum image convention
   - Eliminates surface effects

4. **Statistical Mechanics**:
   - Temperature from kinetic energy: T = KE/N
   - Energy equipartition in 2D

## Expected Output

The simulation produces:
- Real-time animated visualization of particle motion
- Console output showing progress and statistics
- Final statistics including average temperature and energy conservation

Example console output:
```
Starting 2D Molecular Dynamics Simulation...
Number of particles: 50
Box size: 10.00 x 10.00
Time step: 0.010
Total steps: 1000

Step 100/1000: T = 0.982, E_total = -23.456
Step 200/1000: T = 1.003, E_total = -23.452
...
Simulation completed!
Average temperature: 1.001
Average total energy: -23.450
Energy drift: 2.34e-02
```

## Educational Value

This code is designed for:
- Understanding basic molecular dynamics concepts
- Learning numerical integration methods
- Visualizing many-body physics
- Introduction to statistical mechanics
- Practice with MATLAB scientific computing

## Limitations

As educational code, this implementation:
- Is optimized for clarity rather than performance
- Uses nested loops instead of vectorization
- Limited to relatively small systems (< 200 particles)
- Uses reduced units (dimensionless)
- Does not include advanced features like thermostats or barostats

## Further Reading

To learn more about molecular dynamics:
- Frenkel, D., & Smit, B. (2001). *Understanding Molecular Simulation*
- Allen, M. P., & Tildesley, D. J. (2017). *Computer Simulation of Liquids*
- Rapaport, D. C. (2004). *The Art of Molecular Dynamics Simulation*

## License

This code is provided for educational purposes. The original source is based on public domain algorithms (author unknown). Feel free to use and modify for learning and teaching.

## Contributing

This is educational code intended for learning. If you find bugs or have suggestions for improvements that maintain the educational clarity, please feel free to contribute.