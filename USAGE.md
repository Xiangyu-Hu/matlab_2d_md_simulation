# Usage Guide for 2D Molecular Dynamics Simulation

## Quick Start

### Running the Main Simulation

1. Open MATLAB
2. Navigate to the repository directory:
   ```matlab
   cd /path/to/matlab_2d_md_simulation
   ```
3. Run the main script:
   ```matlab
   md_simulation_2d
   ```

This will start a simulation with default parameters and display real-time visualization.

## Understanding the Output

### Visualization Panels

The simulation window shows four panels:

1. **Top-Left: Particle Positions**
   - Blue circles represent particles
   - Shows spatial distribution in the 2D box
   - Demonstrates clustering, movement, and phase behavior

2. **Top-Right: Velocity Vectors**
   - Arrows show particle velocities
   - Arrow length indicates speed
   - Arrow direction shows motion direction

3. **Bottom-Left: Energy Evolution**
   - Red line: Kinetic energy (motion)
   - Blue line: Potential energy (interactions)
   - Black line: Total energy (should be nearly constant)

4. **Bottom-Right: Temperature Evolution**
   - Green line: Current temperature
   - Red dashed line: Target temperature
   - Should equilibrate near target after initial transient

### Console Output

Example output:
```
Starting 2D Molecular Dynamics Simulation...
Number of particles: 50
Box size: 10.00 x 10.00
Time step: 0.010
Total steps: 1000

Step 100/1000: T = 0.982, E_total = -23.456
Step 200/1000: T = 1.003, E_total = -23.452
...
```

## Customizing Parameters

### Basic Parameters

Edit `md_simulation_2d.m` to change:

```matlab
N = 50;          % Number of particles (try 25, 50, 100)
L = 10;          % Box size (try 8, 10, 15)
dt = 0.01;       % Time step (smaller = more accurate, slower)
n_steps = 1000;  % Simulation length (try 500, 1000, 2000)
T_target = 1.0;  % Temperature (try 0.5, 1.0, 2.0)
```

### Effect of Parameters

**Number of Particles (N)**
- Small (16-25): Fast, good for learning
- Medium (50-100): Realistic behavior
- Large (>100): Slower, more statistical

**Temperature (T_target)**
- Low (0.3-0.7): Solid-like, ordered structure
- Medium (0.8-1.5): Liquid-like behavior
- High (2.0+): Gas-like, high mobility

**Box Size (L)**
- Affects density: ρ = N/L²
- Larger box = lower density
- Denser systems show more interactions

**Time Step (dt)**
- Too large: Unstable, energy drift
- Too small: Slow simulation
- Good range: 0.005-0.01

## Example Scenarios

### Scenario 1: Observe Solid-Like Behavior

```matlab
N = 36;
L = 8;
T_target = 0.3;  % Cold system
dt = 0.005;      % Stable integration
n_steps = 1000;
```

Expected: Particles form ordered lattice structure

### Scenario 2: Liquid Simulation

```matlab
N = 64;
L = 10;
T_target = 1.0;  % Moderate temperature
dt = 0.01;
n_steps = 2000;
```

Expected: Fluid motion, particles move but remain clustered

### Scenario 3: Gas Phase

```matlab
N = 36;
L = 12;          % Lower density
T_target = 2.5;  % High temperature
dt = 0.005;
n_steps = 1000;
```

Expected: High mobility, particles spread throughout box

## Additional Scripts

### Visualize Lennard-Jones Potential

```matlab
plot_lj_potential
```

Shows the interaction potential and force curves. Useful for understanding why particles behave the way they do.

### Run Example Scenarios

```matlab
example_runs
```

Automatically runs three different scenarios demonstrating various physical behaviors.

## Troubleshooting

### Simulation is Unstable / Energy Increases

**Problem**: Total energy keeps increasing, particles move too fast

**Solution**: Reduce time step
```matlab
dt = 0.005;  % Instead of 0.01
```

### Particles are Too Slow / Boring

**Problem**: Nothing seems to happen

**Solution**: Increase temperature or decrease density
```matlab
T_target = 2.0;  % Higher temperature
L = 15;          % Larger box (lower density)
```

### Simulation is Too Slow

**Problem**: Takes too long to run

**Solution**: Reduce number of particles or steps
```matlab
N = 25;         % Fewer particles
n_steps = 500;  % Shorter simulation
```

### All Particles Cluster in Corner

**Problem**: Particles collapse to high-density region

**Solution**: This may indicate:
- Temperature too low
- Time step too large
- Try rerunning (random initialization can affect this)

## Understanding the Physics

### Energy Conservation

Good simulation: Total energy varies by < 1%
```
Energy drift: < 1e-2 (good)
Energy drift: > 1e-1 (check time step)
```

### Temperature Equilibration

- Initial: Fluctuates (system equilibrating)
- After ~100 steps: Should stabilize near target
- Fluctuations: Normal (microcanonical ensemble)

### Time and Distance Units

The simulation uses reduced units:
- Length scale: σ (sigma)
- Energy scale: ε (epsilon)
- Time scale: σ√(m/ε)

These are dimensionless. To get real units, multiply by appropriate physical constants.

## Performance Tips

1. **Avoid very large N**: Performance scales as O(N²)
2. **Use appropriate cutoff**: Default 2.5σ is good
3. **Balance dt and n_steps**: Smaller dt needs more steps
4. **Update visualization less often**: Change `mod(step, 10)` to `mod(step, 20)`

## Educational Exercises

1. **Vary temperature**: Observe solid → liquid → gas transitions
2. **Energy conservation**: Check drift with different dt values
3. **Density effects**: Change N and L to vary ρ = N/L²
4. **Equilibration time**: How long until temperature stabilizes?
5. **Radial distribution**: Plot average particle spacing (advanced)

## Further Modifications

Ideas for extending the code:
- Add thermostat (velocity rescaling, Nosé-Hoover)
- Calculate radial distribution function g(r)
- Implement neighbor lists for speed
- Add different potentials (Morse, soft-sphere)
- Save trajectories for later analysis
- Compute pressure from virial
