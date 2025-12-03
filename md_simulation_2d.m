%% 2D Molecular Dynamics Simulation
% Educational code for simulating particle interactions in 2D
% Based on public domain molecular dynamics algorithms
%
% This script simulates particles interacting via Lennard-Jones potential
% in a 2D periodic box using velocity Verlet integration
%
% Author: Unknown (Public Resource)
% Modified for educational purposes

clear; close all; clc;

%% System Parameters
N = 50;                    % Number of particles
L = 10;                    % Box size (L x L)
dt = 0.01;                 % Time step
n_steps = 1000;            % Number of simulation steps
T_target = 1.0;            % Target temperature (reduced units)

% Lennard-Jones parameters (reduced units)
epsilon = 1.0;             % Energy scale
sigma = 1.0;               % Length scale
r_cut = 2.5 * sigma;       % Cutoff radius

%% Initialize positions
% Place particles on a square lattice with small random perturbations
n_side = ceil(sqrt(N));
spacing = L / n_side;
particle_count = 0;
positions = zeros(N, 2);

for i = 1:n_side
    for j = 1:n_side
        if particle_count < N
            particle_count = particle_count + 1;
            positions(particle_count, 1) = (i - 0.5) * spacing + 0.1 * (rand - 0.5);
            positions(particle_count, 2) = (j - 0.5) * spacing + 0.1 * (rand - 0.5);
        end
    end
end

%% Initialize velocities
% Random velocities with zero total momentum
velocities = randn(N, 2);
% Remove center of mass motion
velocities = velocities - mean(velocities, 1);
% Scale to target temperature
KE = 0.5 * sum(sum(velocities.^2));
T_current = KE / N;
velocities = velocities * sqrt(T_target / T_current);

%% Storage for analysis
energy_kinetic = zeros(n_steps, 1);
energy_potential = zeros(n_steps, 1);
energy_total = zeros(n_steps, 1);
temperature = zeros(n_steps, 1);

%% Initial force calculation
forces = calculate_forces(positions, L, epsilon, sigma, r_cut);

%% Visualization setup
figure('Position', [100, 100, 800, 700]);

%% Main simulation loop
fprintf('Starting 2D Molecular Dynamics Simulation...\n');
fprintf('Number of particles: %d\n', N);
fprintf('Box size: %.2f x %.2f\n', L, L);
fprintf('Time step: %.3f\n', dt);
fprintf('Total steps: %d\n\n', n_steps);

for step = 1:n_steps
    % Velocity Verlet integration - Step 1
    % Update positions: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
    positions = positions + velocities * dt + 0.5 * forces * dt^2;
    
    % Apply periodic boundary conditions
    positions = mod(positions, L);
    
    % Calculate forces at new positions
    forces_new = calculate_forces(positions, L, epsilon, sigma, r_cut);
    
    % Velocity Verlet integration - Step 2
    % Update velocities: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    velocities = velocities + 0.5 * (forces + forces_new) * dt;
    
    % Update forces
    forces = forces_new;
    
    % Calculate energy and temperature
    KE = 0.5 * sum(sum(velocities.^2));
    [PE, ~] = calculate_potential_energy(positions, L, epsilon, sigma, r_cut);
    
    energy_kinetic(step) = KE;
    energy_potential(step) = PE;
    energy_total(step) = KE + PE;
    temperature(step) = KE / N;
    
    % Visualization every 10 steps
    if mod(step, 10) == 0
        subplot(2, 2, 1);
        plot(positions(:, 1), positions(:, 2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        xlim([0, L]); ylim([0, L]);
        xlabel('X position'); ylabel('Y position');
        title(sprintf('Particle Positions (Step %d)', step));
        axis square; grid on;
        
        subplot(2, 2, 2);
        quiver(positions(:, 1), positions(:, 2), velocities(:, 1), velocities(:, 2), 0.5);
        xlim([0, L]); ylim([0, L]);
        xlabel('X position'); ylabel('Y position');
        title('Velocity Vectors');
        axis square; grid on;
        
        subplot(2, 2, 3);
        plot(1:step, energy_kinetic(1:step), 'r-', ...
             1:step, energy_potential(1:step), 'b-', ...
             1:step, energy_total(1:step), 'k-', 'LineWidth', 1.5);
        xlabel('Time step'); ylabel('Energy');
        legend('Kinetic', 'Potential', 'Total', 'Location', 'best');
        title('Energy Evolution');
        grid on;
        
        subplot(2, 2, 4);
        plot(1:step, temperature(1:step), 'g-', 'LineWidth', 1.5);
        hold on;
        plot([1, step], [T_target, T_target], 'r--', 'LineWidth', 1);
        hold off;
        xlabel('Time step'); ylabel('Temperature');
        legend('Current', 'Target', 'Location', 'best');
        title('Temperature Evolution');
        grid on;
        
        drawnow;
        
        % Print progress
        if mod(step, 100) == 0
            fprintf('Step %d/%d: T = %.3f, E_total = %.3f\n', ...
                    step, n_steps, temperature(step), energy_total(step));
        end
    end
end

fprintf('\nSimulation completed!\n');
fprintf('Average temperature: %.3f\n', mean(temperature(end-100:end)));
fprintf('Average total energy: %.3f\n', mean(energy_total(end-100:end)));
fprintf('Energy drift: %.3e\n', std(energy_total(end-100:end)));

%% Functions

function forces = calculate_forces(positions, L, epsilon, sigma, r_cut)
    % Calculate forces between particles using Lennard-Jones potential
    % with periodic boundary conditions
    
    N = size(positions, 1);
    forces = zeros(N, 2);
    
    for i = 1:N
        for j = i+1:N
            % Vector from i to j
            rij = positions(j, :) - positions(i, :);
            
            % Apply minimum image convention (periodic boundaries)
            rij = rij - L * round(rij / L);
            
            % Distance
            r = sqrt(sum(rij.^2));
            
            % Only calculate if within cutoff
            if r < r_cut && r > 0
                % Lennard-Jones force magnitude
                % F = 24*epsilon/r * (2*(sigma/r)^12 - (sigma/r)^6)
                sr6 = (sigma / r)^6;
                force_magnitude = 24 * epsilon / r * (2 * sr6^2 - sr6);
                
                % Force vector (direction from i to j)
                force_vector = force_magnitude * rij / r;
                
                % Newton's third law
                forces(i, :) = forces(i, :) + force_vector;
                forces(j, :) = forces(j, :) - force_vector;
            end
        end
    end
end

function [PE, virial] = calculate_potential_energy(positions, L, epsilon, sigma, r_cut)
    % Calculate total potential energy and virial
    
    N = size(positions, 1);
    PE = 0;
    virial = 0;
    
    for i = 1:N
        for j = i+1:N
            % Vector from i to j
            rij = positions(j, :) - positions(i, :);
            
            % Apply minimum image convention
            rij = rij - L * round(rij / L);
            
            % Distance
            r = sqrt(sum(rij.^2));
            
            % Only calculate if within cutoff (r > 0 prevents division by zero)
            if r < r_cut && r > 0
                % Lennard-Jones potential
                % U = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
                sr6 = (sigma / r)^6;
                PE = PE + 4 * epsilon * (sr6^2 - sr6);
                
                % Virial contribution
                virial = virial + 24 * epsilon * (2 * sr6^2 - sr6);
            end
        end
    end
end
