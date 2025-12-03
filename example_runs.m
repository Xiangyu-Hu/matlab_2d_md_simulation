%% Example Runs for 2D Molecular Dynamics Simulation
% This script demonstrates different simulation scenarios
% 
% Each example showcases different physical behaviors by varying
% the system parameters

%% Example 1: Small System - Quick Demonstration
% Fast simulation with few particles for quick testing
fprintf('=== Example 1: Small System ===\n');
N = 16;                    % Fewer particles
L = 5;                     % Smaller box
dt = 0.01;
n_steps = 500;             % Fewer steps
T_target = 1.0;

run_simulation(N, L, dt, n_steps, T_target);
pause(2);

%% Example 2: Cold System - Low Temperature
% Observe crystallization-like behavior at low temperature
fprintf('\n=== Example 2: Cold System ===\n');
N = 36;
L = 8;
dt = 0.005;                % Smaller time step for stability
n_steps = 1000;
T_target = 0.3;            % Low temperature

run_simulation(N, L, dt, n_steps, T_target);
pause(2);

%% Example 3: Hot System - High Temperature
% Observe gas-like behavior at high temperature
fprintf('\n=== Example 3: Hot System ===\n');
N = 36;
L = 10;
dt = 0.005;
n_steps = 1000;
T_target = 3.0;            % High temperature

run_simulation(N, L, dt, n_steps, T_target);

%% Helper Function
function run_simulation(N, L, dt, n_steps, T_target)
    % Run a molecular dynamics simulation with specified parameters
    % This is a simplified version for demonstration
    
    epsilon = 1.0;
    sigma = 1.0;
    r_cut = 2.5 * sigma;
    
    % Initialize positions on lattice
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
    
    % Initialize velocities
    velocities = randn(N, 2);
    velocities = velocities - mean(velocities, 1);
    KE = 0.5 * sum(sum(velocities.^2));
    T_current = KE / N;
    velocities = velocities * sqrt(T_target / T_current);
    
    % Storage
    energy_total = zeros(n_steps, 1);
    temperature = zeros(n_steps, 1);
    
    % Initial forces
    forces = calculate_forces_simple(positions, L, epsilon, sigma, r_cut);
    
    % Setup figure
    figure('Position', [100, 100, 800, 400]);
    
    % Main loop
    for step = 1:n_steps
        % Velocity Verlet
        positions = positions + velocities * dt + 0.5 * forces * dt^2;
        positions = mod(positions, L);
        forces_new = calculate_forces_simple(positions, L, epsilon, sigma, r_cut);
        velocities = velocities + 0.5 * (forces + forces_new) * dt;
        forces = forces_new;
        
        % Calculate properties
        KE = 0.5 * sum(sum(velocities.^2));
        PE = calculate_potential_simple(positions, L, epsilon, sigma, r_cut);
        
        energy_total(step) = KE + PE;
        temperature(step) = KE / N;
        
        % Visualize
        if mod(step, 20) == 0
            subplot(1, 2, 1);
            plot(positions(:, 1), positions(:, 2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
            xlim([0, L]); ylim([0, L]);
            xlabel('X'); ylabel('Y');
            title(sprintf('Step %d, T=%.2f', step, temperature(step)));
            axis square; grid on;
            
            subplot(1, 2, 2);
            plot(1:step, temperature(1:step), 'r-', 'LineWidth', 2);
            hold on;
            plot([1, step], [T_target, T_target], 'k--', 'LineWidth', 1);
            hold off;
            xlabel('Step'); ylabel('Temperature');
            legend('Current', 'Target');
            title('Temperature');
            grid on;
            
            drawnow;
        end
    end
    
    fprintf('Completed: N=%d, T_target=%.2f, Final T=%.3f\n', ...
            N, T_target, mean(temperature(end-50:end)));
end

function forces = calculate_forces_simple(positions, L, epsilon, sigma, r_cut)
    N = size(positions, 1);
    forces = zeros(N, 2);
    
    for i = 1:N
        for j = i+1:N
            rij = positions(j, :) - positions(i, :);
            rij = rij - L * round(rij / L);
            r = sqrt(sum(rij.^2));
            
            if r < r_cut && r > 0
                sr6 = (sigma / r)^6;
                force_magnitude = 24 * epsilon / r * (2 * sr6^2 - sr6);
                force_vector = force_magnitude * rij / r;
                forces(i, :) = forces(i, :) + force_vector;
                forces(j, :) = forces(j, :) - force_vector;
            end
        end
    end
end

function PE = calculate_potential_simple(positions, L, epsilon, sigma, r_cut)
    N = size(positions, 1);
    PE = 0;
    
    for i = 1:N
        for j = i+1:N
            rij = positions(j, :) - positions(i, :);
            rij = rij - L * round(rij / L);
            r = sqrt(sum(rij.^2));
            
            if r < r_cut && r > 0
                sr6 = (sigma / r)^6;
                PE = PE + 4 * epsilon * (sr6^2 - sr6);
            end
        end
    end
end
