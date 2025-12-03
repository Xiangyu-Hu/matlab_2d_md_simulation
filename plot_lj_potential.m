%% Visualize Lennard-Jones Potential and Force
% Educational script to understand the particle interaction model
%
% This script plots the Lennard-Jones potential energy and force
% as functions of particle separation distance

close all; clear all; clc;

%% Parameters
sigma = 1.0;      % Length scale
epsilon = 1.0;    % Energy scale

% Distance range
r_min = 0.8 * sigma;
r_max = 3.0 * sigma;
r = linspace(r_min, r_max, 1000);

%% Calculate Lennard-Jones Potential
% U(r) = 4*epsilon * [(sigma/r)^12 - (sigma/r)^6]
U = zeros(size(r));
for i = 1:length(r)
    sr6 = (sigma / r(i))^6;
    U(i) = 4 * epsilon * (sr6^2 - sr6);
end

%% Calculate Force
% F(r) = 24*epsilon/r * [2*(sigma/r)^12 - (sigma/r)^6]
F = zeros(size(r));
for i = 1:length(r)
    sr6 = (sigma / r(i))^6;
    F(i) = 24 * epsilon / r(i) * (2 * sr6^2 - sr6);
end

%% Find equilibrium distance
% Minimum of potential occurs at r = 2^(1/6) * sigma
r_eq = 2^(1/6) * sigma;
U_eq = -epsilon;

%% Create visualization
figure('Position', [100, 100, 1000, 400]);

% Plot potential
subplot(1, 2, 1);
plot(r/sigma, U/epsilon, 'b-', 'LineWidth', 2);
hold on;
plot([r_eq/sigma, r_eq/sigma], [min(U)/epsilon, 0], 'r--', 'LineWidth', 1.5);
plot(r_eq/sigma, U_eq/epsilon, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot([r_min/sigma, r_max/sigma], [0, 0], 'k--', 'LineWidth', 0.5);
hold off;
xlabel('Distance r/\sigma');
ylabel('Potential Energy U/\epsilon');
title('Lennard-Jones Potential');
grid on;
legend('U(r)', 'Equilibrium', sprintf('r_{eq} = %.3f\\sigma', r_eq/sigma), ...
       'Location', 'best');
ylim([min(U)/epsilon*1.2, 2]);
text(1.5, 0.5, 'Repulsive region', 'FontSize', 10);
text(2.2, -0.5, 'Attractive region', 'FontSize', 10);

% Plot force
subplot(1, 2, 2);
plot(r/sigma, F/epsilon, 'r-', 'LineWidth', 2);
hold on;
plot([r_eq/sigma, r_eq/sigma], [min(F)/epsilon, max(F)/epsilon], 'b--', 'LineWidth', 1.5);
plot([r_min/sigma, r_max/sigma], [0, 0], 'k--', 'LineWidth', 0.5);
hold off;
xlabel('Distance r/\sigma');
ylabel('Force F\sigma/\epsilon');
title('Lennard-Jones Force');
grid on;
legend('F(r)', 'Equilibrium', 'Location', 'best');
text(1.2, max(F)/epsilon * 0.5, 'Repulsive', 'FontSize', 10);
text(2.0, min(F)/epsilon * 0.3, 'Attractive', 'FontSize', 10);

%% Print key values
fprintf('Lennard-Jones Potential Properties:\n');
fprintf('=====================================\n');
fprintf('Equilibrium distance: r_eq = %.4f * sigma = %.4f\n', r_eq/sigma, r_eq);
fprintf('Potential at equilibrium: U(r_eq) = %.4f * epsilon = %.4f\n', U_eq/epsilon, U_eq);
fprintf('Zero crossing (sigma): r = %.4f\n', sigma);
fprintf('\n');
fprintf('Physical Interpretation:\n');
fprintf('- For r < r_eq: Repulsive force (particles push apart)\n');
fprintf('- For r = r_eq: Zero force (equilibrium position)\n');
fprintf('- For r > r_eq: Attractive force (particles pull together)\n');
fprintf('- For r >> r_eq: Weak attraction (van der Waals)\n');
fprintf('\n');
fprintf('In the simulation:\n');
fprintf('- Particles tend to maintain separation near r_eq\n');
fprintf('- Too close → strong repulsion (like hard spheres)\n');
fprintf('- Too far → weak attraction (like gases)\n');
