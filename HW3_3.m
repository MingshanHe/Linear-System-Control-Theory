%% Parameter Definition
A = [0 1 0; 0 0 1; 1 0 0];
B = [0.3; 0; -0.3];
C = [1.9 1.3 1];
D = 0;

% Define the initial conditions
x0 = [0; 0; 0];        % Initial state
u = 0;                  % Input

% Define the simulation parameters
dt = 1;
t_end = 200;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));
x_Noise = zeros(3, length(t));
y_Noise = zeros(1, length(t));

x(:, 1) = x0;
x_Noise(:, 1) = x0 + processGaussianNoise(1);
y(:, 1) = C*x0; 
y_Noise(:, 1) = C*x0 + measurementGaussianNoise(1);
% Parameters
meanValue = 0;  % Mean value
stdDeviation = 0.7;  % Standard deviation
measurementGaussianNoise = meanValue + stdDeviation * randn(t_end+1, 1);
meanVector = [0; 0; 0];
covarianceMatrix = [0.004 0.002 0.001; 0.002 0.004 0.000; 0.001 0.000 0.001];
processGaussianNoise = meanVector + chol(covarianceMatrix)' * randn(3, t_end+1);

for i = 2:length(t)
    if(i<=100)
        u = 1;
    else
        u = -1;
    end
    x(:, i) = (A*x(:, i-1) + B*u);
    x_Noise(:, i) = (A*x(:, i-1) + B*u + processGaussianNoise(i));
    y(:, i) = C*x(:, i-1);
    y_Noise(:, i) = C*x(:, i-1) + measurementGaussianNoise(i);
end

plot_ = true;
if plot_
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), 'b');
    hold on;
    plot(t, x_Noise(1, :), 'r');
    ylabel('State 1');
    legend('True State w/o Noise', 'True State w Noise');
    
    subplot(3, 1, 2);
    plot(t, x(2, :), 'b');
    hold on;
    plot(t, x_Noise(2, :), 'r');
    xlabel('Time');
    ylabel('State 2');
    legend('True State w/o Noise', 'True State w Noise');

    subplot(3, 1, 3);
    plot(t, x(3, :), 'b');
    hold on;
    plot(t, x_Noise(3, :), 'r');
    xlabel('Time');
    ylabel('State 3');
    legend('True State w/o Noise', 'True State w Noise');
end

plot_ = true;
if plot_
    figure;
    plot(t, y(1, :), 'b');
    ylabel('Output');
    legend('True Output w/o Noise');
    hold on;
    plot(t, y_Noise(1, :), 'r');
    ylabel('Output');
    legend('True Output w Noise');
end