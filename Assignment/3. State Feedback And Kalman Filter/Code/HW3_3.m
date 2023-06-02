%% Parameter Definition
A = [0 1 0; 0 0 1; 1 0 0];
B = [0.3; 0; -0.3];
C = [1.9 1.3 1];
D = 0;
% Define the initial conditions
x0 = [0; 0; 0];         % Initial state
u = 0;                  % Input

% Define the simulation parameters
dt = 1;                  % 1：Discrete Time
t_end = 200;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));
x_true = zeros(3, length(t));
y_true = zeros(1, length(t));
x_estimate = zeros(3, length(t));
y_estimate = zeros(1, length(t));
% Parameters
meanValue = 0;  % Mean value
stdDeviation = 0.7;  % Standard deviation
measurementGaussianNoise = meanValue + stdDeviation * randn(t_end+1, 1);
meanVector = [0; 0; 0];
covarianceMatrix = [0.004 0.002 0.001; 0.002 0.004 0.000; 0.001 0.000 0.001];
processGaussianNoise = meanVector + chol(covarianceMatrix)' * randn(3, t_end+1);

x(:, 1) = x0;
x_true(:, 1) = x0 + processGaussianNoise(1);
y(:, 1) = C*x0; 
y_true(:, 1) = C*x0 + measurementGaussianNoise(1);
x_estimate(:, 1) = [0;0;0];
y_estimate(:, 1) = C*x_estimate(:, 1); 

Q = covarianceMatrix; % Process Noise Covariance
R = stdDeviation; % Measurement Noise Covariance
P = eye(3);

for i = 1:length(t)-1
    if(i<=100)
        u = 1;
    else
        u = -1;
    end

    % Predict the next state
    x_predicted = A*x_estimate(:,i) + B*u;
    P_predicted = A*P*A'+Q;
    % Update the state estimate based on the measurement
    K = P_predicted * C' / (C*P_predicted*C' + R);
    x_estimate(:, i+1) = x_predicted + K * (y_true(:, i) - C * x_predicted);
    P = (eye(3) - K * C) * P_predicted;   
    % True state update
    x(:, i+1) = (A*x(:, i) + B*u);
    x_true(:, i+1) = (A*x(:, i) + B*u + processGaussianNoise(i));
    % Simulate measurement
    y(:, i) = C*x(:, i);
    y_true(:, i) = C*x(:, i) + measurementGaussianNoise(i);
    % Store the estimated state
end

plot_ = true;
if plot_
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), 'b');
    hold on;
    plot(t, x_true(1, :), 'r');
    hold on;
    plot(t, x_estimate(1,:), 'g');
    ylabel('State 1');
    legend('True State w/o Noise', 'True State w Noise', 'Estimated State w Noise');
    
    subplot(3, 1, 2);
    plot(t, x(2, :), 'b');
    hold on;
    plot(t, x_true(2, :), 'r');
    hold on;
    plot(t, x_estimate(1,:), 'g');
    xlabel('Time');
    ylabel('State 2');
    legend('True State w/o Noise', 'True State w Noise', 'Estimated State w Noise');

    subplot(3, 1, 3);
    plot(t, x(3, :), 'b');
    hold on;
    plot(t, x_true(3, :), 'r');
    hold on;
    plot(t, x_estimate(1,:), 'g');
    xlabel('Time');
    ylabel('State 3');
    legend('True State w/o Noise', 'True State w Noise', 'Estimated State w Noise');
end

plot_ = false;
if plot_
    figure;
    plot(t, y(1, :), 'b');
    ylabel('Output');
    legend('True Output w/o Noise');
    hold on;
    plot(t, y_true(1, :), 'r');
    hold on;
    plot(t, y_hat(1, :), 'g');
    ylabel('Output');
    legend('True Output w Noise');
end