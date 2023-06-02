% Define the system matrices
A = [0, -0.2; 0.2, 0];  % State matrix
B = [0; 0];  % Input matrix
C = [1, 0];  % Output matrix
D = 0;  % Feedforward matrix

% Define the initial condition: x0 = [1; -1]
x0 = [1; 1];

% Define the time span
tspan = 0:0.1:50;  % Time vector

% Compute the state and output response
sys = ss(A, B, C, D);  % Create the state space system object
[y, t, x] = lsim(sys, zeros(size(tspan)), tspan, x0);

% Plot the state variables
figure;
subplot(2, 1, 1);
plot(t, x(:, 1), 'b');
xlabel('Time');
ylabel('x1');
title('State Variable x1');

subplot(2, 1, 2);
plot(t, x(:, 2), 'r');
xlabel('Time');
ylabel('x2');
title('State Variable x2');
% Plot the output response
figure;
plot(tspan, y);
xlabel('Time');
ylabel('Output (y)');
title('Output Response');
grid on;

