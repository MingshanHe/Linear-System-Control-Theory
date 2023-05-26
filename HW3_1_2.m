% Define the continuous-time state space equations
A = [1 2; -1 0];
B = [1; 0];
C = [1 0];
D = 0;

% Define the observer gain matrix
L = [3; 2];

% Define the control gain matrix
K = [1 0];

% Define the initial conditions
x0 = [1; 1];         % Initial state
x_hat0 = [0; 0];     % Initial estimated state

% Define the simulation parameters
dt = 0.01;           % Time step
t_end = 10;          % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(2, length(t));
x_hat = zeros(2, length(t));

x(:, 1) = x0;
x_hat(:, 1) = x_hat0;

for i = 2:length(t)
    % Update the state using the system equations
    u = -K * x_hat(:, i-1);  % Control law based on the estimated state
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    
    % Update the estimated state using the observer equations
    y = C*x(:, i);
    x_hat(:, i) = x_hat(:, i-1) + dt * ((A - L*C)*x_hat(:, i-1) + B*u + L*y);
end

% Plot the results
figure;
subplot(2, 1, 1);
plot(t, x(1, :), 'b', t, x_hat(1, :), 'r--');
ylabel('State 1');
legend('True State', 'Estimated State');

subplot(2, 1, 2);
plot(t, x(2, :), 'b', t, x_hat(2, :), 'r--');
xlabel('Time');
ylabel('State 2');
legend('True State', 'Estimated State');
