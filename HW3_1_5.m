% System matrices
A = [0.8 0.2; -0.3 0.9];  % State matrix
C = [1 0; 0 1];           % Measurement matrix

% Initial conditions
x0 = [0.5; -0.5];  % Initial state estimate

% Generate input and measurement data
T = 100;   % Number of time steps
u = zeros(2, T);  % Input sequence (zero input)
y = zeros(2, T);  % Measurements
x_true = zeros(2, T);  % True state
x_estimate = zeros(2, T);  % Estimated state

% Initialize the state estimate
x_estimate(:, 1) = x0;

% Open-loop observer loop
for t = 1:T
    % Update the state estimate based on the system dynamics
    x_estimate(:, t+1) = A * x_estimate(:, t);
    
    % Generate measurement noise
    v = mvnrnd([0; 0], eye(2))';
    
    % Simulate measurement
    y(:, t) = C * x_true(:, t) + v;
    
    % Update the true state based on the system dynamics
    x_true(:, t+1) = A * x_true(:, t);
end

% Plot the true state and estimated state
t = 0:T-1;
figure;
subplot(3, 1, 1);
plot(t, x_true(1, 1:100), 'b-', t, x_estimate(1, 1:T), 'r--');
legend('True State', 'Estimated State');
xlabel('Time');
ylabel('X1');
title('Open-Loop Observer: X1');

subplot(3, 1, 2);
plot(t, x_true(2, 1:100), 'b-', t, x_estimate(2, 1:T), 'r--');
legend('True State', 'Estimated State');
xlabel('Time');
ylabel('X2');
title('Open-Loop Observer: X2');

subplot(3,1,3);
plot(t,y);
