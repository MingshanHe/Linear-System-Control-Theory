% Define the system matrices
A = [1 1; 0 1];  % State transition matrix
B = [0.5; 1];    % Control input matrix
C = [1 0];       % Measurement matrix

% Define the process noise and measurement noise covariance
Q = [0.1 0; 0 0.1];  % Process noise covariance
R = 1;               % Measurement noise covariance

% Initialize the state estimate and error covariance
x_hat = [0; 0];     % Initial state estimate
P = eye(2);         % Initial error covariance

% Initialize the true state
x_true = [0; 0];   % Initial true state

% Simulate the system and perform Kalman filtering
T = 100;           % Number of time steps
y = zeros(1, T);   % Measurement vector
x_est = zeros(2, T);  % Estimated state vector

for t = 1:T
    % Generate the true state and measurement
    x_true = A * x_true + B .* sqrtm(Q) * randn(2,1);
    y(t) = C * x_true + sqrt(R) * randn;
    
    % Kalman filter update
    y_pred = C * x_hat;
    e = y(t) - y_pred;
    K = P * C' / (C * P * C' + R);
    x_hat = A * x_hat + K * e;
    P = A * P * A' + Q - A * P * C' / (C * P * C' + R) * C * P * A';
    
    % Store the estimated state
    x_est(:, t) = x_hat;
end

% Plot the true state, measurement, and estimated state
t = 1:T;
figure;
subplot(2, 1, 1);
plot(t, x_est(1, :), 'r-', t, y, 'b-', t, x_true(1, :), 'g--');
legend('Estimated state', 'Measurement', 'True state');
xlabel('Time');
ylabel('State');
title('State 1');

subplot(2, 1, 2);
plot(t, x_est(2, :), 'r-', t, x_true(2, :), 'g--');
legend('Estimated state', 'True state');
xlabel('Time');
ylabel('State');
title('State 2');
