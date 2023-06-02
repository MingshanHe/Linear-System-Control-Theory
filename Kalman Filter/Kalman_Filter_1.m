%%
A = [1.1269 -0.4940 0.1129; 1 0 0; 0 1 0];
B = [-0.3832; 0.5919; 0.5191];
C = [1 0 0];
D = 0;
Plant = ss(A,[B B], C, D, -1, 'inputname',{'u','w'},'outputname','y');

Q = 2.3;
R = 1;

[kalmf, L, ~, M, zeros] = kalman(Plant, Q, R);
kalmf = kalmf(1,:);

a = A;
b = [B B 0*B];
c = [C; C];
d = [0 0 0; 0 0 1];
P = ss(a, b, c, d, -1, 'inputname',{'u','v','w'},'outputname',{'y', 'yv'});

sys = parallel(P, kalmf, 1, 1, [], []);
SimModel = feedback(sys, 1, 4, 2, 1);
SimModel = SimModel([1, 3], [1, 2, 3]);

t = (0:200)';
u = sin(t/5);

rng(10, 'twister');
w = sqrt(Q)*randn(length(t), 1);
v = sqrt(R)*randn(length(t), 1);

out = lsim(SimModel, [w, v, u]);
y = out(:, 1);
ye = out(:, 2);
yv = y + v;
clf;
plot(t,y,'b',t,ye,'r--');

%%
clear;
clc;
% System matrices
A = [1 1; 0 1];                 % State transition matrix
B = [0.5; 1];                   % Input matrix
C = [1 0];                      % Measurement matrix
D = 0;                          % Feedthrough matrix
Q = [0.01 0; 0 0.01];           % Process noise covariance
R = 1;                          % Measurement noise covariance

% Initial state estimate
x0 = [0; 0];                    % Initial state estimate
P0 = eye(2);                    % Initial error covariance

% Generate input and measurement data
T = 100;                        % Number of time steps
u = randn(1, T);                % Input sequence
y = zeros(1, T);                % Measurements
x_true = zeros(2, T);           % True state
x_estimate = zeros(2, T);       % Estimated state

% Initialize the state estimate and error covariance
x_estimate(:, 1) = x0;
P = P0;

% Kalman filter loop
for t = 1:T
    % Predict the next state
    x_predicted = A * x_estimate(:, t) + B * u(:, t);
    P_predicted = A * P * A' + Q;
    
    % Update the state estimate based on the measurement
    K = P_predicted * C' / (C * P_predicted * C' + R);
    x_estimate(:, t+1) = x_predicted + K * (y(:, t) - C * x_predicted);
    P = (eye(2) - K * C) * P_predicted;
    
    % Generate process noise
    w = mvnrnd([0; 0], Q)';
    
    % Generate measurement noise
    v = mvnrnd(0, R);
    
    % True state update
    x_true(:, t+1) = A * x_true(:, t) + B * u(:, t) + w;
    
    % Simulate measurement
    y(:, t) = C * x_true(:, t) + v;
end

% Plot the true state and estimated state
t = 1:T;
figure;
subplot(2, 1, 1);
plot(t, x_true(1, 1:T), 'b-', t, x_estimate(1, 1:T), 'r--');
legend('True State', 'Estimated State');
xlabel('Time');
ylabel('X1');
title('Kalman Filter: X1');

subplot(2, 1, 2);
plot(t, x_true(2, 1:T), 'b-', t, x_estimate(2, 1:T), 'r--');
legend('True State', 'Estimated State');
xlabel('Time');
ylabel('X2');
title('Kalman Filter: X2');

