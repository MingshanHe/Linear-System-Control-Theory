A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
% Define State Space Equation System
sys = ss(A, B, C, D);
eigenvals = eig(A);
char_poly = poly(eigenvals);
% Compute State Feedback Gain K w/o 'place Funtion
C = [B A*B A*A*B];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C*C_bar);
K_bar = [17 92 169];
K = K_bar*P;
% Compute State Feedback Gain K w 'place Funtion
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

display_ = true;
if display_
    disp("Computed By 'place' function: ")
    disp(K_placed);
    disp("Computed without 'place' function: ")
    disp(K);
end
responce_ = true;
if responce_
    A_ = A-B*K;
    x0 = [1; 1; 1];
    sys_ = ss(A_, B, C, D);
    [y,t] = initial(sys_, x0);
    plot(t,y)

    % sys_ = ss(A_, B, C, D);
    % step(sys_);
end

% Open-Loop Observer

desired_poles = [-1 -2 -3];
L = place(A', C', desired_poles)';

x0 = [0;0;0];
x_hat0 = [1;1;1];

% Set simulation parameters
t = 0:0.01:10;          % Time vector
u = zeros(size(t));    % Input vector (zero input)
x = zeros(3, length(t));       % Actual state vector
x_hat = zeros(3, length(t));   % Estimated state vector

% Simulate the closed-loop system
for i = 1:length(t)
    if i == 1
        x(:, i) = x0;    % Set initial state of the real system
        x_hat(:, i) = x_hat0;  % Set initial state of the observer
    else
        x(:, i) = A*x(:, i-1) + B*u(i-1);   % Update the actual state
        x_hat(:, i) = (A - L*C)*x_hat(:, i-1) + B*u(i-1) + L*C*x(:, i-1);   % Update the estimated state
    end
end



% Plot the response
figure;
subplot(3, 1, 1);
plot(t, x(1, :), 'r', 'LineWidth', 2);   % Plot the actual state 1
hold on;
plot(t, x_hat(1, :), '--r', 'LineWidth', 2);   % Plot the estimated state 1
hold on;
legend('Actual State 1', 'Estimated State 1');
subplot(3, 1, 2);
plot(t, x(2, :), 'b', 'LineWidth', 2);   % Plot the actual state 2
hold on;
plot(t, x_hat(2, :), '--b', 'LineWidth', 2);   % Plot the estimated state 2
hold on;
legend('Actual State 1', 'Estimated State 1');
subplot(3, 1, 3);
plot(t, x(2, :), 'b', 'LineWidth', 3);   % Plot the actual state 2
hold on;
plot(t, x_hat(2, :), '--b', 'LineWidth', 3);   % Plot the estimated state 2
hold on;
legend('Actual State 1', 'Estimated State 1');
xlabel('Time');
ylabel('State Value');
title('Closed-Loop Observer Response');

