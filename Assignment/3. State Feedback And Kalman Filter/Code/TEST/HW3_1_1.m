A = [0 1 0 0; 0 0 -1 0; 0 0 0 1; 0 0 5 0];
B = [0; 1; 0; -2];
C = [1 0 0 0];
D = 0;

eigenvalues = eig(A);
display_ = false;

C = [B A*B A*A*B A*A*A*B];
C_bar = [1 0 -5 0; 0 1 0 -5; 0 0 1 0; 0 0 0 1];
P = inv(C*C_bar);
k_bar = [5 15.5 11 5];
K = k_bar*P
p = [-1.5+0.5i -1.5-0.5i -1+i -1-i]
k = place(A,B,p)
if display_
     disp(eigenvalues);  
     disp(C);
end

% Define the system's state-space matrices
A = [1 1; 0 2];
B = [0; 1];
C = [1 0];
D = 0;

% Define the desired observer poles
desired_poles = [-2; -3];

% Compute the observer gain matrix using pole placement
L = place(A', C', desired_poles)';

% Compute the observer dynamics matrix
F = A - L*C;

eig(F)


% Define system matrices
A = [1 2; -1 0.5];    % System matrix
B = [1; 0];           % Input matrix
C = [1 0; 0 1];       % Output matrix

% Define desired observer poles
observer_poles = [-2, -3];

% Compute observer gain matrix
L = place(A', C', observer_poles).';  % Observer gain matrix

% Define initial conditions
x0 = [0; 0];           % Initial state of the real system
x_hat0 = [1; -1];      % Initial state of the observer

% Set simulation parameters
t = 0:0.01:10;          % Time vector
u = zeros(size(t));    % Input vector (zero input)
x = zeros(2, length(t));       % Actual state vector
x_hat = zeros(2, length(t));   % Estimated state vector

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
subplot(2, 1, 1);
plot(t, x(1, :), 'r', 'LineWidth', 2);   % Plot the actual state 1
hold on;
plot(t, x(2, :), 'b', 'LineWidth', 2);   % Plot the actual state 2
plot(t, x_hat(1, :), '--r', 'LineWidth', 2);   % Plot the estimated state 1
plot(t, x_hat(2, :), '--b', 'LineWidth', 2);   % Plot the estimated state 2
legend('Actual State 1', 'Actual State 2', 'Estimated State 1', 'Estimated State 2');
xlabel('Time');
ylabel('State Value');
title('Closed-Loop Observer Response');



