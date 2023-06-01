%% Parameter Definition
mc = 1;
mp = 1;
l = 1;
g = 10;

A = [0 0 1 0; 0 0 0 1; 0 -mp*g/mc 0 0; 0 -(mc+mp)*g/(l*mc) 0 0];
B = [0; 0; 1/mc; 1/(l*mc)];
C = [1 0 0 0; 0 1 0 0];
D = 0;


plot_ = false;
if plot_
    figure;
    t_max = 40;
    subplot(2,1,1)
    x0 = [0; 0.01; 0; 0];
    sys = ss(A, B, C, D);
    [y,t] = initial(sys, x0);
    plot(t(1:t_max),y((1:t_max),1));
    xlabel('time')
    ylabel('x(position)')
    subplot(2,1,2)
    plot(t(1:t_max),y((1:t_max),2))
end

Q = diag([1,1,1,1]);
R = 1;
K = lqr(A, B, Q, R);
% Define the initial conditions
x0 = [0; 0.01; 0; 0];   % Initial state
u = 0;                  % Initial Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(4, length(t));
y = zeros(2, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 


for i = 2:length(t)
    u = -K * x(:, i-1);
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    y(:, i) = C*x(:, i);
end

plot_ = true;
if plot_
    figure;
    % subplot(2, 1, 1);
    % plot(t, x(1, :), 'b', t, x(2,:), 'r');
    % ylabel('State');
    % legend('True State: x', 'True Output: theta');

    % subplot(2, 1, 2);
    plot(t, y(1,:), 'b', t, y(2,:), 'r');
    xlabel('Time');
    ylabel('Output');
    legend('True Output: x', 'True Output: theta');
end