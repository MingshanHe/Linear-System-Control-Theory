% System parameters
m = 1; % mass (kg)
c = 4; % damping coefficient (N s/m)
k = 10; % spring constant (N/m)

% State space equation
A = [0 1; -k/m -c/m];
B = [0; 1/m];
C = [1 0];
D = 0;
sys = ss(A,B,C,D);

% Initial conditions
x0 = [0.1; 0]; % initial position and velocity

% Simulation time
tspan = 0:0.01:10;

% Input
% u = 10000*ones(size(tspan)); % no input force
u = zeros(size(tspan));

% Simulate system response
[y, t, x] = lsim(sys, u, tspan, x0);

% Plot response
plot(t, y);
xlabel('Time (s)');
ylabel('Position (m)');
title('Mass-Spring-Damper System Response');
