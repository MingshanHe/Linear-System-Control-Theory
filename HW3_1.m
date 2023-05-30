%%% HW3_1
% Mainly for the State Feedback and Track Controller
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
% Define State Space Equation System
sys = ss(A, B, C, D);
eigenvals = eig(A);
char_poly = poly(eigenvals);
% Compute State Feedback Gain K w/o 'place Funtion
C_ = [B A*B A*A*B];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
K_bar = [17 92 169];
K = K_bar*P;
% Compute State Feedback Gain K w 'place Funtion
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

display_ = false;
if display_
    disp("Computed By 'place' function: ")
    disp(K_placed);
    disp("Computed without 'place' function: ")
    disp(K);
end
plot_ = false;
if plot_
    sys = ss(A, B, C, D);
    A_ = A-B*K;
    x0 = [1; 1; 1];
    sys_ = ss(A_, B, C, D);
    [y,t] = initial(sys, x0);
    plot(t(1:20),y(1:20));
    hold on;
    [y_,t] = initial(sys_, x0);
    plot(t,y_);
    xlabel('Time');
    ylabel('Output');
    legend('Original System', 'State Feedback System');
end

%% Open-Loop Observer
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

% Define the initial conditions
x0 = [1; 1; 1];         % Initial state
x_hat0 = [0; 0; 0];     % Initial estimated state
u = 0;                  % Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(3, length(t));

x_hat = zeros(3, length(t));
y_hat = zeros(3, length(t));

x_ = zeros(3, length(t));
y_ = zeros(3, length(t));

x(:, 1) = x0;
y(:, 1) = A*x0 + B*u; 
x_hat(:, 1) = x_hat0;
y_hat(:, 1) = A*x_hat0+B*u;
x_(:,1) = x0;


for i = 2:length(t)
    % Update the state using the system equations
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    y(:, i) = C*x(:, i); 
    % Update the estimated state using the observer equations
    x_hat(:, i) = x_hat(:, i-1) + dt * (A*x_hat(:, i-1) + B*u);
    y_hat(:, i) = C*x_hat(:, i); 
    % State Feedback with Open-Loop Observer
    % x_(:, i) = x_(:, i-1) + dt * (A*x_(:, i-1)-(B*K_placed)*x_(:, i-1) + B*u);
    % y_(:, i) = C*x_(:, i); 
end

% Plot the results
plot_ = true;
if plot_
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), 'b');
    hold on;
    plot(t, x_hat(1, :), 'r--');

    ylabel('State 1');
    legend('True State', 'Estimated State','Feedback State');

    subplot(3, 1, 2);
    plot(t, x(2, :), 'b');
    hold on;
    plot(t, x_hat(2, :), 'r--');
    ylabel('State 2');
    legend('True State', 'Estimated State','Feedback State');

    subplot(3, 1, 3);
    plot(t, x(3, :), 'b');
    hold on;
    plot(t, x_hat(3, :), 'r--');
    ylabel('State 3');
    legend('True State', 'Estimated State','Feedback State');
end

%% Closed-Loop Observer
% Compute Observer Gain K w/o 'place Funtion
C_ = [C' A'*C' A'*A'*C'];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
L_bar = [8 19 5];
L = (L_bar*P)';
% Compute State Feedback Gain K w 'place Funtion
desired_poles = [-1 -2 -3];
L_placed = place(A', C', desired_poles)';

display_ = true;
if display_
    disp("Computed By 'place' function: ")
    disp(L_placed);
    disp("Computed without 'place' function: ")
    disp(L);
end

% Define the initial conditions
x0 = [1; 1; 1];         % Initial state
x_hat0 = [0; 0; 0];     % Initial estimated state
u = 0;                  % Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(3, length(t));
x_hat = zeros(3, length(t));
y_hat = zeros(3, length(t));
x_ = zeros(3, length(t));
y_ = zeros(3, length(t));

x(:, 1) = x0;
% y(:, 1) = A*x0 + B*u; 
x_hat(:, 1) = x_hat0;
% y_hat(:, 1) = A*x_hat0+B*u;

for i = 2:length(t)
    u = -K * x_hat(:, i-1);
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    y = C*x(:, i);
    x_hat(:, i) = x_hat(:, i-1) + dt * ((A - L*C)*x_hat(:, i-1) + B*u + L*y);
end

plot_ = true;
if plot_
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
end

%% Reference Tracking Controller
dt = 0.01;              % Time step
t_end = 20;   
K_Ref = [17 35 35];
Ka_Ref = 10;
R = 1;
A_Ref = [A+B*K_Ref, B*Ka_Ref; -C, 0];
B_Ref = [zeros(3,1); 1];
C_Ref = [C 0];
t = 0:dt:t_end;

x = zeros(3, length(t));
y = zeros(1, length(t));
xa = zeros(1, length(t));
x(:, 1) = x0;
y(:, 1) = C*x0; 
xa(:,1) = 0;
for i = 2:length(t)
    X = [x(:, i-1); xa(:, i-1)] + dt * (A_Ref*[x(:, i-1); xa(:, i-1)] + B_Ref*R);
    x(:, i) = X(1:3);
    xa(:, i) = X(4);
    y(:, i) = C_Ref*[x(:, i); xa(:, i)]; 
end
plot_ = false;
if plot_
    figure;
    plot(t, y);
    xlabel('Time');
    ylabel('Output');
end


dt = 0.01;              % Time step
t_end = 20;   
K_Ref = [17 35 35];
Ka_Ref = 10;
R = 1;
A_Ref = A;
B_Ref = B;
C_Ref = C;
t = 0:dt:t_end;

x = zeros(3, length(t));
y = zeros(1, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 
for i = 2:length(t)
    x(:, i) = x(:, i-1) + dt * (A_Ref*x(:, i-1) + B_Ref*(K_Ref*x(:, i-1)-85*R));
    y(:, i) = C_Ref*x(:, i); 
end
plot_ = false;
if plot_
    figure;
    plot(t, y);
    xlabel('Time');
    ylabel('Output');
end

%% Robust Controller
dt = 0.01;              % Time step
t_end = 20;   
xi = 1;
omega = 1;
K_Ref = [17 35 35];
Ka_Ref = 10;
R = 1;
A_Rob = A;
B_Rob = B;
C_Rob = C;
t = 0:dt:t_end;

x = zeros(3, length(t));
y = zeros(1, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 
for i = 2:length(t)
    x(:, i) = x(:, i-1) + dt * (A_Rob*x(:, i-1) + B_Rob*(K_Ref*x(:, i-1))+xi*sin(omega*t(i-1)));
    y(:, i) = C_Ref*x(:, i); 
end
plot_ = false;
if plot_
    figure;
    plot(t, y);
    xlabel('Time');
    ylabel('Output');
end

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));
x_hat = zeros(3, length(t));
y_hat = zeros(3, length(t));
x_ = zeros(3, length(t));
y_ = zeros(3, length(t));

% Compute Observer Gain K w/o 'place Funtion
C_ = [C' A'*C' A'*A'*C'];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
L_bar = [8 19 5];
L = (L_bar*P)';
% Compute State Feedback Gain K w 'place Funtion
desired_poles = [-2 -2+0.5i -2-0.5i];
L_placed = place(A', C', desired_poles)';

x(:, 1) = x0;
y(:, 1) = C_Rob*x0; 
x_hat(:, 1) = x_hat0;


for i = 2:length(t)
    u = -K * x_hat(:, i-1);
    x(:, i) = x(:, i-1) + dt * (A_Rob*x(:, i-1) + B_Rob*u + B_Rob*(K_Ref*x(:, i-1))+xi*sin(omega*t(i-1)));
    y(:, i) = C_Rob*x(:, i); 
    x_hat(:, i) = x_hat(:, i-1) + dt * ((A_Rob - L_placed*C_Rob)*x_hat(:, i-1) + B_Rob*u + L_placed*y(:, i));
end
plot_ = true;
if plot_
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), 'b', t, x_hat(1, :), 'r--');
    ylabel('State 1');
    legend('True State', 'Estimated State');
    
    subplot(3, 1, 2);
    plot(t, x(2, :), 'b', t, x_hat(2, :), 'r--');
    xlabel('Time');
    ylabel('State 2');
    legend('True State', 'Estimated State');

    subplot(3, 1, 3);
    plot(t, x(3, :), 'b', t, x_hat(3, :), 'r--');
    xlabel('Time');
    ylabel('State 3');
    legend('True State', 'Estimated State');
end