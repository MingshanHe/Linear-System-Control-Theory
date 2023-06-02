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
% Mainly for the State Feedback and Track Controller
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

% Define the initial conditions
x0 = [1; 1; 1];         % Initial state
x_estimate0 = [0; 0; 0];% Initial estimated state
u = 0;                  % Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(3, length(t));

x_estimate = zeros(3, length(t));
y_estimate = zeros(3, length(t));

x(:, 1) = x0;
x_estimate(:, 1) = x_estimate0;


for i = 2:length(t)
    % Update the state using the system equations
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    % Update the estimated state using the observer equations
    x_estimate(:, i) = x_estimate(:, i-1) + dt * (A*x_estimate(:, i-1) + B*u);
end

% Plot the results
plot_ = true;
if plot_
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), 'b');
    hold on;
    plot(t, x_estimate(1, :), 'r--');
    ylabel('State 1');
    legend('True State', 'Estimated State','Feedback State');

    subplot(3, 1, 2);
    plot(t, x(2, :), 'b');
    hold on;
    plot(t, x_estimate(2, :), 'r--');
    ylabel('State 2');
    legend('True State', 'Estimated State','Feedback State');

    subplot(3, 1, 3);
    plot(t, x(3, :), 'b');
    hold on;
    plot(t, x_estimate(3, :), 'r--');
    ylabel('State 3');
    legend('True State', 'Estimated State','Feedback State');
end

%% Closed-Loop Observer
% Compute Observer Gain K w/o 'place Funtion
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
% Compute State Feedback Gain L w/o place Function
C_ = [C' A'*C' A'*A'*C'];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
L_bar = [8 19 5];
L = (L_bar*P)';
% Compute State Feedback Gain L and K w 'place Funtion
desired_poles = [-0.01 -0.02 -0.03];
L_placed = place(A', C', desired_poles)';
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

display_ = true;
if display_
    disp("Computed By 'place' function: ")
    disp(L_placed);
    disp("Computed without 'place' function: ")
    disp(L);
end

% Define the initial conditions
x0 = [1; 1; 1];         % Initial state
x_estimate0 = [0; 0; 0];     % Initial estimated state
u = 0;                  % Input

% Define the simulation parameters
dt = 0.01;              % Time step
t_end = 10;             % Simulation end time

% Simulate the system and observer dynamics
t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));
x_estimate = zeros(3, length(t));
y_estimate = zeros(1, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 
x_estimate(:, 1) = x_estimate0;
y_estimate(:, 1) = C*x_estimate0;

for i = 2:length(t)
    u = -K_placed * x_estimate(:, i-1);
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u);
    y(i) = C*x(:, i);
    x_estimate(:, i) = x_estimate(:, i-1) + dt * ((A - L*C)*x_estimate(:, i-1) + B*u + L*y(i));
    y_estimate(i) = C*x_estimate(:, i);
end

plot_ = true;
if plot_
    figure;
    subplot(2, 2, 1);
    plot(t, x(1, :), 'b', t, x_estimate(1, :), 'r--');
    ylabel('State 1');
    legend('True State', 'Estimated State');
    
    subplot(2, 2, 2);
    plot(t, x(2, :), 'b', t, x_estimate(2, :), 'r--');
    xlabel('Time');
    ylabel('State 2');
    legend('True State', 'Estimated State');

    subplot(2, 2, 3);
    plot(t, x(3, :), 'b', t, x_estimate(3, :), 'r--');
    xlabel('Time');
    ylabel('State 3');
    legend('True State', 'Estimated State');

    subplot(2,2,4);
    plot(t, y,'b', t, y_estimate,'r--');
    xlabel('Time');
    ylabel('Output');
    legend('True Output', 'Estimated Output');
end

%% Reference Tracking Controller
% Compute Observer Gain K w/o 'place Funtion
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
% Compute State Feedback Gain L and K w 'place Funtion
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);
dt = 0.01;              % Time step
t_end = 10;   
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
t_end = 10;   
K_Ref = [17 35 35];
R = 1;

t = 0:dt:t_end;
x = zeros(3, length(t));
y = zeros(1, length(t));

x(:, 1) = x0;
y(:, 1) = C*x0; 
for i = 2:length(t)
    x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*(K_Ref*x(:, i-1)-85*R));
    y(:, i) = C*x(:, i); 
end
plot_ = true;
if plot_
    figure;
    plot(t, y);
    xlabel('Time');
    ylabel('Output');
    title('Reference Tracking Controller')
end

%% Robust Controller
A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
Ad = [0, -2; 2, 0];
Cd = [1, 0];
A_aug = [A B*Cd; zeros(2,3) Ad];
B_aug = [B; zeros(2,1)];
C_aug = [C zeros(1,2)];
% Compute State Feedback Gain L and K w 'place Funtion
desired_poles = [-4 -4+0.5i -4-0.5i -4+2i -4-2i];
L_placed = place(A_aug', C_aug', desired_poles)';
desired_poles=[-5 -5+3i -5-3i];
K_placed = place(A, B, desired_poles);

R = 1;
dt = 0.01
t = 0:dt:10;
x = zeros(5, length(t));
x_estimate = zeros(5, length(t));
y = zeros(1, length(t));
y_estimate = zeros(1, length(t));

% Define the initial conditions
x(:, 1) = [1; 1; 1; 5; 5];
x_estimate(:,1) = [0; 0; 0; 5; 5];
y(:, 1) = C_aug*x(:, 1) ; 
y_estimate(:,1) = C_aug*x_estimate(:,1);

for i = 2:length(t)
    u = -K_placed * x_estimate(1:3, i-1);
    x(:, i) = x(:, i-1) + dt * (A_aug*x(:,i-1)+ B_aug*u);
    y(:, i) = C_aug*x(:, i); 
    x_estimate(:, i) = x_estimate(:, i-1) + dt * (A_aug*x_estimate(:, i-1) + B_aug*u + L_placed*(y(i-1)-C_aug*x_estimate(:,i-1)));
end

plot_ = true;
if plot_
    figure;
    plot(t, y);
    xlabel('Time');
    ylabel('Output');
end
plot_ = true;
if plot_
    figure;
    subplot(3, 1, 1);
    plot(t, x(1, :), 'b', t, x_estimate(1, :), 'r--');
    ylabel('State 1');
    legend('True State', 'Estimated State');

    subplot(3, 1, 2);
    plot(t, x(2, :), 'b', t, x_estimate(2, :), 'r--');
    xlabel('Time');
    ylabel('State 2');
    legend('True State', 'Estimated State');

    subplot(3, 1, 3);
    plot(t, x(3, :), 'b', t, x_estimate(3, :), 'r--');
    xlabel('Time');
    ylabel('State 3');
    legend('True State', 'Estimated State');
end

% % Simulate the system and observer dynamics
% t = 0:dt:t_end;
% x = zeros(3, length(t));
% y = zeros(1, length(t));
% x_estimate = zeros(3, length(t));
% y_estimate = zeros(1, length(t));
% 
% % Compute State Feedback Gain K w 'place Funtion
% desired_poles = [-0.1 -0.2 -0.3];
% L_placed = place(A', C', desired_poles)';
% desired_poles=[-10 -10+3i -10-3i];
% K_placed = place(A, B, desired_poles);
% 
% x(:, 1) = x0;
% y(:, 1) = C*x0; 
% x_estimate(:, 1) = x_estimate0;
% 
% 
% for i = 2:length(t)
%     u = -K_placed * x_estimate(:, i-1);
%     x(:, i) = x(:, i-1) + dt * (A*x(:, i-1) + B*u + B*(K_Ref*x(:, i-1))+xi*sin(omega*t(i-1)));
%     y(:, i) = C*x(:, i); 
%     x_estimate(:, i) = x_estimate(:, i-1) + dt * ((A - L_placed*C)*x_estimate(:, i-1) + B*u + L_placed*y(:, i));
%     y_estimate(:,i) = C*x_estimate(:,i); 
% end
% plot_ = false;
% if plot_
%     figure;
%     subplot(2, 2, 1);
%     plot(t, x(1, :), 'b', t, x_estimate(1, :), 'r--');
%     ylabel('State 1');
%     legend('True State', 'Estimated State');
% 
%     subplot(2, 2, 2);
%     plot(t, x(2, :), 'b', t, x_estimate(2, :), 'r--');
%     xlabel('Time');
%     ylabel('State 2');
%     legend('True State', 'Estimated State');
% 
%     subplot(2, 2, 3);
%     plot(t, x(3, :), 'b', t, x_estimate(3, :), 'r--');
%     xlabel('Time');
%     ylabel('State 3');
%     legend('True State', 'Estimated State');
% 
%     subplot(2,2,4);
%     plot(t, y,'b', t, y_estimate,'r--');
%     xlabel('Time');
%     ylabel('Output');
%     legend('True Output', 'Estimated Output');
% end