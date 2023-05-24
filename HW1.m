
%%

% HOMWORK1-1(A)
syms t t0 real
assume(t > t0 > 0)
A = [0 1/(1+t); -1/(1+t) 0];
Phi = expm(int(A,t0, t));
display(Phi);

%%

% HOMEWORK1-1(B)
syms t t0 real
assume(t0 > 5)
assume(t > 5)
A1 = [0 0; 0 1]; % for 0 <= t <2
A2 = [0 1; 0 1]; % for 2<= t
A = piecewise(t<2, A1, t>=2, A2);
Phi = expm(int(A, t0, t));
display(Phi);

%%

% HOMEWORK1-2(A)
syms t
h = (1+t);
L1_norm = int(h, 0, Inf);
if isfinite(L1_norm)
    disp('The system g(t) = 1/(1+t) is BIBO stable.');
else
    disp('The system g(t) = 1/(1+t) is not BIBO stable.');
end

%%

% HOMEWORK1-2(B)
% Define the system matrices
A = [-1 10; 0 1];
B = [-2; 0];
C = [-2 3];
D = 2;

sys = ss(A,B,C,D);


h = impulse(sys);

tf_sys = tf(sys);

w = linspace(0,10,1000);

H = freqresp(tf_sys,w);

mag = abs(H);

% Check if the system is BIBO stable
if max(mag(:)) < Inf
    disp('System is BIBO stable');
else
    disp('System is not BIBO stable');
end

%%

% HOMEWORK1-2(C)
% Define the system matrices
A = [-1 0 1; 0 0 1; 0 0 0];
B = [0; 0; 0];
C = [0 0 0];
D = 0;
sys = ss(A, B, C, D);

% Compute the system eigenvalues
e = eig(A);

% Check for marginal stability
if all(real(e) == 0)
    disp('System is marginally stable');
else
    disp('System is not marginally stable');
end


% Check for asymptotic stability
if all(real(e) < 0)
    disp('System is asymptotically stable');
else
    disp('System is not asymptotically stable');
end

%%

% HOMEWORK1-2(D)
% Define the LTI system
A = [-1 0 1; 0 0 1; 0 0 0];
B = [0; 0; 0];
C = [0 0 0];
D = 0;
sys = ss(A, B, C, D);

% Compute the Jordan form of A
J = jordan(A);
disp(J)
% Check if the system is marginally stable
is_marginally_stable = any(diag(J) == 0)

%%

% HOMEWORK1-2(E)
% Define the system matrices
A = [-1 0 1; 0 0 1; 0 0 0];
B = [0; 0; 0];
C = [0 0 0];
D = 0;
sys = ss(A, B, C, D);

% Compute the system eigenvalues
e = eig(A);

% Check for asymptotic stability
if all(real(e) < 0)
    disp('System is asymptotically stable');
else
    disp('System is not asymptotically stable');
end

%%

% HOMEWORK1-3(A)
% Define the system
A = [0.5 1; 0 -0.5];
x0 = [1; -1];

% Compute the zero-input response for the first state
n = 0:20;
x1 = zeros(2,length(n));
x1(:,1) = x0;
for i = 2:length(n)
    x1(:,i) = A*x1(:,i-1);
end

% Plot the two states in the same figure
figure;
stem(n,x1(1,:));
hold on;
stem(n,x1(2,:),'--');
legend('State 1', 'State 2');
xlabel('n');
ylabel('Amplitude');
title('Zero-Input Response of the System');



%%

% HOMEWORK1-3(B)
% Define the system
A = [0.5 1; 0 -0.5];
B = [0; 1];
x0 = [0; 0];

% Compute the zero-input response for the first state
n = 0:20;
x = zeros(2,length(n));
u = (9/2)*n;
x(:,1) = x0;
for i = 2:length(n)
    x(:,i) = A*x1(:,i-1)+B*u(:,i-1);
end

% Plot the two states in the same figure
figure;
stem(n,x(1,:));
hold on;
stem(n,x(2,:),'--');
legend('State 1', 'State 2');
xlabel('n');
ylabel('Amplitude');
title('Zero-State Response of the System');




