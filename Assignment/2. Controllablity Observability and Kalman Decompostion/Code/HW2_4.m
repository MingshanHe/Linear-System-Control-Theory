
num = [4 -2 -6];     % numerator coefficients of transfer function
den = [2 6 8 7 4 1]; % denominator coefficients of transfer function

% Convert transfer function to state-space model
[A, B, C, D] = tf2ss(num, den);
disp(A)
disp(B)
disp(C)
disp(D)

disp(eig(A))
sys_state_space = ss(A,B,C,D);

r = 5;
% Calculate the controllability and observability Gramians
Wc = gram(sys_state_space,'c');
Wo = gram(sys_state_space,'o');
T = chol(Wo,'lower');
Tinv = inv(T);

A_4 = Tinv*A*T;
B_4 = Tinv*B;
C_4 = C*T;
D_4 = D;

A_4 = A_4(1:r,1:r);
B_4 = B_4(1:r);
C_4 = C_4(1:r);
D_4 = D;
sys_state_space_4 = ss(A_4,B_4,C_4,D_4);
W_c = gram(sys_state_space,'c');
W_o = gram(sys_state_space,'o');

[U_c, S_c, V_c] = svd(W_c);
R = sqrtm(S_c)*U_c';
[U,S2,V_c] = svd(R*W_o*R');
mSigma = sqrtm(S2);
disp(mSigma);





% Define the system transfer function
s = tf('s');
G = (4*s^2-2*s-6)/(2*s^5+6*s^4+8*s^3+7*s^2+4*s+1);

% Convert to state-space representation
sys = ss(G);

% Balance the system
sys_bal = balancmr(sys);

% Perform balanced truncation with desired order
order = 3; % Choose the desired order
sys_bal_red = balred(sys_bal,order);

% Compute step response of the original and reduced systems
t = linspace(0,60,1000); % Time vector for simulation
[y,t] = step(sys,t); % Step response of original system
[y_red,t] = step(sys_bal_red,t); % Step response of reduced system
[y_bal,t] = step(sys_state_space, t);
[y_bal_4,t] = step(sys_state_space_4,t);

% Plot step response comparison
figure;
plot(t,y,'b-.',t,y_red,'r--');
legend('Original','Reduced','Balanced');
xlabel('Time (s)');
ylabel('Output');
title('Step Response Comparison');

figure;
plot(t,y_bal,'k-',t,y_bal_4,'b-');
legend('Original','Reduced','Balanced');
xlabel('Time (s)');
ylabel('Output');
title('Step Response Comparison');