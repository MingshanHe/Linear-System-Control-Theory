
%Method1:
num = [4 -2 -6];     % numerator coefficients of transfer function
den = [2 6 8 7 4 1]; % denominator coefficients of transfer function
% Convert transfer function to state-space model
[A, B, C, D] = tf2ss(num, den);
sys = ss(A, B, C, D);
% Compute the controllability and observability gramians
Wc = gram(sys, 'c');
Wo = gram(sys, 'o');
% Compute the SVD of the gramians
[Uc,Sc,Vc] = svd(Wc);
R = sqrtm(Sc)*Uc';
[Uo,So,Vo] = svd(R*Wo*R');
mSigma = sqrtm(So);
% Determine the balancing transformation matrix
P = mSigma*Uo'*inv(R');
% Compute the balanced realization
Ab = P * A * inv(P);
Bb = P * B;
Cb = C * inv(P);
% Define the desired reduced order
% Truncate the balanced realization to obtain the reduced-order system\
A_red_3 = Ab(1:3,1:3);
B_red_3 = Bb(1:3,:);
C_red_3 = Cb(:,1:3);

A_red_2 = Ab(1:2,1:2);
B_red_2 = Bb(1:2,:);
C_red_2 = Cb(:,1:2);
% Convert the reduced-order system back to standard state-space form
sys_red_3 = ss(A_red_3,B_red_3,C_red_3,D);
sys_red_2 = ss(A_red_2,B_red_2,C_red_2,D);
% Compare the step response of the original and reduced systems
t = 0:0.1:60;
u = ones(size(t));
[y,t] = step(sys,t);
[y_red_3,~] = step(sys_red_3,t);
[y_red_2,~] = step(sys_red_2,t);
figure;
plot(t,y,'k',t,y_red_3,'r--',t,y_red_2,'b-');
legend('Original system','Reduced system 3','Reduced system 2');
xlabel('Time (s)');
ylabel('Output');
title('Method1');
%Method2:
% Define the system transfer function
s = tf('s');
G = (4*s^2-2*s-6)/(2*s^5+6*s^4+8*s^3+7*s^2+4*s+1);
% Convert to state-space representation
sys = ss(G);
% Balance the system
sys_bal_3 = balancmr(sys,3);
% Perform balanced truncation with desired order
order = 3; % Choose the desired order
sys_red_3 = balred(sys_bal_3,order);
order = 2; % Choose the desired order
sys_bal_2 = balancmr(sys,2);
sys_red_2 = balred(sys_bal_2,order);
% Compute step response of the original and reduced systems
t = linspace(0,60,1000); % Time vector for simulation
[y,t] = step(sys,t); % Step response of original system
[y_red_3,t] = step(sys_red_3, t);
[y_red_2,t] = step(sys_red_2,t);

% Plot step response comparison
figure;
plot(t,y,'k',t,y_red_3,'r--',t,y_red_2,'b-');
legend('Original system','Reduced system 3','Reduced system 2');
xlabel('Time (s)');
ylabel('Output');
title('Method2');
