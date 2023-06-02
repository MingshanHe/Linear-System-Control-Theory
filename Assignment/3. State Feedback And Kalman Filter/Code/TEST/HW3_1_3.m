A = [1 -2 1; 1  1 3; -1  4 0];
B = [1; -1; 0];
C = [1 0 2];
D = 0;
% Define State Space Equation System
sys = ss(A, B, C, D);
eigenvals = eig(A);
char_poly = poly(eigenvals);
% Compute State Feedback Gain K w/o 'place Funtion
C_ = [C' A'*C' A'*A'*C'];
C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
P = inv(C_*C_bar);
L_bar = [8 19 5];
L = (L_bar*P)'


% C = [B A*B A*A*B];
% C_bar = [1 -2 -8; 0 1 -2; 0 0 1];
% P = inv(C*C_bar);
% K_bar = [17 92 169];
% K = K_bar*P;
% Compute State Feedback Gain K w 'place Funtion
desired_poles = [-1 -2 -3];
% Compute the observer gain matrix using pole placement
L = place(A', C', desired_poles)'
inv(P)*L
eig(A-L*C)