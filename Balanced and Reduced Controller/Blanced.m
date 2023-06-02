% given state space model
num = [4 -2 -6];     % numerator coefficients of transfer function
den = [2 6 8 7 4 1]; % denominator coefficients of transfer function

% Convert transfer function to state-space model
[A, B, C, D] = tf2ss(num, den);
sys = ss(A, B, C, D);

% compute controllability and observability Gramians
Wc = gram(sys, 'c');
Wo = gram(sys, 'o');

% calculate transformation matrix T
T = chol(Wc,'lower') \ inv(sqrt(Wo));

% apply transformation to state space model
A_bal = T * A * inv(T);
B_bal = T * B;
C_bal = C * inv(T);

% choose reduced order for balanced realization
r =4;

% reduce order of balanced realization
Ab = A_bal(1:r, 1:r);
Bb = B_bal(1:r, :);
Cb = C_bal(:, 1:r);

% compute step response of original and reduced systems
t = 0:1:60;
y = step(sys, t);
y_bal = step(ss(Ab, Bb, Cb, D), t);

% plot step response of original and reduced systems
plot(t, y, t, y_bal);
legend('Original', 'Reduced');
