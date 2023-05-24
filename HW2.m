A = [5 3 5 4; -8 -6 0 -8; 3 3 -2 3; -1 -1 -5 0];
B = [1; -2; 0; 1];
C = [2 1 0 2];
disp(A*B)
disp(A*A*B)
disp(A*A*A*B)

C_ = [1 3 7 15; -2 -4 -8 -16; 0 0 0 0; 1 1 1 1];
disp(rank(C_))
C_2 = [1  0 0 0; 
      -2  2 0 0; 
       0  0 0 1;
       1  -2 1 0];
disp(rank(C_2))
disp(inv(C_2))

disp(inv(C_2)*A*C_2)
disp(inv(C_2)*B)

disp(C*A)
disp(C*A*A)
disp(C*A*A*A)
O = [2 1 0 2; 0 -2 0 0; 16 12 0 16; -32 -40 0 -32];
rank(O)
O_2 = [2    1 0 1; 
       0   -2 0 0; 
       16  12 1 0; 
      -32 -40 0 1];
rank(O_2)

% Define transfer function
num = [4 -2 -6];
den = [2  6  8 7 4 1];
sys_tf = tf(num, den);

% Find minimal realization
sys_minreal = minreal(sys_tf);

% Display transfer functions
disp('Original transfer function:')
tf(sys_tf)

disp('Minimal realization:')
tf(sys_minreal)
