A = [0 1; -10 -20];
B = [0; 10];
C = [1 0];
D = 0;

Q = [1 0; 0 1];
R = 1;
K = lqr(A,B,Q,R);
disp(K);

sys = ss(A, B, C, D);

% Check stability
isStable = isstable(sys);

% Display the result
if isStable
    disp('The system is stable.');
else
    disp('The system is not stable.');
end