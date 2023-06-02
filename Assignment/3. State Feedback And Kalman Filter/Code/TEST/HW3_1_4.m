% Define system matrices in controllable canonical form
A_cont = [1 1; -1 -1];
B_cont = [0; 1];
C_cont = [1 0];
D_cont = 0;

% Compute observability matrix
O = obsv(A_cont, C_cont);

% Check observability
if rank(O) < size(A_cont, 1)
    disp('System is not completely observable.');
    return;
end

% Compute transformation matrix
T = inv(O);

% Transform system matrices to observable canonical form
A_obs = T * A_cont * inv(T);
B_obs = T * B_cont;
C_obs = C_cont * inv(T);
D_obs = D_cont;

% Display the transformed system matrices
disp('Observable Canonical Form Matrices:');
disp('A_obs:');
disp(A_obs);
disp('B_obs:');
disp(B_obs);
disp('C_obs:');
disp(C_obs);
disp('D_obs:');
disp(D_obs);
