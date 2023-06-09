% Define the transfer function
num = [1];
den = [1, 0.5, 0.1];
sys = tf(num, den);

% Check stability
isStable = isstable(sys);

if isStable
    disp('The transfer function is stable.');
else
    disp('The transfer function is unstable.');
end

clear;
% Define the polynomial coefficients
coeffs = [1, 2, 3, 4];

% Compute the roots
roots = roots(coeffs);

% Display the roots
disp('Roots:');
disp(roots);
