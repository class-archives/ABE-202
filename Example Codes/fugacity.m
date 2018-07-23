% Enter Tc, Pc, w into EOS
Tc = 190.6; % K
Pc = 4.641; % MPa
w = 0.0115;

% There's a smart way to guess P -- calculate the zeroes of the derivative
% of the curve

peng(Tc,Pc,w);

% Enter T, initial guess P
T = input('Input a value of T: ');
P = input('Input an initial guess for P: ');

% Compute EOS specific parameters (a, b, etc.)

% Compute A, B; 
A = (a * P) / (R^2 * T^2);
B = (P * b) / (R * T);
              
% Solve EOS, cubic equation for Zv, Zl

% Compute fV

% Compute fL

% Test fugacity [0.0001 > abs(fL/fV - 1)]
% This is testing error tolerance; it can be changed if needed

% If yes, equilibrium is achieved (output results)

% If no, choose Pnew = Pold * (fL/fV)
% Re-compute A, B