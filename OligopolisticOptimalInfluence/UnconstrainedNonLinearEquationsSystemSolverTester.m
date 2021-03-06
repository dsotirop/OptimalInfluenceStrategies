% This is a script file the testing the process of solving the underlying
% unconstrained system of nonlinear equations that is considered in order 
% to obtain the optimal influence strategies wihtin an oligopolistic 
% environment.

% Initialize the external optimizations variables.
Lambda_A_1 = 0.25;
Lambda_A_2 = 0.25;
Lambda_B_1 = 0.25;
Lambda_B_2 = 0.25;
Theta1 = 0.2;
Theta2 = 0.2;
P_A_1 = 0.1;
P_A_2 = 0.1;
P_B_1 = 0.1;
P_B_2 = 0.1;
M = 0.01;
K = 0.01;
C = 0.01;
Gamma = 0.01;

% Set the dimensionality of the search space.
Dimensionality = 4;
% Set the number of different initial points to be considered.
N = 300;
% Run the solver.
[Solutions,Fvals,ExitFlags] = UnconstrainedNonLinearEquationsSystemSolver(Dimensionality,N,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);