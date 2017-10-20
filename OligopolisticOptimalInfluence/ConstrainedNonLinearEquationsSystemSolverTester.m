% This is a script file the testing the process of solving the underlying
% bounded system of nonlinear equations that is considered in order 
% to obtain the optimal influence strategies wihtin an oligopolistic 
% environment.

% Clear workspace and command window
clc
clear all

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
M = 0.25;
K = 0.25;
C = 0.25;
Gamma = 10;

% Set the dimensionality of the search space.
Dimensionality = 4;
% Set the number of different initial points to be considered.
N = 1000;
% Set lower and upper bounds for optimization variables.
lb = zeros(1,Dimensionality);
ub = ones(1,Dimensionality);
% Run the solver.
[Solutions,Fvals,ExitFlags] = ConstrainedNonLinearEquationsSystemSolver(Dimensionality,N,lb,ub,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);

% Get obtained optimal solutions per variable.
T1_A_opt = Solutions(:,1);
T2_A_opt = Solutions(:,2);
T1_B_opt = Solutions(:,3);
T2_B_opt = Solutions(:,4);

% Plot obtained optimal solutions optimal
figure('Name','Optimal Solutions Points');
subplot(2,1,1);
plot(T1_A_opt,T2_A_opt,'*r','LineWidth',1.4);
xlabel('T_1');
ylabel('T_2');
legend({'Firm A'});
grid on
subplot(2,1,2);
plot(T1_B_opt,T2_B_opt,'*g','LineWidth',1.4);
xlabel('T_1');
ylabel('T_2');
legend({'Firm B'});
grid on