% This is a script file the testing the process of solving the underlying
% bounded / constrained system of nonlinear equations  
% that is considered in order to obtain the optimal influence strategies  
% wihtin an oligopolistic environment. Mind that solving the system of
% nonlinear equations is internally handled as a bounded / constrained
% optimization problem.

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
P_A_1 = 0.7;
P_A_2 = 0.7;
P_B_1 = 0.7;
P_B_2 = 0.7;
M = 0.3; % Good values are 0.2 0.3 0.4
K = 0.06; % Good values are 0.2 0.3 0.4
C = 0.0001;
Gamma = 0.5;

% Set the dimensionality of the search space.
Dimensionality = 4;
% Set the number of different initial points to be considered.
N = 500;
% Set lower and upper bounds for optimization variables.
lb = zeros(1,Dimensionality);
ub = ones(1,Dimensionality);
% Set the tolerance value for the minimizer. (Preferable value = 1e-10)
Tolerance = 1e-10;

% Run the solver.
[Solutions,Fvals,ExitFlags] = NonLinearObjectiveFunctionMinimizer(Tolerance,Dimensionality,N,lb,ub,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);

fprintf('Checking second order optimality conditions...\n');
% Check Second Order Optimality Conditions.
[CHF_A,CHF_B,HF_A_Eigs,HF_B_Eigs] = SecondOrderConditionsTester(Solutions,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);
[HF_A_Flags,HF_B_Flags,HF_Flags] = LocalOptimalCharacterization(HF_A_Eigs,HF_B_Eigs);
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

% Plot obtained minimum values.
figure('Name','Optimal F Values');
plot([1:N],Fvals,'*b','LineWidth',1.4);
xlabel('Solution ID');
ylabel('Fobj');
grid on

% Plot the eigen values for both Hessian matrices.
figure('Name','F_A Hessian Matrix Eigenvalues');
plot(HF_A_Eigs(:,1),HF_A_Eigs(:,2),'*r','LineWidth',1.4);
xlabel('L_{A,1}');
ylabel('L_{A,2}');
grid on

figure('Name','F_B Hessian Matrix Eigenvalues');
plot(HF_B_Eigs(:,1),HF_B_Eigs(:,2),'*g','LineWidth',1.4);
xlabel('L_{B,1}');
ylabel('L_{B,2}');
grid on
