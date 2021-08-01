% This script file provides the fundamental testing functionality for
% evaluating the process of solving the underlying constrained optimization 
% problem. 

% Clear workspace and command window
clc
clear all

% Initialize the external optimizations variables.
LA = 0.50;
LB = 0.50;
PA = 0.02;
PB = 0.50;
M = 0.5;
K = 0.5;
C = 0.00;
G = 0.2; % Mind that G is the Gamma parameter defined elsewhere.

% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);

% Set the dimensionality of the search space.
Dimensionality = 2;
% Set the number of different in3itial points to be considered.
N = 1000;
% Set lower and upper bounds for optimization variables.
lb = zeros(1,Dimensionality);
ub = ones(1,Dimensionality);
% Set the tolerance value for the minimizer. (Preferable value = 1e-10)
Tolerance = 1e-15;

% Run the solver.
[Solutions,Fvals,ExitFlags] = NonLinearObjectiveFunctionMinimizer(Tolerance,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Run the First and Second Order Derivatives Checker.
[FD,SD] = OptimalityConditionsTester(Solutions,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Keep the solutions for which the corresponding exit flags are
% non-negative.
GoodSolutions = Solutions(ExitFlags>=0,:);

% Get the mean point.
mean_point = mean(GoodSolutions);
% Plot GoodSolutions.
figure('Name','Optimal Solutions');
hold on
plot(GoodSolutions(:,1),GoodSolutions(:,2),'.r','LineWidth',1.4);
plot(mean_point(1),mean_point(2),'*k','LineWidth',1.4);
xlabel('T_A');
ylabel('T_B');
axis([0 1 0 1]);
grid on
% Extract the representative solutions.
%[RepresentativeSolutions,DigitsAccuracy,~,Weight] = ExtractClusterRepresentativeSolutions(GoodSolutions);
