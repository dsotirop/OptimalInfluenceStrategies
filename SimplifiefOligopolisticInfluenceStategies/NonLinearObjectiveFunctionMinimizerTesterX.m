% This script file provides the fundamental testing functionality for
% evaluating the process of solving the underlying constrained optimization 
% problem which takes into consideration the upper-limit-related Lagrangian
% multipliers.

% Clear workspace and command window
clc
clear all

% Initialize the external optimizations variables.
LA = 0.25;
LB = 0.25;
PA = 0.1;
PB = 0.1;
M = 0.3;
K = 0.3;
C = 0.0001;
G = 0.2; % Mind that G is the Gamma parameter defined elsewhere.

% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);

% Set the dimensionality of the search space.
Dimensionality = 4;
% Set the number of different initial points to be considered.
N = 300;
% Set lower and upper bounds for optimization variables.
lb = [0 0 0 0];
ub = [1 1 Inf Inf];
% Set the tolerance value for the minimizer. (Preferable value = 1e-10)
Tolerance = 1e-15;

% Run the solver.
[Solutions,Fvals,ExitFlags] = NonLinearObjectiveFunctionMinimizerX(Tolerance,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Run the First and Second Order Derivatives Checker.
[FD,SD] = OptimalityConditionsTester(Solutions,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Keep the solutions for which the corresponding exit flags are
% non-negative.
GoodSolutions = Solutions(ExitFlags>=0,:);

% Extract the representative solutions.
%[RepresentativeSolutions,DigitsAccuracy,~,Weight] = ExtractClusterRepresentativeSolutions(GoodSolutions);
