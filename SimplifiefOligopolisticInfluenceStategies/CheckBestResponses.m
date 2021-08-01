% This script file tests the computation procedure of Best Response
% Influence levels according to the individual optimization procedures:
% (i): FirmAProfitMaximizer and (ii): FirmBProfitMaximizer for a given
% setiing of external optimization parameters.

% Clear workspace and command window
clc
clear all

% Initialize the external optimizations variables.
LA = 0.25;
LB = 0.25;
PA = 0.1;
PB = 0.1;
M = 0.5;
K = 0.5;
C = 0.0001;
G = 0.2; % Mind that G is the Gamma parameter defined elsewhere.

% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);

% Set the dimensionality of the search space.
Dimensionality = 1;

% Set the number of different initial points to be considered.
N = 1000;

% Set the tolerance value for the minimizer. (Preferable value = 1e-10)
Tolerance = 1e-10;

% Set the value of TA for which TBopt will be computed.
TA = 0.1440;
% Set the value of TB for which TAopt will be computed.
TB = 0.1440;

% Run the solver for TAopt.
[Solutions_TA,Fvals_TA,ExitFlags_TA] = FirmAProfitMaximizer(Tolerance,Dimensionality,N,TB,C,G,LA,LB,PA,PB,alpha,beta,gamma);
% Run the solver for TBopt.
[Solutions_TB,Fvals_TB,ExitFlags_TB] = FirmBProfitMaximizer(Tolerance,Dimensionality,N,TA,C,G,LA,LB,PA,PB,alpha,beta,gamma);