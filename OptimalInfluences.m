function [T1_opt,T2_opt,Fval] = OptimalInfluences(P1,P2,Theta,Delta,Gamma)

% This function provides the optimal influences T1_opt, T2_opt for a given
% set of the optimal influence model parameters that are described by the 
% input variables.

% Set boundary values for the free parameters of the optimization
% problem.
T1_min = 0;
T1_max = 1/2;
T2_min = 0;
T2_max = 1 - (Theta/2);

% Set lower and upper bounds for the free parameters of the optimization.
lb = [T1_min,T2_min];
ub = [T1_max,T2_max];

% Run objective function minimizer.
N = 20;
Dimensionality = 2;
[Solution,Fval] = NonLinearObjectiveMinimizer(Dimensionality,N,lb,ub,P1,P2,Theta,Delta,Gamma);
T1_opt = Solution(1);
T2_opt = Solution(2);
Fval = Fval(1);

end

