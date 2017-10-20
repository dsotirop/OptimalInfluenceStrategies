%function [T1_opt,T2_opt,Fval,Flag] = GeneralOptimalInfluences(P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma)
function [T1_opt,T2_opt,Fval] = GeneralOptimalInfluences(P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma)

% This function provides the optimal influences T1_opt, T2_opt for a given
% set of the general optimal influence model parameters that are described  
% by the input variables.

% Set boundary values for the free parameters of the optimization
% problem.
T1_min = 0;
T1_max = 1-Theta2;
T2_min = 0;
T2_max = 1 - Theta1;



% The following set of lower and upper bounds should be enforced when
% the equality constraint T1 + T1 = 1 is active.
% T1_min = Theta1;
% T1_max = 1-Theta2;
% T2_min = Theta2;
% T2_max = 1 - Theta1;

% Set lower and upper bounds for the free parameters of the optimization.
lb = [T1_min,T2_min];
ub = [T1_max,T2_max];

% Run objective function minimizer.
N = 2;
Dimensionality = 2;
[Solution,Fval] = NonLinearGeneralObjectiveMinimizer(Dimensionality,N,lb,ub,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma);
T1_opt = Solution(1);
T2_opt = Solution(2);
Fval = Fval(1);

% Characterize optimal solutions when the equality constraint T1+T2=1 is
% active.
% if (T1_opt==Theta1)
%     Flag = 1;
% end;
% if(T1_opt==(1-Theta2))
%     Flag = 2;
% end;
% if(and((T1_opt>=Theta1),(T1_opt<=(1-Theta2))))
%     Flag = 0;
% end;

end
