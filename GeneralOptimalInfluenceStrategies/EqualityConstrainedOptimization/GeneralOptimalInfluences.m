function [T1_opt,T2_opt,Fval,Flag1,Flag2] = GeneralOptimalInfluences(P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma,K)

% This function provides the optimal influences T1_opt, T2_opt for a given
% set of the general optimal influence model parameters that are described  
% by the input variables.

% Set boundary values for the free parameters of the optimization
% problem when the equality constraint T1 + T1 = K is active.
T1_min = max(0,Theta1+K-1);
T1_max = min(1-Theta2,K);
T2_min = max(0,Theta2+K-1);
T2_max = min(1-Theta1,K);

% Set boundary values for the free parameters of the optimization
% problem when the equality constraint is de-activated.
% T1_min = 0;
% T1_max = 1 - Theta2;
% T2_min = 0;
% T2_max = 1 - Theta1;

% Set lower and upper bounds for the free parameters of the optimization.
lb = [T1_min,T2_min];
ub = [T1_max,T2_max];

% Run objective function minimizer.
N = 2;
Dimensionality = 2;
[Solution,Fval] = NonLinearGeneralObjectiveMinimizer(Dimensionality,N,lb,ub,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma,K);
T1_opt = Solution(1);
T2_opt = Solution(2);
Fval = Fval(1);

% Characterize internal and external optimal solutions when the equality 
% constraint T1+T2=K is active.
% Flag1 indicates the status of the T1 optimal solution while Flag2
% indicates the status of the T2 optimal solution. 
% Thus, {Flag1,Flag2} in {-1,0,+1} where [-1] indicates a minimum external
% solutions, [0] indicates an internal solution and [+1] indicates a
% maximum external solution.

if(T1_opt==T1_min)
    Flag1 = -1;
else
    if(T1_opt==T1_max)
        Flag1 = +1;
    else
        Flag1 = 0;
    end;
end;

if(T2_opt==T2_min)
    Flag2 = -1;
else
    if(T2_opt==T2_max)
        Flag2 = +1;
    else
        Flag2 = 0;
    end;
end;


end
