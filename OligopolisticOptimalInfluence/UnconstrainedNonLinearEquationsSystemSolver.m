function [Solutions,Fvals,ExitFlags] = UnconstrainedNonLinearEquationsSystemSolver(Dimensionality,N,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma)

% The main functionality provided by this function is to obtain a solution
% to the unconstraned system of nonlinear nonlinear equations. Handling the
% unconstrained system of nonlinear equations may be conducted by
% utlilizing the "fsolve" solver.
                                                                    
% Dimensionality: a parameter that corresponds to the number of 
% optimization variables.
% N: is the number of different random initial poitns to be considered.
% That is, number N, corresponds to the number of nonlinear equations
% systems to be internally solved.

% Ensuring reproducibility of the results. 
rng default

% Uniform sampling of the initial points.
initial_points = 100*rand(N,Dimensionality);
% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,Dimensionality);
ExitFlags = zeros(N,1);
% Set optimization options and algorithm.
%options = optimoptions('fsolve','Display','iter','Algorithm','trust-region-dogleg');
%options = optimoptions('fsolve','Display','iter','Algorithm','trust-region-reflective');
options = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt','InitDamping',0.005);

% Set the handle to the function that defines the system of nonlinear 
% equations. 
Fobj = @(T)NonLinearEquationsSystem(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);

% Solve the system of nonlinear equations for each quadraplet of initial 
% points.

for k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:)] = fsolve(Fobj,initial_points(k,:),options);
end;

end

