function [Solutions,Fvals,ExitFlags] = NonLinearObjectiveFunctionMinimizer(Tolerance,MaxIterations,MaxFunctionEvaluations,Dimensionality,N,lb,ub,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma)

% The main functionality provided by this function is to obtain a solution
% to the bounded and constrained system of nonlinear equations that underly
% the problem of determining the optimal influence strategies within an
% oligopolistic enviroment. The system of nonlinear equations, however,
% will be reformulated as a bounded / constrained optimization problem as
% described in the function m-file "NonLinearObjectiveFunction". Therefore,
% optimization will be performed by utilizing the "fmincon" solver.

% Dimensionality: a parameter that corresponds to the number of 
% optimization variables.
% N: is the number of different random initial poitns to be considered.
% That is, number N, corresponds to the number of nonlinear equations
% systems to be internally solved.
% lb: is the vector of lower bounds.
% ub: is the vector of upper bounds.
% Tolerance: is an internal optimization parameter that affects the
% convergence of the minimizer as well as the "similarity" of the obtained
% solutions. Therefore, lower values for the "Tolerance" parameter are
% preferable which are, however, associated with a slower rate of
% convergence. For larger values of the "Tolerance" parameter the obtained
% solutions may be falsly regarded as "different". 

% Ensuring reproducibility of the results. 
rng default

% Uniform sampling of the initial points.
initial_points = rand(N,Dimensionality);

% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,1);
ExitFlags = zeros(N,1);

% More accurate optimization results may be obtained by setting each one of
% the tolerance parameters equal to 1.0e-10. However, faster convergence
% for a large number of initial random solutions demands for a much higher
% tolerance value.

tolerance_value = Tolerance;
max_iterations_value = MaxIterations;
max_function_evaluations_value = MaxFunctionEvaluations;

% Set optimization options and algorithm.
%options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter');
%options = optimoptions(@fmincon,'Algorithm','active-set','Display','iter');
options = optimoptions(@fmincon,'Algorithm','sqp','Display','off','ConstraintTolerance',tolerance_value,'MaxIterations',max_iterations_value,'MaxFunctionEvaluations',max_function_evaluations_value,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);

% Mind that there exist no nonlinear equality constraints for this
% particular optimization problem.
Cobj = [];

% Set the handle to the function that defines the reformulated optimization problem.
Fobj = @(T)NonLinearObjectiveFunction(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);
%Fobj = @(T)NonLinearObjectiveFunctionX(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);
%Fobj = @(T)NonLinearObjectiveFunctionXX(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);

% Set the linear equality-related constraints matrices Aeq and Beq.
Aeq = [];
Beq = [];

% Set the linear inequality-related constraints matrices A and B.
% Mind that there exist two linear inequality constraints on the
% optimization variables T1_A, T2_A, T1_B and T2_B with respect to the
% external variables Theta1 and Theta2:
%                   T1_A + T1_B <= 1 - Theta2 [I]
%                   T2_A + T2_B <= 1 - Theta1 [II]
%
% Given, that T = [T1_A;T2_A;T1_B;T2_B], matrices A and B will be of the
% following form:
%                |1 0 1 0|         |1 - Theta2|
%            A = |0 1 0 1| and B = |1 - Theta1| 

A = [1 0 1 0;0 1 0 1];   % A should be a [2 x 4] matrix.
B = [1-Theta2;1-Theta1]; % B should be a [2 x 1] vector.

% A = [];
% B = [];

% Solve the system of nonlinear equations for each quadraplet of initial 
% points.
for k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:)] = fmincon(Fobj,initial_points(k,:),A,B,Aeq,Beq,lb,ub,Cobj,options);
end;

end

