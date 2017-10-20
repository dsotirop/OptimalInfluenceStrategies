function [Solutions,Fvals,ExitFlags] = ConstrainedNonLinearEquationsSystemSolver(Dimensionality,N,lb,ub,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma)

% The main functionality provided by this function is to obtain a solution
% to the constrained system of nonlinear nonlinear equations. Handling the
% constrained system of nonlinear equations may be conducted by
% utlilizing the "fmincon" solver which internally reformulats the original
% problem as constrained optimization problem.

% Dimensionality: a parameter that corresponds to the number of 
% optimization variables.
% N: is the number of different random initial poitns to be considered.
% That is, number N, corresponds to the number of nonlinear equations
% systems to be internally solved.
% lb: is the vector of lower bounds.
% ub: is the vector of upper bounds.

% Ensuring reproducibility of the results. 
rng default

% Uniform sampling of the initial points.
initial_points = rand(N,Dimensionality);

% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,Dimensionality);
ExitFlags = zeros(N,1);

% Set optimization options and algorithm.
%options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter');
%options = optimoptions(@fmincon,'Algorithm','active-set','Display','iter');
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','ConstraintTolerance',1.0e-10,'MaxIterations',1000,'MaxFunctionEvaluations',1000,'StepTolerance',1.0e-10);


% Set the handle to the function that defines the system of nonlinear 
% equality constraints.
Cobj = @(T)NonLinearConstraints(T,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);

% Set the handle to the function that defines a dummy optimization problem.
Fobj = @(T)0;

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

% Solve the system of nonlinear equations for each quadraplet of initial 
% points.
for k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:)] = fmincon(Fobj,initial_points(k,:),A,B,Aeq,Beq,lb,ub,Cobj,options);
end;

end

