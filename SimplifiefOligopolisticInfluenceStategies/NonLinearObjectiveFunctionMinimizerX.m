function [Solutions,Fvals,ExitFlags] = NonLinearObjectiveFunctionMinimizerX(Tolerance,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function is an extended version of
% "NonLinearObjectiveFunctionMinimizer". Its main purpose is to provide a
% solution to the combined minimization problem which arises by merging the KKT
% conditions associated with the pair of constrained optimization problems  
% that underlie the two-player continuous game which governs the proposed  
% oligopolistic framework between the two firms A and B which attempt to  
% optimally influence the single consumer C.

% The main differentation of this function is that it takes into
% consideration the Lagrange-related optimization variables bA and bB 
% which are associated with the upper bounds of the influence-related
% optimization variables.

% Dimensionality: a parameter that corresponds to the number of 
% optimization variables.
% N: is the number of different random initial poitns to be considered.
% That is, number N, corresponds to the number of nonlinear equations
% systems to be internally solved.
% lb: is the vector of lower bounds.
% ub: is the vector of upper bounds.
% Tolerance: is an internal optimization parameter that affects the
% convergence of the minimizer. 

% Ensuring reproducibility of the results. 
rng default

% Uniform sampling of the initial points.
initial_points = rand(N,Dimensionality);

% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,1);
ExitFlags = zeros(N,1);

tolerance_value = Tolerance;

% Set optimization options and algorithm.
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',250,'MaxFunctionEvaluations',100000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);
%options = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed');
%options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',100000,'MaxFunctionEvaluations',100000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);


% The following line of code should be commented out in case no non-linear
% inequality constraints should be imposed on the combined optimization
% problem.
Cobj = @(T)NonLinearConstraintFunctionX(T,C,G,LA,LB,PA,PB,alpha,beta,gamma);
%Cobj = [];

% Set the handle to the function that defines the reformulated optimization problem.
Fobj = @(T)NonLinearObjectiveFunctionX(T,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Set the linear equality-related constraints matrices Aeq and Beq.
Aeq = [];
Beq = [];

% Set the linear inequality-related constraints matrices A and B.
% Mind that there exists a single linear inequality constraint on the pair
% of optimization variables TA and TB which emerges from the stochasticity
% of the social interaction matrix T:
%                   TA + TB + 0 * bA + 0 * bB <=1 [I]
%
% Given that T = [TA;TB;bA;bB], matrices A and B will be of the following form:
%
% A = [1 1 0 0] and B = [1].

A = [1 1 0 0]; % A should be a [1 x 4] matrix.
B = [1]; % B should be a [1 x 1] matrix.

% Solve the combined minimization problem for each pair of initial points.
for k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:)] = fmincon(Fobj,initial_points(k,:),A,B,Aeq,Beq,lb,ub,Cobj,options);
end;

end

