function [Solutions,Fvals,ExitFlags] = NonLinearObjectiveFunctionMinimizer(Tolerance,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% The main functionality provided by this function is to obtain a solution
% for the combined minimization problem which arises by merging the KKT
% conditions associated with the pair of constrained optimization problems  
% that underlie the two-player continuous game which governs the proposed  
% oligopolistic framework between the two firms A and B which attempt to  
% optimally influence the single consumer C. 

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
% Filter the randomly generated initial solution points so that they fall 
% within the TA + TB <= 1 region.
% initial_points = initial_points(sum(initial_points,2)<=1,:);
% % Re-estimate the number of points.
% N = size(initial_points,1);

% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,1);
ExitFlags = zeros(N,1);

tolerance_value = Tolerance;

% Set optimization options and algorithm.
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',250,'MaxFunctionEvaluations',10000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);
%options = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',250,'MaxFunctionEvaluations',10000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);
%options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',500,'MaxFunctionEvaluations',100000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);

% The following line of code should be commented out in case no non-linear
% inequality constraints should be imposed on the combined optimization
% problem.
Cobj = @(T)NonLinearConstraintFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Set the handle to the function that defines the reformulated optimization problem.
Fobj = @(T)NonLinearObjectiveFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Set the linear equality-related constraints matrices Aeq and Beq.
Aeq = [];
Beq = [];

% Set the linear inequality-related constraints matrices A and B.
% Mind that there exists a single linear inequality constraint on the pair
% of optimization variables TA and TB which emerges from the stochasticity
% of the social interaction matrix T:
%                   TA + TB <=1 [I]
%
% Given that T = [TA;TB], matrices A and B will be of the following form:
%
% A = [1 1] and B = [1].

A = [1 1]; % A should be a [1 x 2] matrix.
B = [1]; % B should be a [1 x 1] matrix.

% Solve the combined minimization problem for each pair of initial points.
for k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:)] = fmincon(Fobj,initial_points(k,:),A,B,Aeq,Beq,lb,ub,Cobj,options);
end;


end

