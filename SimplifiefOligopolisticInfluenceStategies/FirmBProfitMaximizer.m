function [Solutions,Fvals,ExitFlags,Lambdas] = FirmBProfitMaximizer(Tolerance,MaxIterations,MaxFunctionEvaluations,Dimensionality,N,TA,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% The main functionality provided by this function is to obtain the set of
% solutions for the profit maximization problem of Firm B given the
% influence value TA which was selected by the Firm A.

% Dimensionality: a parameter that corresponds to the number of 
% optimization variables.
% N: is the number of different random initial poitns to be considered.
% That is, number N, corresponds to the number of nonlinear equations
% systems to be internally solved.
% lb: is the vector of lower bounds.
% ub: is the vector of upper bounds.
% Tolerance: is an internal optimization parameter that affects the
% convergence of the minimizer.
% MaxIterations: is an internal optimization parameter that controls the
% maximum number of iterations to be conducted by the optimizer.
% MaxFunctionEvaluations: is an internal optimization parameter that
% controls the maximum number of function evaluations to be conducted by
% the optimizer.

% Ensuring reproducibility of the results. 
rng default

% Set the lower and upper maximimization bounds.
lb = zeros(1,Dimensionality);
ub = ones(1,Dimensionality);

% Uniform sampling of the initial points.
initial_points = rand(N,Dimensionality);

% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,1);
ExitFlags = zeros(N,1);
% Preallocate cell array for storing the Lagrangian multipliers.
Lambdas = cell(N,1);

% Set optimization options and algorithm.
%options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',10000,'MaxFunctionEvaluations',100000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);
%options = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed');
options = optimoptions(@fmincon,'Algorithm','sqp','Display','off','ConstraintTolerance',Tolerance,'MaxIterations',MaxIterations,'MaxFunctionEvaluations',MaxFunctionEvaluations,'StepTolerance',Tolerance,'OptimalityTolerance',Tolerance);

% There exists a single non linear constraints which aims at ensuring the
% non-negative of the optimal profit for the second firm.
Cobj = @(TB)FirmBNonLinearConstraint(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Set the handle to the function that defines the objective function to be
% maximized by the second firm (Firm B). Take into consideration the fact
% that the maximization problm should be transformed into a minimization
% one by negative the objective function.
Fobj = @(TB)(-FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB));

% Set the linear equality-related constraints matrices Aeq and Beq.
Aeq = [];
Beq = [];

% Set the linear inequality-related constraints matrices A and B.
% Given that T = [TB] the set of linear constraints may be expressed by an
% inequality of the form: A * T <= B.
% -------------------------------------------------------------------------
% CASE I (beta >= 0):
% -------------------------------------------------------------------------
% For the case where the competition-related parameter beta is non-negative
% there exists a unique linear inequality constraint which according to the
% previously defined matrix form can be encoded by setting matrics A and B 
% such that:
% A = [1] and B = [1-TA]
% -------------------------------------------------------------------------
% CASE II (beta < 0):
% -------------------------------------------------------------------------
% For the case where the competition-related parametr beta is negative
% there exists an additional linear inequality constraint in order to
% ensure the positivity of the optimal demand function for Firm B which
% finally yields that:
%      -       -            -                                 -
%     |  1      |          |1-TA                               |
% A = |-alpha*LA|  and B = |beta*LB*TA+LA*LB*(alpha*PB+beta*PA)| 
%      -       -            -                                 -
% -------------------------------------------------------------------------

% Case I: (beta >=0)
if(beta>=0)
    A = 1;    % A is a [1 x 1] matrix.
    B = 1-TA; % B is a [1 x 1] matrix.
% Case II: (beta < 0)
else
    A = [1;-alpha*LA];                             % A is a [2 x 1] matrix.
    B = [1-TA;beta*LB*TA+LA*LB*(alpha*PB+beta*PA)];% B is a [2 x 1] matrix.
end

% Solve the constrained optimization problem for the first firm for each
% initial point.
parfor k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:),~,Lambdas{k},~,~] = fmincon(Fobj,initial_points(k,:),A,B,Aeq,Beq,lb,ub,Cobj,options);
end

end

