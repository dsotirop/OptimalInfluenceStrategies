 function [Solutions,Fvals,ExitFlags,Cineq,Ceq] = NonLinearObjectiveFunctionMinimizer(DisplayFlag,Tolerance,MaxIterations,MaxFunctionEvaluations,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma)

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
% DisplayFlag: is a string taking values within the discrete set 
%              {'on','iter'} indicating whether the individual optimization
%              steps will be displayed.
% Tolerance: is an internal optimization parameter that affects the
% convergence of the minimizer.
% MaxIterations: is an internal optimization parameter that controls the
% maximum number of iterations to be conducted by the optimizer.
% MaxFunctionEvaluations: is an internal optimization parameter that
% controls the maximum number of function evaluations to be conducted by
% the optimizer.

% Ensuring reproducibility of the results. 
rng default

% Uniform sampling of the initial points.
initial_points = rand(N,Dimensionality);

% Preallocate matrices of solutions, fvals and exit_flags.
Solutions = zeros(N,Dimensionality);
Fvals = zeros(N,1);
ExitFlags = zeros(N,1);
% Preallocate matrices storing the evaluation of the inequality and
% equality constraints at the solution points. Mind that the second
% argument is actually determined in function NonLinearConstraintFunction.m
Cineq = zeros(N,2);
% Case I (beta >= 0):
if(beta>=0)
    Ceq = zeros(N,1);
% Case II: (beta < 0):
else
    Ceq = zeros(N,3);
end

tolerance_value = Tolerance;

% Set optimization options and algorithm.
%options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',250,'MaxFunctionEvaluations',10000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);
%options = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed','ConstraintTolerance',tolerance_value,'MaxIterations',250,'MaxFunctionEvaluations',10000,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);
options = optimoptions(@fmincon,'Algorithm','sqp','Display',DisplayFlag,'FunctionTolerance',tolerance_value,'ConstraintTolerance',tolerance_value,'MaxIterations',MaxIterations,'MaxFunctionEvaluations',MaxFunctionEvaluations,'StepTolerance',tolerance_value,'OptimalityTolerance',tolerance_value);

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
% Mind that the number of linear inequality constraints depends on the sign
% of the parameter beta.
% Case I:  beta > 0 ==> 1 Linear Inequality Constraint ([I])
% Case II: beta < 0 ==> 3 Linear Inequality Constraints ([I],[II],[III])
% The fundamental linear inequality constraint on the pair of optimization
% variables TA and TB which emerges from the stochasticity of the social
% interaction matrix T:
%                   TA + TB + 0 * MUa + 0 * MUb + 0 * Muab <= 1 [I]
%
% There exist two additional linear inequality constraints which are
% associated with enforcing the non-negativity of the optimal demand
% functions (Qa >= 0 and Qb >= 0) for both firms when the competition  
% related parameter beta is negative (beta < 0):
%
% - alpha * LB * TA - beta * LA * TB <= LA * LB * (alpha * PA + beta * PB) [II]
% - beta * LB * TA - alpha * LA * TB <= LA * LB * (alpha * PB + beta * PA) [III] 

% Case I (beta >= 0):
% Given that T = [TA;TB;MUa;MUb;MUab], matrices A and B will be of the 
% following form:
%                A = [1 1 0 0 0] and B = [1].
% Case II (beta <0):
% Given that T = [TA;TB;MUa;MUb;MUab1;MUab2;MUab3]
%           |1             1      0  0  0 0 0|     |1                       |
% A =       |-alpha*LB -beta*LA   0  0  0 0 0| B = |LA*LB*(alpha*PA+beta*PB)|
%           |-beta*LB  -alpha*LB  0  0  0 0 0|     |LA*LB*(alpha*PB+beta*PA)|

if(beta >= 0)
    A = [1 1 0]; % A should be a [1 x 3] matrix.
    B = 1; % B should be a [1 x 1] matrix.
else
    A = [1 1 0 0 0;-alpha*LB -beta*LA 0  0  0;-beta*LB -alpha*LA 0 0 0]; % A should be a [3 x 5] matrix.
    B = [1;LA*LB*(alpha*PA+beta*PB);LA*LB*(alpha*PB+beta*PA)]; % B should be a [3 x 1] matrix.
end
% A = [];
% B = [];

% Solve the combined minimization problem for each pair of initial points.
parfor k = 1:1:N
    [Solutions(k,:),Fvals(k,:),ExitFlags(k,:)] = fmincon(Fobj,initial_points(k,:),A,B,Aeq,Beq,lb,ub,Cobj,options);
    % Extract the solution point Topt = [TAopt TBopt].
    Topt = Solutions(k,:);
    % Evaluate the non-linear contraint function at the extracted solution
    % point.
    [Cineq_opt,Ceq_opt] = NonLinearConstraintFunction(Topt,C,G,LA,LB,PA,PB,alpha,beta,gamma);
    % Populate matrices Cineq and Ceq.
    if(~isempty(Cineq_opt))
        Cineq(k,:) = Cineq_opt';
    end
    Ceq(k,:) = Ceq_opt';
end

end