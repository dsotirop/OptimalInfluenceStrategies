function [TAopt,TBopt,FilterFlag,DigitsAccuracy] = SimplifiedOligopolisticOptimalInfluences(N,DisplayFlag,Tolerance,FvalsTolerance,DerivativeTolerance,MaxIterations,MaxFunctionEvaluations,MinimumDigitsAccuracy,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime)

% This function controls the nonlinear minimization process implemented by
% the following functions:
% (i):   NonLinearObjectiveFunctionMinimizer.
% (ii):  OptimalityConditionsTester.
% (iii): RetrieveOptimalModelParameters.
% (iv):  FilterSolutions.
% (vi):  ExtractClusterRepresentativeSolutions.

% Besides the external model parameters that need to be specified there
% exist some extra regulatory input arguments which are required, such as:
% (i):   N (number of random initial points) (Preferable value = 200).
% (ii):  Tolerance (internal optimization stopping threshold) (Preferable value = 1e-15).
% (iii): FvalsTolerance (specifies acceptable values for the overall
%                        minimimization objective) (Preferable value = 1e-15).
% (iv):  DerivativeTolerance (threshold value indicating that first order 
%                             or second order derivatives are actually zero)
%                            (Preferable value = 1e-10).
% (v):   MaxIterations (internal optimization stopping criterion)
%                      (Preferable value = 1000).
% (vi):  MaxFunctionEvaluations (internal optimization stopping criterion) 
%                               (Preferable value = 10000).
% (vii): DisplayFlag (controls whether internal optimization output will be
%                     displayed) (Preferable value = 'off').

% FilterFlag: characterizes the status of the obtained solutions. 
%             (0) indicates that the optimization process resulted into
%             valid solution points. Negative values indicate either that
%             the optimization process failed within the optimizer (-1) or
%             the filtering process eliminated the complete set of the
%             obtained solutions {(-2),(-3),(-4),(-5),(-6)}.

% IMPORTANT NOTE!!!
% Mind that the underlying continuous game may not have an equlibrium point
% for any given configuration of the external parameters which is
% indicated by a negative FilterFlag value. In such a case, the rest of the 
% output arguments {TAopt,TBopt,DigitsAccuracy} will be assigned the value
% of(-1).

% Set the number of firms.
FirmsNumber = 2;

% Set the dimensionality of the search space for the case beta >= 0:
if(beta>=0)
    Dimensionality = 3;
else
% Set the dimensionality of the search space for the case beta < 0:
    Dimensionality = 5;
end

% Set lower bounds for the optimization variables:
lb = zeros(1,Dimensionality);
% Set upper bounds for the optimization variables for the case beta < 0:
if(beta>=0)
    ub = [1 1 Inf];
% Set upper bounds for the optimization variables for the case beta > 0:
else
    ub = [1 1 Inf Inf Inf];
end

% Run the solver.
[Solutions,Fvals,ExitFlags,Cineq,Ceq] = NonLinearObjectiveFunctionMinimizer(DisplayFlag,Tolerance,MaxIterations,MaxFunctionEvaluations,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Filter Solutions.
[~,~,Solutions,Fvals,FilterFlag,~,~,~,~,~,~,~,~,~,~] = FilterSolutions(Solutions,ExitFlags,Fvals,Ceq,Cineq,FvalsTolerance,DerivativeTolerance,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);

% Obtain optimal investment levels.
% Check whether the underlying continuous game has a NE.
if(FilterFlag == 0)
    % Extract Representative Solutions.
    [RSOL,DigitsAccuracy,~] = ExtractClusterRepresentativeSolutions(Solutions,Fvals,MinimumDigitsAccuracy,FirmsNumber);
    % Separate Optimal Solutions.
    TAopt = RSOL(1);
    TBopt = RSOL(2);
else
    % The underlying continuous game does not possess a NE.
    warning('Optimization procedure failed: FilterFlag %d',FilterFlag);
    TAopt = -1;
    TBopt = -1;
    DigitsAccuracy = -1;
end
    

end

