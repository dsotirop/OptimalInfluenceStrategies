function [T1_A_opt,T2_A_opt,T1_B_opt,T2_B_opt,Fval,OptFlag,DigitFlag] = OligopolisticOptimalInfluences(Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma)

% This function contols the nonlinear minimization process provided by
% the functions:
% (i):   NonLinearObjectiveFunctionMimimizer
% (ii):  SecondOrderConditionsTester
% (iii): LocalOptimaCharecterization

% The main functionality provided by this function is to determine the
% optimal quadraplet of influences (T1_A_opt,T2_A_opt,T1_B_opt,T2_B_opt)
% for a given set of the oligopolistic model paremeters that are passsed as
% input arguments.

% For the rest of the output arguments the following hold:
% [Fval]: is the value of the composite optimization metric which transforms
% the underlying system of nonlinear equations into a single objective
% nonlinear optimization problem with linear inequality constraints.
% Therfore, Fval, combines the nonlinear equalities that are assciated with 
% the first order conditions of the two optimization problems into a sigle 
% objective whose values are given by the function
% "NonLinearObjectiveFuncion". Ideal values for Fval are the ones that are
% infinitisimally close to zero.
% [OptFlag]: is an indicator of the number of "clusters" of different
% solutions that were find during the optimization process. The ideal value
% for OptFlag is 1 suggesting that all of the obtained solutions are
% actually members of the same "cluster" sharing a maximal number of same
% digits. When OptFlag > 1, the function returns a warning so that further
% investigation must take place.
% [DigitFlag]: is the number of common digits amongst the solutions
% pertaining to the "cluster" of the obtained solution.


% The following requirements are met:
% [1]: Only the subset of solutions for which ExitFlag = 1 should be taken 
%      into consideration for further processing. ExitFlags are returned by 
%      the "NonLinearObjectiveFunctionMinimizer" function where the value
%      of 1 indicates that the optimization process terminated without any
%      problems.
%
% [2]: The "clusters" of slightly dissimilar solutions should be identified.
%      Ideally, all valid solutions must pertain to the same "cluster".
%
% [3]: Each different cluster should be represented by the value which
%      corresponds to the maximal series of common digits amongst all the
%      solutions within a certain cluster.
%
% [4]: For each distinct representative solution second order optimality
%      tests should be made in order to ensure its validity. Valid
%      solutions are the ones for which the eigenvalues of the associated
%      Hessian matrix are all negative, suggesting the presence of a local
%      maximum for the profit functions F_A and F_B. In case, there exist
%      more than one "clusters" of valid solutions a warning should be
%      returned. If no valid solutions remain after checking the second
%      order optimality conditions then an error should be returned
%      enforcing the termination of the code execution.
%
% [5]: At the extreme case where more than one valid solutions are retained 
%      after ensuring that the second order conditions indicate the existence
%      of a local maximum then solution with the higher weigh value is to
%      be finally returned. (Of course, this is a rather arbitrary decision 
%      to be made and in case this is a very frequent event additional 
%      measures should be taken).

% Set the dimensionality of the search space.
Dimensionality = 4;

% Set lower and upper bounds for optimization variables.
lb = zeros(1,Dimensionality);
ub = ones(1,Dimensionality);

% Set the number of different initial points to be considered.
N = 100;
% Set the tolerance value for the minimizer. (Secure value = 1e-10)
Tolerance = 1e-10;
% Set the maximum iterations value for the minimizer. (Secure value = 1000)
MaxIterations = 1000;
% Set the maximum function evaluations value for the minimizer. (Secure value = 10000)
MaxFunctionEvaluations = 10000;

% Run the solver.
[Solutions,~,ExitFlags] = NonLinearObjectiveFunctionMinimizer(Tolerance,MaxIterations,MaxFunctionEvaluations,Dimensionality,N,lb,ub,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);

% Deterine solutions for which ExitFlag = 1.
% GoodSolutions = Solutions(ExitFlags==1,:);
% if(isempty(GoodSolutions))
%     GoodSolutions = Solutions;
% end;
GoodSolutions = Solutions;

% Identify the representative solutions.
[RepresentativeSolutions,DigitsAccuracy,~,Weight] = ExtractClusterRepresentativeSolutions(GoodSolutions);
% Check Second Order Optimality Conditions.
% Mind that if the pair of obtained Hessian eigenvalues with respect to both 
% F_A and F_B are negative the corresponding value of the vector HF_Flags 
% will be -2. This is indicative of the fact that both optimization
% function exhibit a local maximum at the given solution quadraplet. 
[HF_A_Eigs,HF_B_Eigs] = SecondOrderConditionsTester(RepresentativeSolutions,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);
[~,~,HF_Flags] = LocalOptimalCharacterization(HF_A_Eigs,HF_B_Eigs);
% Get the subset of valid solutions according to the second order conditions
% examination. Also, obtain the corresponding vectors for the variables 
% DigitsAccuracy and Weight.
ValidSolutions = RepresentativeSolutions(HF_Flags == -2,:);
DigitsAccuracy = DigitsAccuracy(HF_Flags == -2);
Weight = Weight(HF_Flags == -2);

% Get the number of final valid solutions.
valid_solutions_number = length(Weight);

% In case the number of valid solutions is zero, then an error should be
% returned which will enforce an unexpected ternitation of the code
% execution.
if(valid_solutions_number == 0)
    error('No valid solutions where found for the current setting of external parameters');
else
    % In case more than one valid solutions where obtained, then throw a 
    % warning without terminating the execution process and keep the
    % solution which is associated with the higher weight value (that is,
    % the most frequent one).
    if(valid_solutions_number > 1)
        warning('More than one valid solutions were found.');
        % Get the maximum weight index.
        [~,max_weight_index] = max(Weight);
        % Obtain the corresponding entry from matrix ValidSolutions.
        valid_solution = ValidSolutions(max_weight_index,:);
        % Obtain the corresponding entry from vector DigitsAccuracy.
        digits_accuracy = DigitsAccuracy(max_weight_index);
    else
        % Only one valid solutions was retained.
        valid_solution = ValidSolutions;
        digits_accuracy = DigitsAccuracy;
    end;
end;

% Get the individual solutions for the optimal influence strategies.
T1_A_opt = valid_solution(1);
T2_A_opt = valid_solution(2);
T1_B_opt = valid_solution(3);
T2_B_opt = valid_solution(4);
% Set the optimal influence quadraplet.
T_opt = [T1_A_opt,T2_A_opt,T1_B_opt,T2_B_opt];

% Set the OptFlag output argument.
OptFlag = valid_solutions_number;

% Set the DigitFlag output argument.
DigitFlag = digits_accuracy;

% Set the Fval output argument.
Fval = NonLinearObjectiveFunction(T_opt,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);




end

