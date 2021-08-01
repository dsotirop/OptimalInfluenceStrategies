function [HF_A_Eigs,HF_B_Eigs] = SecondOrderConditionsTester(Solutions,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma)

% This code function checks the second order optimality conditions for the
% objective functions F_A and F_B to be maximized for the given quadraplets
% of solutions. Second Order Optimality tests are conducted by computing
% the Hessian matrix for each configuration of input parameters.

% Mind that for a local maximum to be identified the corresponding Hessian
% matrix needs to be Negative Definite. Negative Definiteness can be
% tested by checking (ensuring) that the sign of all eigenvalues is negative.

% Computation of the Hessian matrix for both profit functions (F_A and F_B)
% will not be performed based on symbolic computational operations.
  
% Get the number of solutions to be tested for second order optimality.
N = size(Solutions,1);


% Initialize the output variables. Since each Hessian will be computed with
% respect to a two-variable function, the resulting matrices will be of
% size 2 x 2. Therefore, each Hessian matrix will be associated with 2
% eigenvalues.
HF_A_Eigs = zeros(N,2);
HF_B_Eigs = zeros(N,2);

% Loop through the various solutions in order to compute the final forms
% for the Hessian matrices and the corresponding eigenvalues.
for k = 1:1:N
    % Get the obtained optimal solutions per variable.
    T1_A_opt = Solutions(k,1);
    T2_A_opt = Solutions(k,2);
    T1_B_opt = Solutions(k,3);
    T2_B_opt = Solutions(k,4);
    % Get the obtained optimal solutions associated with Firm A and Firm B.
    TA_opt = [T1_A_opt,T2_A_opt];
    TB_opt = [T1_B_opt,T2_B_opt];
    % Set the function handles to the profit funtions F_A and F_B.
    F_A = @(TA)FirmANonLinearObjectiveFunction(TA,T1_B_opt,T2_B_opt,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);    
    F_B = @(TB)FirmBNonLinearObjectiveFunction(TB,T1_A_opt,T2_A_opt,Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2,Theta1,Theta2,P_A_1,P_A_2,P_B_1,P_B_2,M,K,C,Gamma);
    % Compute the Hessian function of F_A at the point TA_opt.
    HF_A = hessian(F_A,TA_opt);
    % Compute the Hessian function of F_B at the point TB_opt.
    HF_B = hessian(F_B,TB_opt);
    % Compute the corresponding eigenvalues.
    HF_A_Eigs(k,:) = eigs(HF_A);
    HF_B_Eigs(k,:) = eigs(HF_B);
end;

end

