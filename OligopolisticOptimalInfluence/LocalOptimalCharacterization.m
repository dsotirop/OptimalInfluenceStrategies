function [HF_A_Flags,HF_B_Flags,HF_Flags] = LocalOptimalCharacterization(HF_A_Eigs,HF_B_Eigs)


% The main purpose of this function is to characterize the kind of the
% local optima obtained during the optimization process. Actually, the 
% optimization is not performed per se but the optimal points are obtained 
% by solving the system of nonlinear equations which emerge by the First 
% Order Conditions of the two optimization problems. In this case, there 
% exist two profit functiobs F_A and F_B that are to be maximized. 

% Characterizatio of the obtained local optima points will be performed on 
% the basis of the associated eigenvalues of the corresponding Hessian matrices
% HF_A and and HF_B which are evaluated on each solution quadraplet. These
% values are stored in the N x 2 matrices HF_A and HF_B.

% Since profit functions F_A and F_B are optimized over the pairs of variables
% T1_A, T2_A and T1_B and T2_B, the corresponding Hessian matrices will be
% 2 x 2 matrices. Therefore, each Hessian matrix HF_A and HF_B will be
% associated with two eigenvalues. If both eigenvalues are negative for a
% given configuration of the rest of the problem variables then the
% associated Hessian matrix is negative definite and the local optimum is
% a maximum. Otherwise, if both eigenvalues are positive then the
% associated Hessian matrix is positive definite and the local optimum is a
% minimum.

% FLAGS INTERPRETATION:
% FLAG = 1  ==> Both Eigenvalues are Positive ==> Local Minimum.
% FLAG = -1 ==> Both Eigenvalues are Negative ==> Local Maximum.
% FLAG = 0 ==> One Eigenvalue is Positive and the other is Negetive.

% Get the number of obtained solutions.
N = size(HF_A_Eigs,1);

% Initialize flag matrices.
HF_A_Flags = zeros(N,1);
HF_B_Flags = zeros(N,1);

% Find the positions where the first and second eigenvalues of HF_A are 
% positive.
HF_A_1_pos = find(HF_A_Eigs(:,1)>0);
HF_A_2_pos = find(HF_A_Eigs(:,2)>0);
% Find the positions where the first and second eigenvalues aof HF_A are
% simultaneously positive.
HF_A_pos = intersect(HF_A_1_pos,HF_A_2_pos);

% Find the positions where the first and second eigenvalues of HF_A are 
% negative.
HF_A_1_neg = find(HF_A_Eigs(:,1)<0);
HF_A_2_neg = find(HF_A_Eigs(:,2)<0);
% Find the positions where the first and second eigenvalues of HF_A are
% simulataneously negative.
HF_A_neg = intersect(HF_A_1_neg,HF_A_2_neg);

% Update the HF_A_Flags values.
HF_A_Flags(HF_A_neg) = -1;
HF_A_Flags(HF_A_pos) = 1;


% Find the positions where the first and second eigenvalues of HF_B are 
% positive.
HF_B_1_pos = find(HF_B_Eigs(:,1)>0);
HF_B_2_pos = find(HF_B_Eigs(:,2)>0);
% Find the positions where the first and second eigenvalues aof HF_B are
% simultaneously positive.
HF_B_pos = intersect(HF_B_1_pos,HF_B_2_pos);

% Find the positions where the first and second eigenvalues of HF_A are 
% negative.
HF_B_1_neg = find(HF_B_Eigs(:,1)<0);
HF_B_2_neg = find(HF_B_Eigs(:,2)<0);
% Find the positions where the first and second eigenvalues of HF_A are
% simulataneously negative.
HF_B_neg = intersect(HF_B_1_neg,HF_B_2_neg);

% Update the HF_A_Flags values.
HF_B_Flags(HF_B_neg) = -1;
HF_B_Flags(HF_B_pos) = 1;

% Sum the output of both Flags.
% Mind that zero values in HF_Flags will be indicative of reversed minimum
% and maximum optima for the two underlying optimization problems.
HF_Flags = HF_A_Flags + HF_B_Flags;

end

