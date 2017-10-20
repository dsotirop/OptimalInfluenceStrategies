% This script file symbolically computes the optimal influences and
% beliefs.

%-------------------------------------------------------------------------
%                           IMPORTANT NOTE!!!
%-------------------------------------------------------------------------
% Mind that this script wil be called inside the scipt:
% GeneralOptimalInfluencePricesQuantities. Therefore, network-structure 
% related symbolic variables will be defined in the latter script.
% Also, mind that variables XXA and XXB will be substituting variables XA
% and XB.
%------------------------------------------------------------------------- 

% Setup symbolic variables.
syms Lambda_A_1 Lambda_A_2 Lambda_A
syms Lambda_B_1 Lambda_B_2 Lambda_B

% The following two lines of code should be uncommented in case this script
% is not to be executed inside the GeneralOptimalInfluencePricesQuantities 
% script.
% syms T1_A T1_B
% syms T2_A T2_B

syms Theta1 Theta2
syms VSA VS1 VS2 VSB
syms P_A_A P_A_B
syms P_A_1 P_A_2
syms P_B_B P_B_A
syms P_B_1 P_B_2

% Setup social interaction matrix components.
T = [1-Lambda_A_1-Lambda_A_2 Lambda_A_1 Lambda_A_2 0;...
     T1_A 1-T1_A-T1_B-Theta2 Theta2 T1_B;...
     T2_A Theta1 1-T2_A-T2_B-Theta1 T2_B;...
     0 Lambda_B_1 Lambda_B_2 1-Lambda_B_1-Lambda_B_2];
 
% Setup the components of the final influence vector. 
VS = [VSA VS1 VS2 VSB];
% Setup the corresponding system initiating from the
% left eigenvalues computation problem.
Y = VS * (T - eye(4));
System1 = [Y(1)==0;Y(2)==0;Y(3)==0;VSA+VS1+VS2+VSB==1];
% Solve the corresponding linear system in order to get the limiting
% influence vectors. (The will be equal for both products!!!).
Solution1 = solve(System1,VSA,VS1,VS2,VSB);
SA = Solution1.VSA;
S1 = Solution1.VS1;
S2 = Solution1.VS2;
SB = Solution1.VSB;
S = [SA;S1;S2;SB];
% Check that the sum of S components equals 1.
Ssum = simplify(collect(expand(sum(S))));
if(Ssum==1)
    fprintf('Elements of S sum up to 1\n');
end;
% Accumulate initial beliefs for product A.
PA = [P_A_A P_A_1 P_A_2 P_A_B];
% Accumulate initial beliefs for product B.
PB = [P_B_A P_B_1 P_B_2 P_B_B];
% Compute the limiting belief for product A.
XXA = PA * S;
% Compute the limiting belief for product A.
XXB = PB * S;
% Simplify expressions for XA and XB when P_A_A = 1, P_A_B = 0, P_B_A  = 0
% and P_B_B = 1.
XXA = subs(XXA,[P_A_A,P_A_B,P_B_A,P_B_B],[1 0 0 1]);
XXB = subs(XXB,[P_A_A,P_A_B,P_B_A,P_B_B],[1 0 0 1]);
