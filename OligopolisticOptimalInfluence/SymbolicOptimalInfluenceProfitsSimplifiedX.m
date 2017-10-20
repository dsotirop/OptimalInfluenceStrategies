% This script file symbolically computes the firms' objective functions
% F_A and F_B under the following assumptions:
% (i):   Lambda_A_1 = Lambda_A_2 = Lambda_B_1 = Lambda_B_2 = Lambda.
% (ii):  Theta1 = Theta2 = Theta.
% (iii): P_A_1 = P_A_2 = P_B_1 = P_B_2 = P.
% Limiting beliefs vectors are assumed to be already computed through 
% calling the script file "SymbolicOptimalInfluenceBeliefs.m" which defines
% the additional symbolic variables.

% Clear workspace and command window.
clc
clear all

% Setup symbolic variables.
syms X_A_1 X_A_2 XA
syms X_B_1 X_B_2 XB
syms Q_A_1 Q_A_2 Q_A
syms Q_B_1 Q_B_2 Q_B
syms M K %These are the sensitivity coefficients.
syms pA pB % These are the prices for products A and B.
syms C % Marginal cost.
syms Gamma % Marginal influence cost.

% Setup symbolic variables concerning the network structure.
syms T1_A T1_B
syms T2_A T2_B

% Setup symbolic variables required for simplification.
syms Lambda Theta P

% Define demand functions for consumers 1 and 2. 
Q_A_1 = X_A_1 - pA + M*pB - K*X_B_1;
Q_B_1 = X_B_1 - pB + M*pA - K*X_A_1;

Q_A_2 = X_A_2 - pA + M*pB - K*X_B_2;
Q_B_2 = X_B_2 - pB + M*pA - K*X_A_2;

% Compute demand functions for products A and B.
Q_A = Q_A_1 + Q_A_2;
Q_B = Q_B_1 + Q_B_2;

% Define profit functions for products A and B.
% F_A = Q_A * (pA - C) - Gamma*(T1_A + T2_A);
% F_B = Q_B * (pB - C) - Gamma*(T1_B + T2_B);

F_A = Q_A * (pA - C) - Gamma*(T1_A^2 + T2_A^2);
F_B = Q_B * (pB - C) - Gamma*(T1_B^2 + T2_B^2);


% Compute first derivative of F_A with respect to pA.
DF_A_pA = diff(F_A,pA);
% Solve DF_A_pA==0 with respect to pA (best responce function).
pA_star = solve(DF_A_pA==0,pA);

% Compute first derivative of F_B with respect to pB.
DF_B_pB = diff(F_B,pB);
% Solve DF_A_pB==0 with respect to pB (best responce function).
pB_star = solve(DF_B_pB==0,pB);

% Equate X_A_1 and X_A_2 to X_A and X_B_1 and X_B_2 to X_B
pA_star = subs(pA_star,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);
pB_star = subs(pB_star,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);

pA_star_copy = pA_star;
pB_star_copy = pB_star;

pA_star = subs(pA_star,pB,pB_star_copy);
pB_star = subs(pB_star,pA,pA_star_copy);

pA_opt = solve(pA==pA_star,pA);
pB_opt = solve(pB==pB_star,pB);

% Substitute variables X_A_1 and X_A_2 with XA and X_B_1 and X_B with XB
% inside the expressions for Q_A_1, Q_A_2, Q_B_1 and Q_B_2.
Q_A_1 = subs(Q_A_1,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);
Q_A_2 = subs(Q_A_2,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);
Q_B_1 = subs(Q_B_1,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);
Q_B_2 = subs(Q_B_2,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);

% Substitute variables pA and pB with pA_opt and pB_opt respectively in the
% expressions for Q_A_1, Q_A_2, Q_B_1 and Q_B_2.
Q_A_1_opt = subs(Q_A_1,[pA,pB],[pA_opt,pB_opt]);
Q_A_2_opt = subs(Q_A_2,[pA,pB],[pA_opt,pB_opt]);
Q_B_1_opt = subs(Q_B_1,[pA,pB],[pA_opt,pB_opt]);
Q_B_2_opt = subs(Q_B_2,[pA,pB],[pA_opt,pB_opt]);

% Simplify expressions for Q_A_1_opt, Q_A_2_opt, Q_B_1_opt and Q_B_2_opt.
Q_A_1_opt = simplify(collect(expand(Q_A_1_opt)));
Q_A_2_opt = simplify(collect(expand(Q_A_2_opt)));
Q_B_1_opt = simplify(collect(expand(Q_B_1_opt)));
Q_B_2_opt = simplify(collect(expand(Q_B_2_opt)));

% Compute optimal demand functions for products A and B.
Q_A_opt = Q_A_1_opt + Q_A_2_opt;
Q_B_opt = Q_B_1_opt + Q_B_2_opt;

% Substitute variables X_A_1 and X_A_2 with XA and X_B_1 and X_B with XB
% inside the expressions for F_A and F_B.
F_A = subs(F_A,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);
F_B = subs(F_B,[X_A_1,X_A_2,X_B_1,X_B_2],[XA,XA,XB,XB]);

% Substitute the optimal expressions for quantities and prices in the
% profit functions F_A and F_B for the products A and B respectively.
F_A = subs(F_A,[pA,pB,Q_A,Q_B],[pA_opt,pB_opt,Q_A_opt,Q_B_opt]);
F_B = subs(F_B,[pA,pB,Q_A,Q_B],[pA_opt,pB_opt,Q_A_opt,Q_B_opt]);

% Compute Limiting beliefs by executing script
% GeneralOptimalInfluenceBeliefs.
SymbolicOptimalInfluenceBeliefs;

% Reduce the complexity of terms XXA and XXB before substituting in
% expressions for F_A and F_B by utilizing the collect() function.
XXA = collect(XXA);
XXB = collect(XXB);

% Substitute the agents' limiting beliefs per product (XA and XB) with the
% corresponding expressions XXA and XXB that have been computed within the
% GeneralOptimalInfluenceBeliefs script inside the expressions for the 
% optimal profits F_A and F_B.
F_A = subs(F_A,[XA,XB],[XXA,XXB]);
F_B = subs(F_B,[XA,XB],[XXA,XXB]);

% Substitute the agents' limiting beliefs per profuct (XA and XB) with the
% corresponding expressions XXA and XXB that have been computed within the
% GeneralOptimalInfluenceBeliefs script inside the expressions for the
% optimal price levels pA_opt and pB_opt.
pA_opt = subs(pA_opt,[XA,XB],[XXA,XXB]);
pB_opt = subs(pB_opt,[XA,XB],[XXA,XXB]);

% Substitute the agents' limiting beliefs per profuct (XA and XB) with the
% corresponding expressions XXA and XXB that have been computed within the
% GeneralOptimalInfluenceBeliefs script inside the expressions for the
% optimal price levels Q_A_opt and Q_B_opt.
Q_A_opt = subs(Q_A_opt,[XA,XB],[XXA,XXB]);
Q_B_opt = subs(Q_B_opt,[XA,XB],[XXA,XXB]);

% -------------------------------------------------------------------------
% SIMPLIFY EXPRESSIONS FOR THE OPTIMAL INFLUENCE LEVELS AT EQUILIBRIUM
% -------------------------------------------------------------------------
% Enforce internal parameters equalities for the quantities SA,S1,S2 and SB
% indicated by the point (i).
S = subs(S,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
% Enforce internal parameter equalities for the quantities XXA and XXB 
% indicated by the point (ii).
S = subs(S,[Theta1,Theta2],[Theta,Theta]);
% Enforce internal parameter equalities for the quantities XXA and XXB 
% indicated by the point (iii).
S = subs(S,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);
% Further simplification of the previously computed quantities.
S = simplify(collect(expand(S)));

% -------------------------------------------------------------------------
% SIMPLIFY EXPRESSIONS FOR THE OPTIMAL BELIEFS AT EQUILIBRIUM
% -------------------------------------------------------------------------
% Enforce internal parameters equalities for the quantities XXA and XXB
% indicated by the point (i).
XXA = subs(XXA,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
XXB = subs(XXB,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
% Enforce internal parameter equalities for the quantities XXA and XXB 
% indicated by the point (ii).
XXA = subs(XXA,[Theta1,Theta2],[Theta,Theta]);
XXB = subs(XXB,[Theta1,Theta2],[Theta,Theta]);
% Enforce internal parameter equalities for the quantities XXA and XXB 
% indicated by the point (iii).
XXA = subs(XXA,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);
XXB = subs(XXB,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);
% Further simplification of the previously computed quantities.
XXA = simplify(collect(expand(XXA)));
XXB = simplify(collect(expand(XXB))); 

% -------------------------------------------------------------------------
% SIMPLIFY EXPRESSIONS FOR THE OPTIMAL PRICE LEVELS AT EQUILIBRIUM
% -------------------------------------------------------------------------
% Enforce internal parameter equalities for the quantities F_A and F_B 
% indicated by the point (i).
pA_opt = subs(pA_opt,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
pB_opt = subs(pB_opt,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
% Enforce internal parameter equalities for the quantities F_A and F_B 
% indicated by the point (ii).
pA_opt = subs(pA_opt,[Theta1,Theta2],[Theta,Theta]);
pB_opt = subs(pB_opt,[Theta1,Theta2],[Theta,Theta]);
% Enforce internal parameter equalities for the quantities F_A and F_B 
% indicated by the point (iii).
pA_opt = subs(pA_opt,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);
pB_opt = subs(pB_opt,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);

% -------------------------------------------------------------------------
% SIMPLIFY EXPRESSIONS FOR THE OPTIMAL PROFITS AT EQUILIBRIUM
% -------------------------------------------------------------------------
% Enforce internal parameter equalities for the quantities F_A and F_B 
% indicated by the point (i).
F_A = subs(F_A,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
F_B = subs(F_B,[Lambda_A_1,Lambda_A_2,Lambda_B_1,Lambda_B_2],[Lambda,Lambda,Lambda,Lambda]);
% Enforce internal parameter equalities for the quantities F_A and F_B 
% indicated by the point (ii).
F_A = subs(F_A,[Theta1,Theta2],[Theta,Theta]);
F_B = subs(F_B,[Theta1,Theta2],[Theta,Theta]);
% Enforce internal parameter equalities for the quantities F_A and F_B 
% indicated by the point (iii).
F_A = subs(F_A,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);
F_B = subs(F_B,[P_A_1,P_A_2,P_B_1,P_B_2],[P,P,P,P]);

% -------------------------------------------------------------------------
% The expressions for the optimal profits for both firms (F_A and F_B) 
% expressed in terms of the underlying fundamental parameters may be
% decomposed in the following forms:
% F_A = F_A_Rev - F_A_Cost
% F_B = F_B_Rev - F_B_Cost
% such that:
% F_A_Cost = Gamma*(T1_A^2 + T2_A^2)
% F_B_Cost = Gamma*(T1_B^2 + T2_B^2)

% Define the symbolic variables for the expressions F_A_Cost and F_B_Cost
F_A_Cost = Gamma*(T1_A^2 + T2_A^2);
F_B_Cost = Gamma*(T1_B^2 + T2_B^2);

% Define the symbolic variables for the expressions F_A_Rev and F_B_Rev.
F_A_Rev = F_A + F_A_Cost;
F_B_Rev = F_B + F_B_Cost;

% Decompose the symbolic variables for the revenues F_A_Rev and F_B_Rev
% into the corresponding numerator denumerator expressions.
[DF_A_Rev,NF_A_Rev] = numden(F_A_Rev);
[DF_B_Rev,NF_B_Rev] = numden(F_B_Rev);

% Compute the derivatives of the expressions F_A and F_B with respect to
% the primary optimization variables {T1_A,T2_A} and {T1_B,T2_B}.
DF_A_T1_A = diff(F_A,T1_A);
DF_A_T2_A = diff(F_A,T2_A);
DF_B_T1_B = diff(F_B,T1_B);
DF_B_T2_B = diff(F_B,T2_B);
% Re-express the above derivatives-related quantities as fractions of the
% form:
%                           PA1
%            DF_A_T1_A = ---------
%                           QA1
%
%                           PA2
%            DF_A_T2_A = ---------
%                           QA2
%
%                           PB1
%            DF_B_T1_B = ---------
%                           QB1
%
%                           PB2
%            DF_B_T2_B = ---------
%                           QB2
[PA1,QA1] = numden(DF_A_T1_A);
[PA2,QA2] = numden(DF_A_T2_A);
[PB1,QB1] = numden(DF_B_T1_B);
[PB2,QB2] = numden(DF_B_T2_B);
% We need to compute the Q's differences in order to ensure that all
% derivatives-related quantities share the same denumerator.
Q = [QA1;QA2;QB1;QB2];
LQ = length(Q);
DQ = ones(LQ,LQ);
for i = 1:1:LQ
    for j = 1:1:LQ
        DQ(i,j) = Q(i) - Q(j);
    end;
end;
SQ = sum(sum(DQ));
if(SQ==0)
    fprintf('All derivative-related denumerators are equal\n');
end;

% Since the equilibrium point for this particular instance of the
% Oligopolistic game may be obtained by setting all the above
% derivative-related quantities to zero, we need to express the numerators
% as polynomials with respect to the corresponding primary optimization
% variables.
[c_a1,t_a1] = coeffs(PA1,T1_A);
[c_a2,t_a2] = coeffs(PA2,T2_A);
[c_b1,t_b1] = coeffs(PB1,T1_B);
[c_b2,t_b2] = coeffs(PB2,T2_B);