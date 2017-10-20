% This script file symbolically computes the optimal prices.
% Limiting beliefs vectors are assumed to be already computed.

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

% Store the expressions for the profit functions F_A and F_B before
% substituting the final expressions for the limiting beliefs XXA and XXB.
% Let these expressions be named as Fa and Fb.
Fa = F_A;
Fb = F_B;

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

% Compute the derivatives of F_A with respect to T1_A and and T2_A
DF_A_T1_A = diff(F_A,T1_A);
DF_A_T2_A = diff(F_A,T2_A);

% Compute the derivatives of F_B with respect to T1_B and T2_B.
DF_B_T1_B = diff(F_B,T1_B);
DF_B_T2_B = diff(F_B,T2_B);

% Simplify expressions for the derivatives of the profit functions F_A and
% F_B for products A and B respectively with respect to the corresponding
% optimization variables T1_A, T2_A and T1_B and T2_B.
DF_A_T1_A = simplify(DF_A_T1_A);
DF_A_T2_A = simplify(DF_A_T2_A);
DF_B_T1_B = simplify(DF_B_T1_B);
DF_B_T2_B = simplify(DF_B_T2_B);

% Compute the Hessian function of F_A with respect to T1_A and T2_A.
HF_A = hessian(F_A,[T1_A,T2_A]);
% Compute the Hessian function of F_B with respect to T1_B and T2_B.
HF_B = hessian(F_B,[T1_B,T2_B]);

% -------------------------------------------------------------------------
% OBJECTIVE FUNCTIONS SIMPLIFICATION ATTEMPT!!!
% -------------------------------------------------------------------------
% Performe some further symbolic operations on the expressions of Fa and Fb
% before substituting the final expressions for the limiting beliefs XXA
% (XA) and XXB (XB).

% 1st Step: Collect the expressions of Fa and Fb with respect to the
% combined set of variables XA and XB.
Fa = collect(Fa,[XA,XB]);
Fb = collect(Fb,[XA,XB]);

% 2nd Step: Decompose expressions of Fa and Fb into the corresponding
% numerators and denumerators such as:
%        Ra              Rb
% Fa = ------ and Fb = ------
%        Qa              Qb
[Ra,Qa] = numden(Fa);
[Rb,Qb] = numden(Fb);

% 3rd Step: Collect the expressions of Ra and Rb with respect to the
% combined set of variables XA and Xb.
Ra = collect(Ra,[XA,XB]);
Rb = collect(Rb,[XA,XB]);

% 4th Step: Decompose expressions XXA (XA) and XXB (XB) into the corresponding
% numerators and denumerators such as:
%        Ua              Ub
% XA = ------ and XB = ------
%        Va              Vb
[Ua,Va] = numden(XXA);
[Ub,Vb] = numden(XXB);

% 5th Step: Collect the expressions of Ua, Va and Ub, Vb with respect to
% the main optimization variables T1_A, T2_A, T1_B and T2_B.
Ua = collect(Ua,[T1_A, T2_A,T1_B,T2_B]);
Va = collect(Va,[T1_A, T2_A,T1_B,T2_B]);
Ub = collect(Ub,[T1_A, T2_A,T1_B,T2_B]);
Vb = collect(Vb,[T1_A, T2_A,T1_B,T2_B]);
% -------------------------------------------------------------------------
% The following code segment has been intentionally commented out.
% % This code segment needs to be revisited since the formulas of the
% % objectives functions have been modified.
% % First Order Conditions System Simplification.
% % Let Ga1, Ga2, Gb1 and Gb2 be the corresponding functions that define the
% % system of nonlinear equations which is to be solved such that:
% %                        Ga1 = 0
% %                        Ga2 = 0
% %                        Gb1 = 0
% %                        Gb2 = 0
% Ga1 = DF_A_T1_A;
% Ga2 = DF_A_T2_A;
% Gb1 = DF_B_T1_B;
% Gb2 = DF_B_T2_B;
% % It turns out that each one of the above functionals is of the following
% % form:
% %                        P1
% %                Ga1 = ------ - Gamma
% %                        Q1
% %
% %                        P2
% %                Ga2 = ------  - Gamma
% %                        Q2
% %                        
% %                        P3
% %                Gb1 = ------ - Gamma
% %                        Q3
% %
% %                        P4
% %                Gb2 = ------  - Gamma
% %                        Q4
% %
% % By assigning each fraction Pi / Qi to the variable Si, for i in {1,2,3,4}
% % we may obtain the Si expressios according to:
% S1 = solve(Ga1==0,Gamma);
% S2 = solve(Ga2==0,Gamma);
% S3 = solve(Gb1==0,Gamma);
% S4 = solve(Gb2==0,Gamma);
% % Decompose each fraction to the corresponding numerator and denumerator.
% [P1,Q1] = numden(S1);
% [P2,Q2] = numden(S2);
% [P3,Q3] = numden(S3);
% [P4,Q4] = numden(S4);
% % Check that all denumerators are equal (pairwise equality).
% % Form a matrix with all the Q's.
% Q = [Q1;Q2;Q3;Q4];
% % Initialize a 4x4 with all the possible pairwise differences.
% QDiffs = ones(4,4);
% for k = 1:1:4
%     for l = 1:1:4
%         QDiffs(k,l) = eval(Q(k)-Q(l));
%     end;
% end;
% Qsum = sum(sum(QDiffs));
% if(Qsum==0)
%     fprintf('All Q denumerators were found to be equal\n');
% end;
% % Let Qo be the common expression for all Q's.
% Qo = Q1;
% % Get coefficients and corresponding monomial terms for each one of the
% % quantities P1, P2, P3, P4 and Qo.
% [CP1,TP1] = coeffs(P1,[T1_A,T2_A,T1_B,T2_B]);
% [CP2,TP2] = coeffs(P2,[T1_A,T2_A,T1_B,T2_B]);
% [CP3,TP3] = coeffs(P3,[T1_A,T2_A,T1_B,T2_B]);
% [CP4,TP4] = coeffs(P4,[T1_A,T2_A,T1_B,T2_B]);
% [CQo,TQo] = coeffs(Qo,[T1_A,T2_A,T1_B,T2_B]);
% % Having demonstrated that Q1 = Q2 = Q3 = Q4 = Qo, the system of fractional
% % polynomial equations above can be reformulated as follows:
% %                       P1 = Qo * Gamma [1]
% %                       P2 = Qo * Gamma [2]
% %                       P3 = Qo * Gamma [3]
% %                       P4 = Qo * Gamma [4]
% % which entails that:    
% %                       P1 = P2 = P3 = P4 [5]
% % According to [5] we could reformulate the system of nonlinear equations 
% % [1],[2],[3],[4] in the following way:
% %                       R1 = P1 - P2 = 0 [6]
% %                       R2 = P1 - P3 = 0 [7]
% %                       R3 = P1 - P4 = 0 [8]
% % No other difference of the form (Pi - Pj) can be used as an additional 
% % equation in order to solve the system since it can be produced by
% % performing pairwise subtraction from the set of equations [6], [7] and
% % [8].                     
% % So, the last equation that needs to be takem into consideration is one of
% % the original ones ([1],[2],[3] and [4]). So the last equation of the
% % system could be the following:
% %                       R4 = P1 - Qo * Gamma = 0 [9]
% R1 = P1 - P2;
% R2 = P1 - P3;
% R3 = P1 - P4;
% R4 = P2 - P3;
% R5 = P2 - P4;
% R6 = P3 - P4;
% %R4 = P1 - Qo * Gamma;
% % An alternative system of equations may be build by considering the
% % following set of equations.
% W1 = P1 - Qo * Gamma;
% W2 = P2 - Qo * Gamma;
% W3 = P3 - Qo * Gamma;
% W4 = P4 - Qo * Gamma;