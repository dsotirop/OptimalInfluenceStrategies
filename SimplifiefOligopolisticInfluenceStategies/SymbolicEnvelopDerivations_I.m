% This script file evaluates the derivatives of all necessary quantities in 
% order to apply the Envelop Theorem for the two-player continuous game
% which underlies the Simplified Oligopolistic Optimal Influence model.
clc
clear
% Define fundamental symbolic quantities which will be initially treated as
% constants.
syms alpha beta gamma_prime
% Auxiliary parameters alpha, beta are given by
%
%           K*M - 2                    2*K - M                     
% alpha = ----------- [1] and beta = ----------- [2] 
%           M^2 - 4                    M^2 - 4                  
%

% Define fundamental symbolic quantities of the optimization problem for
% the two firms (j) and (_j).
syms Sj(Tj,T_j) S_j(Tj,T_j) % Sj and S_j denote the limiting influences.
syms Xj X_j                 % Xj and X_j denote the limiting beliefs.
syms Poj Po_j               % Poj and Po_j denote the initial beliefs.
syms Lj L_j                 % Direct influences exerted by the consumer 
                            % the towards two firms.
syms G                      % Marginal influence cost or Gamma.
                            
% Define the expressions for the equlibrium quantities (Qj) and (Q_j) at
% the first stage of backward induction procedure. Mind that the initial
% expressions for the quantities have been modified in order to ensure that
% they are non-negative quantities.
Qj  = alpha * Xj + beta * X_j - gamma_prime;
Q_j = alpha * X_j + beta * X_j - gamma_prime;

% Symbolic quantities (Xj) and (X_j) will be replaced by their equivalent
% expressions (XXj) and (XX_j).
XXj  =  Poj  + (1-Poj)  * Sj  - Poj  * S_j;
XX_j =  Po_j + (1-Po_j) * S_j - Po_j * Sj;
% Substitute quantities (XXj) and (XX_j) in the corresponding expressions
% for (Qj) and (Q_j).
Qj  =  subs(Qj,[Xj,X_j],[XXj,XX_j]);
Q_j =  subs(Q_j,[Xj,X_j],[XXj,XX_j]);
% Express quantities Qj and Q_j as multi-variate polynomials of Sj and S_j.
[CQj,TQj]   =   coeffs(Qj,[Sj,S_j]);
[CQ_j,TQ_j] =   coeffs(Q_j,[Sj,S_j]);

% Define the quantities (SSj) and (SS_j) that will provide the expressions   
% for the limiting influences as functions of the investment levels Tj and 
% T_j.
SSj = L_j * Tj / (Lj*L_j + L_j*Tj + Lj*T_j);
SS_j = Lj * T_j / (Lj*L_j + L_j*Tj + Lj*T_j);

% -------------------------------------------------------------------------
% Uncomment the following block of code in order to enforce Lj = L_j = 1.
% SSj = subs(SSj,[Lj L_j],[1 1]);
% SS_j = subs(SS_j,[Lj L_j],[1 1]);
% -------------------------------------------------------------------------

% Substitute quantities (SSj) and (SS_j) in the corresponding
% expressions for (Qj) and (Q_j).
Qj = subs(Qj,[Sj S_j],[SSj SS_j]);
Q_j = subs(Q_j,[Sj S_j],[SSj SS_j]);

% Substitute quantities (SSj) and (SS_j) in the corresponding expressions
% for (XXj) and (XX_j).
XXj = subs(XXj,[Sj S_j],[SSj SS_j]);
XX_j = subs(XX_j,[Sj S_j],[SSj SS_j]);


% Define the expressions for the profit functions Fj and F_j.
Fj = Qj^2 - G*SSj^2;
F_j = Q_j^2 - G*SS_j^2;


% -------------------------------------------------------------------------
% Uncomment the following lines of code in order to enforce gamma_prime = 0.
% Fj = subs(Fj,gamma_prime,0);
% F_j = subs(F_j,gamma_prime,0);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Uncomment the following lines of code in order to enforce beta = 0.
% Fj = subs(Fj,beta,0);
% F_j = subs(F_j,beta,0);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Uncomment the following lines of code in order to enforce alpha = 0.
% Fj = subs(Fj,alpha,0);
% F_j = subs(F_j,alpha,0);
% -------------------------------------------------------------------------

% Compute the first order derivatives of Fj and F_j with respect to Tj and
% T_j.
DFjTj = diff(Fj,Tj);
DF_jT_j = diff(F_j,T_j);
% Decompose the quantities DFjTj and DF_jT_j as fractions.
[DFjTj_n,DFjTj_d] = numden(DFjTj);
[DF_jT_j_n,DF_jT_j_d] = numden(DF_jT_j);
% Quantities DFjTj_d and DF_jT_j_d are always positive.
% Express the numerators of the above expressions as multivariate
% polynomials of T_j and Tj respectively.
[DFjTj_n_c,DFjTj_n_t] = coeffs(DFjTj_n,T_j);
[DF_jT_j_n_c,DF_jT_j_n_t] = coeffs(DF_jT_j_n,Tj);


% Obtain the best response functions Tj_star(T_j) and T_j_star(Tj).
Tj_star = solve(DFjTj_n==0,Tj);
T_j_star = solve(DF_jT_j_n==0,T_j);
% Try to obtain the equilibrium points (Tj_opt,T_j_opt) by taking into
% consideration the following clarifications:
% Tj_star = fj(T_j) : expresses the value of the best response Tj given as
%                     a function of T_j.
% T_j_star = f_j(Tj): expresses the value of the best response T_j given as
%                     a function of Tj.
%
% In this setting: 
% (i)  the equilibrium value for Tj (Tj_opt) may be obtained by solving the 
%      following equation: Tj = fj(f_j(Tj)) ==> Tj_opt = ...
%  
% (ii) the equilibrium value for T_j (T_j_opt) may be obtained by solving 
%      the following equation: T_j = f_j(fj(T_j)) ==> T_j_opt = ...  
%
% Make copies of the original variables in order to perform the necessary
% substitutions.
% Tj_star_copy provides the best response function for Tj as:
% Tj_star = BR(T_j).
% T_j_star_copy provides the best response function for T_j as:
% T_j_star = BR(Tj);
Tj_star_copy = Tj_star; 
T_j_star_copy = T_j_star;
% Express quantities Tj_star_copy and T_j_star_copy as fractions.
[Tj_star_copy_n,Tj_star_copy_d] = numden(Tj_star_copy);
[T_j_star_copy_n,T_j_star_copy_d] = numden(T_j_star_copy);
% Express numerator and denominator of Tj_star_copy as polynomials of T_j.
[Tj_star_copy_n_c,Tj_star_copy_n_t] = coeffs(Tj_star_copy_n,T_j);
[Tj_star_copy_d_c,Tj_star_copy_d_t] = coeffs(Tj_star_copy_d,T_j);
% Expreess numerator and denominator of T_j_star_copy as polynomials of Tj.
[T_j_star_copy_n_c,T_j_star_copy_n_t] = coeffs(T_j_star_copy_n,Tj);
[T_j_star_copy_d_c,T_j_star_copy_d_t] = coeffs(T_j_star_copy_d,Tj);
% Express variable T_j in the fj(T_j) expression of Tj_star as f_j(Tj).
Tj_star = subs(Tj_star,T_j,T_j_star_copy);
% Solve the equation Tj = fj(f_j(Tj)) with respect to Tj and assign the
% solution value to to Tj_opt.
Tj_opt_eql = Tj_star - Tj;
% Re-express Tj_opt_eql as a fraction.
[Tj_opt_eql_n,Tj_opt_eql_d] = numden(Tj_opt_eql);
% Express the numerator of the expression Tj_opt_eql as a polynomial of Tj.
[Tj_opt_eql_n_c,Tj_opt_eql_n_t] = coeffs(Tj_opt_eql_n,Tj);
[Tj_opt_eql_d_c,Tj_opt_eql_d_t] = coeffs(Tj_opt_eql_d,Tj);
Tj_opt = solve(Tj_star==Tj,Tj);
% Express variable Tj in the f_j(Tj) expression of T_j_star as fj(T_j).
T_j_star = subs(T_j_star,Tj,Tj_star_copy);
% Solve the equation T_j = f_j(fj(T_j)) with respect to T_j and assign the
% solution value to the T_j_opt.
T_j_opt_eql = T_j_star - T_j;
% Re-express T_j_opt_eql as a fraction.
[T_j_opt_eql_n,T_j_opt_eql_d] = numden(T_j_opt_eql);
% Express the numerator of the expression T_j_opt_eql as a polynomial of
% T_j.
[T_j_opt_eql_n_c,T_j_opt_eql_n_t] = coeffs(T_j_opt_eql_n,T_j);
[T_j_opt_eql_d_c,T_j_opt_eql_d_t] = coeffs(T_j_opt_eql_d,T_j);
% Re-express T_j_opt_eql as a fraction.
[T_j_opt_eql_n,T_j_opt_eql_d] = numden(T_j_opt_eql);
T_j_opt = solve(T_j_star==T_j,T_j);
% Obtain the NE solutions by determining the roots of the numerator of
% Tj_opt_eql and T_j_opt_eql.
Tj_opt_sols = roots(Tj_opt_eql_n_c);
T_j_opt_sols = roots(T_j_opt_eql_n_c);

% Compute the second order derivatives of Fj and F_j with respect to Tj and
% T_j accordingly.
D2FjTjTj = diff(diff(Fj,Tj),Tj);
D2F_jT_jT_j = diff(diff(F_j,T_j),T_j);
% Re-express the above derivatives as fractions.
[D2FjTjTj_n,D2FjTjTj_d] = numden(D2FjTjTj);
[D2F_jT_jT_j_n,D2F_jT_jT_j_d] = numden(D2F_jT_jT_j);
% Simplify the expressions for the numerators.
D2FjTjTj_n = simplify(collect(expand(D2FjTjTj_n)));
D2F_jT_jT_j_n = simplify(collect(expand(D2F_jT_jT_j_n)));

% Compute the second order derivatives with respect to the other firm's 
% investment level.
D2FjTjT_j = diff(diff(Fj,Tj),T_j);
D2F_jTjT_j = diff(diff(F_j,Tj),T_j);
% Re-express the above derivatives as fractions.
[D2FjTjT_j_n,D2FjTjT_j_d] = numden(D2FjTjT_j);
[D2F_jTjT_j_n,D2F_jTjT_j_d] = numden(D2F_jTjT_j);
% Simplify the expressions for the numerators.
D2FjTjT_j_n = simplify(collect(expand(D2FjTjT_j_n)));
D2F_jTjT_j_n = simplify(collect(expand(D2F_jTjT_j_n)));

% Compute the second order derivatives with respect to Poj.
D2FjTjPoj = diff(diff(Fj,Tj),Poj);
D2F_jT_jPoj = diff(diff(F_j,T_j),Poj);
% Re-express the above derivatives as fractions.
[D2FjTjPoj_n,D2FjTjPoj_d] = numden(D2FjTjPoj);
[D2F_jT_jPoj_n,D2F_jT_jPoj_d] = numden(D2F_jT_jPoj);
% Simplify the expressions for the numerators.
D2FjTjPoj_n = simplify(collect(expand(D2FjTjPoj_n)));
D2F_jT_jPoj_n = simplify(collect(expand(D2F_jT_jPoj_n)));

% -------------------------------------------------------------------------
% Explicitly express the Best Response Functions for the two firms such
% that: 
%       (i):  Tj_star = Rj(T_j) 
%       (ii): T_j_star = R_j(Tj)
Rj = Tj_star_copy;
R_j = T_j_star_copy;
% Express the above quantities as fractions.
[Rj_n,Rj_d] = numden(Rj);
[R_j_n,R_j_d] = numden(R_j);
% We know that Rj_n and Rj_d are polynomial expressions of T_j.
[Rj_n_c,Rj_n_t] = coeffs(Rj_n,T_j);
[Rj_d_c,Rj_d_t] = coeffs(Rj_d,T_j);
% We also know that R_j_n and R_j_d are polynomial expressions of Tj.
[R_j_n_c,R_j_n_t] = coeffs(R_j_n,Tj);
[R_j_d_c,R_j_d_t] = coeffs(R_j_d,Tj);
% -------------------------------------------------------------------------