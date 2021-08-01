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
SSj = subs(SSj,[Lj L_j],[1 1]);
SS_j = subs(SS_j,[Lj L_j],[1 1]);
% -------------------------------------------------------------------------

% Substitute the quantities (SSj) and (SS_j) in the corresponding
% expressions for (Qj) and (Q_j).
Qj = subs(Qj,[Sj S_j],[SSj SS_j]);
Q_j = subs(Q_j,[Sj S_j],[SSj SS_j]);

% Define the expressions for the profit functions Fj and F_j where the 
% associated cost functions will be linear.
Fj = Qj^2 - G*SSj;
F_j = Q_j^2 - G*SS_j;

% -------------------------------------------------------------------------
% Uncomment the following lines of code in order to enforce gamma_prime = 0.
Fj = subs(Fj,gamma_prime,0);
F_j = subs(F_j,gamma_prime,0);
% -------------------------------------------------------------------------

% Define the expressions for the revenues of each firm.
Rj = Qj^2;
R_j = Q_j^2;
% Compute the first order derivatives of the revenues Rj and R_j with
% respect to Tj and T_j respectively.
DRjTj = diff(Rj,Tj);
DR_jT_j = diff(R_j,T_j);
% Decompose the expressions DRjTj and DR_jT_j as fractions.
[DRjTj_n,DRjTj_d] = numden(DRjTj);
[DR_jT_j_n,DR_jT_j_d] = numden(DR_jT_j);
% Re-express the numerators of the above expressions as multivariate
% polynomials of Tj and T_j.
[DRjTj_n_c,DRjTj_n_t] = coeffs(DRjTj_n,[Tj T_j]);
[DR_jT_j_n_c,DR_jT_j_n_t] = coeffs(DR_jT_j_n,[Tj T_j]);
% Compute the first order derivatives of Fj and F_j with respect to Tj and
% T_j.
DFjTj = diff(Fj,Tj);
DF_jT_j = diff(F_j,T_j);
% Decompose the quantities DFjTj and DF_jT_j as fractions.
[DFjTj_n,DFjTj_d] = numden(DFjTj);
[DF_jT_j_n,DF_jT_j_d] = numden(DF_jT_j);
% Express the numerators of the above expressions as multivariate
% polynomials of Tj and T_j.
[DFjTj_n_c,DFjTj_n_t] = coeffs(DFjTj_n,[Tj T_j]);
[DF_jT_j_n_c,DF_jT_j_n_t] = coeffs(DF_jT_j_n,[Tj T_j]);


% Compute the second order derivatives of the profit functions Fj and F_j
% with respect to the corresponding strategic variables Tj and T_j.
D2FjTjTj = diff(diff(Fj,Tj),Tj);
D2F_jT_jT_j = diff(diff(F_j,T_j),T_j);
% Express the above quantities as fractions.
[D2FjTjTj_n,D2FjTjTj_d] = numden(D2FjTjTj);
[D2F_jT_jT_j_n,D2F_jT_jT_j_d] = numden(D2F_jT_jT_j);
% Taking into account that the denominators are strictly positive we can
% express the numerators as multi-variate polynomials of Tj and T_j.

% Obtain the best response functions Tj_star(T_j) and T_j_star(Tj).
Tj_star = solve(DFjTj_n==0,Tj);
T_j_star = solve(DF_jT_j_n==0,T_j);
% Try to obtain the equilibrium points (Tj_opt,T_j_opt) by taking into
% consideration the following points:
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

% Important note: Due the particular form of the first derivatives DFjTj 
% and DF_jT_j the quantities Tj_star and T_j_star actually hold the final
% solutions with respect to Tj and T_j.


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
% Re-express T_j_opt_eql as a fraction.
[T_j_opt_eql_n,T_j_opt_eql_d] = numden(T_j_opt_eql);
T_j_opt = solve(T_j_star==T_j,T_j);


% Compute the second order derivatives.
D2FjTjT_j = diff(diff(Fj,Tj),T_j);
D2F_jTjT_j = diff(diff(F_j,Tj),T_j);
% Re-express the above derivative as a fraction.
[D2FjTjT_j_n,D2FjTjT_j_d] = numden(D2FjTjT_j);
% Express the numerator of the fraction as a polynomial with respect to Tj.
[D2FjTjT_j_n_c,D2FjTjT_j_n_t] = coeffs(D2FjTjT_j_n,Tj);
% The quantity D2FjTjT_j_d is always positive.
% The quantity D2FjTjT_j_n is of the following form:
% D2FjTjT_j_d = ALPHA * Tj^2 + BETA * Tj + GAMMA.
% Compute the quantities ALPHA, BETA and GAMMA.
ALPHA = D2FjTjT_j_n_c(1);
BETA = D2FjTjT_j_n_c(2);
GAMMA = D2FjTjT_j_n_c(3);
% Simplify the above expressions.  
ALPHA = simplify(collect(expand(ALPHA)));
BETA = simplify(collect(expand(BETA)));
GAMMA = simplify(collect(expand(GAMMA)));
% Compute the discriminant of the polynomial.
DELTA = BETA^2 - 4*ALPHA*GAMMA;
DELTA = simplify(collect(expand(DELTA)));
% Compute the two roots of the polynomial.
Tj_pos = (-BETA + sqrt(DELTA)) / (2*ALPHA);
Tj_neg = (-BETA - sqrt(DELTA)) / (2*ALPHA);
% Simplify the above expressions.
Tj_pos = simplify(collect(expand(Tj_pos)));
Tj_neg = simplify(collect(expand(Tj_neg)));
% Rewrite the solutions of the polynomial as fractions.
[Tj_pos_n,Tj_pos_d] = numden(Tj_pos);
[Tj_neg_n,Tj_neg_d] = numden(Tj_neg);

% Compute first and second derivatives of (Qj) respect to (Tj) and (T_j).
DQjTj  = diff(Qj,Tj);
DQjT_j = diff(Qj,T_j);
D2QjTjTj = diff(DQjTj,Tj);
D2QjTjT_j = diff(DQjTj,T_j);
D2QjT_jT_j = diff(DQjT_j,T_j);
% Simplify the previously computed expressions.
DQjTj = simplify(collect(expand(DQjTj)));
DQjT_j = simplify(collect(expand(DQjT_j)));
D2QjTjTj = simplify(collect(expand(D2QjTjTj)));
D2QjTjT_j = simplify(collect(expand(D2QjTjT_j)));
D2QjT_jT_j = simplify(collect(expand(D2QjT_jT_j)));
% Re-write the above expressions in fractional form.
[DQjTj_n,DQjTj_d] = numden(DQjTj);
[DQjT_j_n,DQjT_j_d] = numden(DQjT_j);
[D2QjTjTj_n,D2QjTjTj_d] = numden(D2QjTjTj);
[D2QjTjT_j_n,D2QjTjT_j_d] = numden(D2QjTjT_j);
[D2QjT_jT_j_n,D2QjT_jT_j_d] = numden(D2QjT_jT_j);