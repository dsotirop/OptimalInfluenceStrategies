% This script file evaluates the derivatives of all necessary quantities in 
% order to apply the Envelop Theorem for the two-player continuous game
% which underlies the Simplified Oligopolistic Optimal Influence model.

% We need to determine the behaviour of the Nash Equilibrium point 
% p = (Tj,T_j) as a function of the exogenous parameters of the model such
% as: Poj, Po_j, Lj L_j, K, M, C, and G. Thus, if we let 
% W in {Poj, Po_j, Lj, L_j, K, M, C,G}, we need to determine the behaviour of 
% p(W) = (Tj(W),T_j(W)).

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
                            
% Declare the additional symbolic variables dTA and dTB.
syms dTj dT_j

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
% Substitute the quantities (SSj) and (SS_j) in the corresponding
% expressions for (Qj) and (Q_j).
Qj = subs(Qj,[Sj S_j],[SSj SS_j]);
Q_j = subs(Q_j,[Sj S_j],[SSj SS_j]);

% Set the expressions for the payoff functions Fj and F_j with respect to
% the above parameters.
Fj = Qj^2 - G*SSj^2;
F_j = Q_j^2 - G*SS_j^2;

% Compute the first order derivatives with respect to Tj and T_j
% accordingly.
Dj = diff(Fj,Tj); %Dj = DFjTj
D_j = diff(F_j,T_j); %D_j = DF_jT_j

% If DTjW (--> dTj) is the total derivative of Tj with respect to parameter 
% w and DT_jW (--> dT_j)is the total derivative of T_j with respect to  
% parameter w then the following equations hold:
%                               Aj *  dTj  + Bj *  dT_j = - Gj [1]
%                               A_j * dTj  + B_j * dT_j = - G_j [2]
% where:
%        Aj =  DDjTj   = D2FjTjTj    [3]
%        Bj =  DDjT_j  = D2FjTjT_j   [4]
%        Gj =  DDjW    = D2FjTjW     [5]
%        A_j = DD_jTj  = D2F_jTjT_j  [6]
%        B_j = DD_jT_j = D2F_jT_jT_j [7]
%        G_j = DD_jW   = D2F_jT_jW   [8]

% Compute the quantities Aj and Bj as well as the quantities A_j and B_j.
Aj = diff(Dj,Tj);
Bj = diff(Dj,T_j);
A_j = diff(D_j,Tj);
B_j = diff(D_j,T_j);

% Compute the quantities Gj and G_j where the parameter W is Poj.
W = Poj;
Gj = diff(Dj,W);
G_j = diff(D_j,W);

% Re-express the previously defined quantities as fractions.
[NAj,DAj] = numden(Aj);
[NBj,DBj] = numden(Bj);
[NGj,DGj] = numden(Gj);

[NA_j,DA_j] = numden(A_j);
[NB_j,DB_j] = numden(B_j);
[NG_j,DG_j] = numden(G_j);

% Define the system of linear equations with respect to dTj and dT_j.
S_Pj = [Aj * dTj + Bj * dT_j + Gj==0;A_j * dTj + B_j * dT_j + G_j==0];
% Solve the linear system with respect to dTA and dTB.
S_Pj_solution = solve(S_Pj,[dTj dT_j]);

% Provide an alternative solution by solving the system of linear equations
% with respect to the auxiliary variables Aj, Bj, Gj, A_j, B_j and G_j.

% Declare the auxiliary symbolic variables.
syms AAj BBj GGj AA_j BB_j GG_j

R = [AAj BBj;AA_j BB_j];
Q = [-GGj;-GG_j];
So = R\Q;

% Keep in mind that the exact value of So is given by the following
% equation:
%        |(BBj*GG_j - BB_j*GGj)/(AAj*BB_j - AA_j*BBj) |   |dTj |   |Wj / Rj  |
% So  =  |                                            | = |    | = |         | [9]
%        |-(AAj*GG_j - AA_j*GGj)/(AAj*BB_j - AA_j*BBj)|   |dT_j|   |W_j / R_j|
%
% According to the previous definitions we may write that:
% Wj = Bj * G_j - B_j * Gj = BGj_j - BG_jj            [10]
% W_j = -(Aj * G_j - A_j * Gj) = -(AGj_j - AG_jj)     [11]
% Rj = R_j = R = Aj * B_j  - A_j * Bj = ABj_j - AB_jj [12]

% Compute the previously defined intermediate quantities:
ABj_j = Aj * B_j;
AB_jj = A_j * Bj;
BGj_j = Bj * G_j;
BG_jj = B_j * Gj;
AGj_j = Aj * G_j;
AG_jj = A_j * Gj;

% Decompose the above quantities into numerator - denumerator fractions.
[NABj_j,DABj_j] = numden(ABj_j);
[NAB_jj,DAB_jj] = numden(AB_jj);
[NBGj_j,DBGj_j] = numden(BGj_j);
[NBG_jj,DBG_jj] = numden(BG_jj);
[NAGj_j,DAGj_j] = numden(AGj_j);
[NAG_jj,DAG_jj] = numden(AG_jj);

% Solve the systems of differential equations.
Soo = subs(So,[AAj BBj GGj AA_j BB_j GG_j],[Aj Bj Gj A_j B_j G_j]);
Soo = simplify(Soo);
% syms x(Poj) y(Poj)
% Soo = subs(Soo,[Tj T_j],[x y]);
% eqn1 = diff(x,Poj) == Soo(1);
% eqn2 = diff(y,Poj) == Soo(2);
% eqns = [eqn1;eqn2];
% Doo = dsolve(eqns);