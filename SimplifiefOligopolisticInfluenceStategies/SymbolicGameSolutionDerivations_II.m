% This script file performs auxiliary symbolic calculations in order to
% obtain an analytical solution for the continous game underlying the
% Simplified Oligopolistic Optimal Influence model. Experiment with
% alternative cost function formulations.

clc
clear

% Define fundamental symbolic variables.
syms alpha beta gamma_prime G
syms Poj Po_j Po
syms Sj S_j So
syms Lj L_j

% Define the expressions for the limiting beliefs Xj and X_j.
Xj = Poj + (1 - Poj)*Sj - Poj*S_j;
X_j = Po_j + (1 - Po_j)*S_j - Po_j*Sj;

% Define the expressions for the quantities Qj and Q_j.
QQj = alpha * Xj + beta * X_j - gamma_prime;
QQ_j = alpha * X_j + beta * Xj - gamma_prime;

% The expressions Qj and Q_j may be alternatively defined through the
% following symbolic variables:
%
% G(x,y) = alpha * Pox + beta * Poy, where x,y in {j,_j}
% A(x,y) = alpha - alpha * Pox - beta * Poy = alpha - Gxy, where x,y in {j,_j}.
% B(x,y) = beta - beta * Pox - alpha * Poy = beta - Gyx, where x,y in {j,_j}.
%
% We could use only one index in the following way:
% G(j,-j) = alpha * Poj + beta * Po_j = Gj.
% G(-j,j) = alpha * Po_j + beta * Poj = G_j.
% A(j,-j) = alpha - G(j,-j) = Aj = alpha - Gj.
% A(-j,j) = alpha - G(-j,j) = A_j = alpha - G_j.
% B(j,-j) = beta - beta * Poj - alpha * Po_j = Bj = beta - G_j.
% B(-j,j) = beta - beta * Po_j - alpha * Poj = B_j = beta - Gj.

% Define these additional symbolic variables.
syms Aj Bj Gj
syms A_j B_j G_j
syms Tj T_j

% Define the alternative expressions for the Qj and Q_j.
Qj =  Aj*Sj + B_j*S_j + Gj - gamma_prime;
Q_j = A_j*S_j + Bj*Sj + G_j - gamma_prime;

% Define the expressions for (S)So, (S)Sj and (S)S_j.
SSo = Lj*L_j + L_j * Tj + Lj * T_j;
SSj = L_j * Tj / SSo;
SS_j = Lj * T_j / SSo;

% Uncomment the following lines of code in order to rewrite quantities Aj,
% A_j, Bj and B_j as functions of the expressions Gj and G_j.
% Qj = subs(Qj,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
% Q_j = subs(Q_j,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);

% Uncomment the following lines of code in order to set gamma_prime = 0.
Qj = subs(Qj,gamma_prime,0);
Q_j = subs(Q_j,gamma_prime,0);

% Form the expressions for the revenues Rj and R_j of the two firms.
Rj = Qj^2;
R_j = Q_j^2;

% Form the expressions for the cost functions Cj and C_j of the two firms.
% Cj = G*Sj*(1-S_j);
% C_j = G*S_j*(1-Sj);
Cj = G*Sj;
C_j = G*S_j;

%**************************************************************************
%                   MAIN TRICK!!!
%**************************************************************************
% Define all the derivatives of Sj and S_j with respect to both Tj and T_j.
DSjTj = (L_j/So) * (1 - Sj);
DS_jTj = -(L_j/So) * S_j;
DSjT_j = -(Lj/So) * Sj;
DS_jT_j = (Lj/So) * (1 - S_j);

% Define the first derivatives of Qj and Q_j with respect to Tj and T_j
% accordingly.
DQjTj = Aj * DSjTj + B_j * DS_jTj;
DQ_jT_j = A_j * DS_jT_j + Bj * DSjT_j;

% Define the first derivatives of Rj and R_j with respect to Tj and T_j.
DRjTj = 2 * Qj * DQjTj;
DR_jT_j = 2 * Q_j * DQ_jT_j;
% Express the first derivatives of revenues as fractions.
[N_DRjTj,D_DRjTj] = numden(DRjTj);
[N_DR_jT_j,D_DR_jT_j] = numden(DR_jT_j);
% Simplify the numerators.
N_DRjTj = simplify(collect(expand(N_DRjTj)));
N_DR_jT_j = simplify(collect(expand(N_DR_jT_j)));
% Express quantities N_DRjTj and N_DR_jT_j as polynomials of Sj and S_j.
[C_N_DRjTj,T_N_DRjTj] = coeffs(N_DRjTj,[Sj S_j]);
[C_N_DR_jT_j,T_N_DR_jT_j] = coeffs(N_DR_jT_j,[Sj S_j]);

% Define the first derivatives of Cj and C_j with respect to Tj and T_j.
% DCjTj = G*(DSjTj*(1-S_j) - Sj*DS_jTj);
% DC_jT_j = G*(DS_jT_j*(1-Sj) - S_j*DSjT_j);
DCjTj = G*DSjTj;
DC_jT_j = G*DS_jT_j;

% Define the first derivatives of Fj and F_j with respect to Tj and T_j.
DFjTj = DRjTj - DCjTj;
DF_jT_j = DR_jT_j - DC_jT_j;

% Express the previous quantities as fractions.
[Nj,Dj] = numden(DFjTj);
[N_j,D_j] = numden(DF_jT_j);

% Substitute expressions for Aj,A_j and Bj and B_j as a function of the
% expressions for Gj and G_j in Nj and N_j.
Nj = subs(Nj,[Aj,A_j,Bj,Bj],[alpha - Gj,alpha - G_j,beta - G_j,beta - Gj]);
N_j = subs(N_j,[Aj,A_j,Bj,Bj],[alpha - Gj,alpha - G_j,beta - G_j,beta - Gj]);

% Substitute expressions for Gj and G_j within the expressions for Nj and
% N_j.
Nj = subs(Nj,[Gj,G_j],[alpha * Poj + beta * Po_j,alpha * Po_j + beta * Poj]);
N_j = subs(N_j,[Gj,G_j],[alpha * Poj + beta * Po_j,alpha * Po_j + beta * Poj]);