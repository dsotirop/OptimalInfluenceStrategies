% This script file performs auxiliary symbolic calculations in order to
% obtain an analytical solution for the continous game underlying the
% Simplified Oligopolistic Optimal Influence model.

clc
clear

% Define fundamental symbolic variables.
syms alpha beta gamma_prime G
syms Poj Po_j
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

% Define the alternative expressions for the Qj and Q_j.
Qj =  Aj*Sj + B_j*S_j + Gj - gamma_prime;
Q_j = A_j*S_j + Bj*Sj + G_j - gamma_prime;

% Uncomment the following lines of code in order to rewrite quantities Aj,
% A_j, Bj and B_j as functions of the expressions Gj and G_j.
% Qj = subs(Qj,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
% Q_j = subs(Q_j,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);

% Uncomment the following lines of code in order to set gamma_prime = 0.
Qj = subs(Qj,gamma_prime,0);
Q_j = subs(Q_j,gamma_prime,0);

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

% Given that the profit functions for the two firms are defined as:
% Fj = Qj^2 - G * Sj^2
% F_j = Q_j^2 - G * S_j^2
% we may compute their first derivatives with respect to Tj and T_j as:
% DFjTj = 2 * Qj * DQjTj - 2 * G * Sj * DSjTj
% DF_jT_j = 2 * Q_j * DQ_jT_j - 2 * G * DS_jT_j

% Define the first derivatives of Fj and F_j with respect to Tj and T_j
% accordingly.
DFjTj = 2 * Qj * DQjTj - 2 * G * Sj * DSjTj;
DF_jT_j = 2 * Q_j * DQ_jT_j - 2 * G * DS_jT_j;

% Express quantities DFjTj and DF_jT_j as fractions.
[DFjTj_n,DFjTj_d] = numden(DFjTj);
[DF_jT_j_n,DF_jT_j_d] = numden(DF_jT_j);

% Re-express the numerator of the fractions DFjTj and DF_jT_j as
% multivariate polynomials of both Sj and S_j.
[DFjTj_n_c,DFjTj_n_t] = coeffs(DFjTj_n,[Sj,S_j]);
[DF_jT_j_n_c,DF_jT_j_n_t] = coeffs(DF_jT_j_n,[Sj,S_j]);

% Re-express the numerators of the fractions DFjTj and DF_jT_j as
% polynomials of Sj and S_j respectively.
[CNj,TNj] = coeffs(DFjTj_n,Sj);
[CN_j,TN_j] = coeffs(DF_jT_j_n,S_j);

% Compute the best response functions Sj_star(S_j) and S_j_star(Sj) as the
% roots of the polynomials DFjTj_n and DF_jT_j_n respectively:
% In other words, solve equations DFjTj_n(Sj,S_j) = 0 and DF_jT_j_n(Sj,S_j) = 0 with respect 
% to Sj and S_j respectively.
Sj_star = roots(CNj);
S_j_star = roots(CN_j);

% The Nash Equilibrium solution points with respect to Sj and S_j are the
% following:
% Sj_opt_1, Sj_opt2, Sj_opt_3, Sj_opt_4 and
% S_j_opt_1, S_j_opt2, S_j_opt_3, S_j_opt_4

% For the first pair of NE points (Sj_opt_1,S_j_opt_1):
Sj_star_1 = subs(Sj_star(1),S_j,S_j_star(1));
Sj_eql_1 = Sj_star_1 - Sj;

% For the second pair of NE points (Sj_opt_2,S_j_opt_2):
Sj_star_2 = subs(Sj_star(1),S_j,S_j_star(2));
Sj_eql_2 = Sj_star_2 - Sj;

% For the third pair of NE points (Sj_opt_3,S_j_opt_3):
Sj_star_3 = subs(Sj_star(2),S_j,S_j_star(1));
Sj_eql_3 = Sj_star_3 - Sj;

% For the fourth pair of NE points (Sj_opt_4,S_j_opt_4):
Sj_star_4 = subs(Sj_star(2),S_j,S_j_star(2));
Sj_eql_4 = Sj_star_4 - Sj;

% ************************************************************************
%           ALTERNATIVE COMPUTATION OF THE EQUILIBRIUM POINTS
% ************************************************************************
% Quantities DFjTj_n and DF_jT_j are multivariate polynomials with respect 
% to Sj = x and S_j = y that can be expressed in the following form:
% P(x,y) = DFjTj_n = Axx * x^2 + Ayy * y^2 + Axy * x * y + Ax * x + Ay * y + Ao   [I]
% Q(x,y) = DF_jT_j_n = Bxx * x^2 + Byy * y^2 + Bxy * x * y + Bx * x + By * y + Bo [II]
%
% Given the polynomial representation of the quantities DFjTj_n and DF_jT_j_n
% we may write that:
%                   DFjTj_n = <Cx,Tx>  = <DFjTj_n_c,DFjTj_n_t>      [III]
%                   DF_jT_j_n = <Cy,Ty> = <DF_jT_j_n_c,DF_jT_j_n_t> [IV]
%
% We know that Tx = Ty = [x^2 x*y x y^2 y 1] [V]
% Therefore, the correct set of polynomial coefficients may be derived as:
% Axx = Cx(1) Axy = Cx(2) Ax = Cx(3) Ayy = Cx(4) Ay = Cx(5) Ao = Cx(6)
% Bxx = Cy(1) Bxy = Cy(2) Bx = Cy(3) Byy = Cy(4) By = Cy(5) Bo = Cy(6)

% Get the above quantities:
Axx = simplify(collect(expand(DFjTj_n_c(1))));
Axy = simplify(collect(expand(DFjTj_n_c(2))));
Ax = simplify(collect(expand(DFjTj_n_c(3))));
Ayy = simplify(collect(expand(DFjTj_n_c(4))));
Ay = simplify(collect(expand(DFjTj_n_c(5))));
Ao = simplify(collect(expand(DFjTj_n_c(6))));

Bxx = simplify(collect(expand(DF_jT_j_n_c(1))));
Bxy = simplify(collect(expand(DF_jT_j_n_c(2))));
Bx = simplify(collect(expand(DF_jT_j_n_c(3))));
Byy = simplify(collect(expand(DF_jT_j_n_c(4))));
By = simplify(collect(expand(DF_jT_j_n_c(5))));
Bo = simplify(collect(expand(DF_jT_j_n_c(6))));

% Taking into consideration the fact that the above coefficients define the
% following bivariate polynomials P(x,y) and Q(x,y), we may write that:
% P(x,y) = Py(x) = Axx * x^2 + (Axy * y + Ax) * x + (Ayy * y^2 + Ay * y + Ao)
% P(x,y) = Px(y) = Ayy * y^2 + (Axy * x + Ay) * y + (Axx * x^2 + Ax * x + Ao)
% Q(x,y) = Qy(x) = Bxx * x^2 + (Bxy * y + Bx) * x + (Byy * y^2 + By * y + Bo)
% Q(x,y) = Qx(y) = Byy * y^2 + (Bxy * x + By) * y + (Bxx * x^2 + Bx * x + Bo)

% One line of thought towards acquring the solution points to the above
% system of quadratic equations is to eliminate y-related coefficients from
% Py(x) so that it reduces to a polynomial P(x) and to eliminate x-related 
% coefficients from Qx(y) so that it reduces to Q(y). The y-related
% coefficients of Py(x) are {Ayy,Axy,Ay}, whereas the x-related
% coefficients of Qx(y) are {Bxx,Bxy,Bx}. Specifically, the previously
% discussed coefficients have the following values:


% Define the auxiliary symbolic variables x any y in order to form the
% above expressions.
syms x y 
Po = sum([Axx Ayy Axy Ax Ay Ao] .* [x^2 y^2 x*y x y 1]);
Qo = sum([Bxx Byy Bxy Bx By Bo] .* [x^2 y^2 x*y x y 1]);

% Get the coefficients and terms for the polynomials Po and Qo with respect
% to x and y.
[CPox,TPox] = coeffs(Po,x);
[CPoy,TPoy] = coeffs(Po,y);
[CQox,TQox] = coeffs(Qo,x);
[CQoy,TQoy] = coeffs(Qo,y);
% Compute the discriminants of the polynomials Po and Qo with respect to x
% and y.
DPox = CPox(2)^2 - 4 * CPox(1) * CPox(3);
DPoy = CPoy(2)^2 - 4 * CPoy(1) * CPoy(3);
DQox = CQox(2)^2 - 4 * CQox(1) * CQox(3);
DQoy = CQoy(2)^2 - 4 * CQoy(1) * CQoy(3);
% Simplify the expressions above.
DPox = simplify(collect(expand(DPox)));
DPoy = simplify(collect(expand(DPoy)));
DQox = simplify(collect(expand(DQox)));
DQoy = simplify(collect(expand(DQoy)));
% Substitute quantities Aj, A_j and Bj, B_j into the above expressions.
DPx = subs(DPox,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
DPy = subs(DPoy,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
DQx = subs(DQox,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
DQy = subs(DQoy,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
% Substitute quantities Gj and G_j in the above expressions.
DPx = subs(DPx,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
DPy = subs(DPy,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
DQx = subs(DQx,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
DQy = subs(DQy,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
% Simplify expressions for the discriminants DPx, DPy, DQx and DQy.
DPx = simplify(collect(expand(DPx)));
DPy = simplify(collect(expand(DPy)));
DQx = simplify(collect(expand(DQx)));
DQy = simplify(collect(expand(DQy)));
% Substitute expressions Aj, A_j, Bj, B_j, Gj and G_j in the corresponding
% polynomials Po and Qo.
P = subs(Po,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
P = subs(P,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
Q = subs(Qo,[Aj,A_j,Bj,B_j],[alpha-Gj,alpha-G_j,beta-G_j,beta-G_j]);
Q = subs(Q,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
% Get the coefficients and terms for the polynomials P and Q with respect
% to x and y.
[CPx,TPx] = coeffs(P,x);
[CPy,TPy] = coeffs(P,y);
[CQx,TQx] = coeffs(Q,x);
[CQy,TQy] = coeffs(Q,y);

% Uncomment the following lines of code in order to impose B_j = Bj = 0.
CPox = transpose(subs(CPox,[Bj B_j],[0 0]));
CPoy = transpose(subs(CPoy,[Bj B_j],[0 0]));
CQox = transpose(subs(CQox,[Bj B_j],[0 0]));
CQoy = transpose(subs(CQoy,[Bj B_j],[0 0]));
DPox = subs(DPox,[Bj B_j],[0 0]);
DPoy = subs(DPoy,[Bj B_j],[0 0]); % Apparently this is equal to 0 since all coefficients associated with variable y have been eliminated.
DQox = subs(DQox,[Bj B_j],[0 0]); % Apparently this is equal to 0 since all coefficients associated with variable x have been eliminated.
DQoy = subs(DQoy,[Bj B_j],[0 0]);
% Once again simplify expressions for the remaining discriminants DPox and
% DQoy.
DPox = simplify(collect(expand(DPox)));
DQoy = simplify(collect(expand(DQoy)));

% Since MatLab does not simplify sqrt(x^2) into x, it will be done manualy
% for the expressions DPox and DPQoy.
DPoy_copy = 2*L_j*(G + Aj*gamma_prime - Aj^2 - Aj*Gj);
DQoy_copy = 2*Lj*(G + A_j*gamma_prime - A_j^2 - A_j*G_j);

% Get the solutions with respect to x.
Sox_pos = (-CPox(2) + DPoy_copy) / (2*CPox(1));
Sox_neg = (-CPox(2) - DPoy_copy) / (2*CPox(1));

% Simplify expressions for Sox_pos and Sox_neg.
Sox_pos = simplify(collect(expand(Sox_pos))); % Sox_pos = 1.
Sox_neg = simplify(collect(expand(Sox_neg)));

% Get the solutions with respect to y.
Soy_pos = (-CQoy(2) + DQoy_copy) / (2*CQoy(1));
Soy_neg = (-CQoy(2) - DQoy_copy) / (2*CQoy(1));

% Simplify expressions for Soy_pos and Soy_neg.
Soy_pos = simplify(collect(expand(Soy_pos))); % Soy_pos = 1.
Soy_neg = simplify(collect(expand(Soy_neg)));

% Since Sox_pos ==> Sj and Soy_pos ==> S_j and given that 
% Sj + S_j + Sc = 1, we have to reject the pair of solutions Sox_pos == 1 
% and Soy_pos == 1.

% Substitute quantities Aj, A_j, Gj and G_j into the previously
% defined solutions Sox_neg and Soy_neg.
Sox_neg = subs(Sox_neg,[Aj,A_j],[alpha-Gj,alpha-G_j]);
Sox_neg = subs(Sox_neg,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);
Soy_neg = subs(Soy_neg,[Aj,A_j],[alpha-Gj,alpha-G_j]);
Soy_neg = subs(Soy_neg,[Gj G_j],[alpha * Poj + beta * Po_j alpha * Po_j + beta * Poj]);