% This script file will explore the derivation of some alternative
% formulations for the Oligopolistic Optimal Influence model.

clc
clear

% Define the fundamental symbolic quantities.
syms a b % These are the alhpa and beta quantities.
syms Lx Ly % These are the LA and LB quantities (or Lj and L_j).
syms Px Py % These are the PA and Pb quantities (or Pjo and Po_j).
syms Tx Ty % These are the TA and Tb quantities (or Tj and T_j).
syms Sx Sy % These are the Sj and S_j quantities.
syms Qx Qy % These are the Qj and Q_j quantities.
syms x y   % These are the Xj and X_j quantities.

% Define additional symbolic variables.
So = Lx*Ly + Ly*Tx + Lx*Ty;

% The following derivations aim at expressing the quantities Tx and Ty as
% functions of the corresponding equilibrium quantities Qx and Qy. That is,
% the ultimate goal is to derivation the following expressions:
%
%                      Tx = Ux(Qx,Qy) [1]
%                      Ty = Uy(Qx,Qy) [2]
%
% It is known that:
%
%            Tx = Fx(Sx,Sy) [3] Ty = Fy(Sx,Sy) [4]
%            Sx = Gx(x,y)   [5] Sy = Gy(x,y)   [6]
%            x = Hx(Qx,Qy)  [7] y = Hy(Qx,Qy)  [8]

% Derive the expressions for Tx and Ty:
Fx = Ly * Tx - Sx * So == 0;
Fy = Lx * Ty - Sy * So == 0;
F = [Fx;Fy];
Fsol = solve(F,[Tx Ty]);
TTx = Fsol.Tx;
TTy = Fsol.Ty;

% Derive the expressions for Sx and Sy.
Gx = x == Px + (1 - Px) * Sx - Px * Sy;
Gy = y == Py + (1 - Py) * Sy - Py * Sx;
G = [Gx;Gy];
Gsol = solve(G,[Sx Sy]);
SSx = Gsol.Sx;
SSy = Gsol.Sy;

% Derive the expressions for x and y.
Hx = Qx == a * x + b * y;
Hy = Qy == b * x + a * y;
H = [Hx;Hy];
Hsol = solve(H,[x y]);
xx = Hsol.x;
yy = Hsol.y;

% Obtain the final expressions for Tx and Ty.
TTx = subs(TTx,[Sx Sy],[SSx SSy]);
TTy = subs(TTy,[Sx Sy],[SSx SSy]);
TTx = subs(TTx,[x y],[xx yy]);
TTy = subs(TTy,[x y],[xx yy]);

% At this point TTx and TTy are expressed as functions of the equilibrium
% demand functions Qx and Qy and the rest of the exogenous parameters of
% the model.

% Re-express quantities TTx and TTy as fractions.
[TTx_n,TTx_d] = numden(TTx);
[TTy_n,TTy_d] = numden(TTy);

% Rewrite the numerators and denominators of the expressions TTx and TTy as
% polynomials of Qx and Qy.
[TTx_n_c,TTx_n_t] = coeffs(TTx_n,[Qx,Qy]);
[TTx_d_c,TTx_d_t] = coeffs(TTx_d,[Qx,Qy]);
[TTy_n_c,TTy_n_t] = coeffs(TTy_n,[Qx,Qy]);
[TTy_d_c,TTy_d_t] = coeffs(TTy_d,[Qx,Qy]);

% At this point TTx and TTy are fractions of first degree polynomials with
% respect to both Qx and Qy.

% Transpose the coefficient matrices.
TTx_n_c = transpose(TTx_n_c);
TTx_d_c = transpose(TTx_d_c);
TTy_n_c = transpose(TTy_n_c);
TTy_d_c = transpose(TTy_d_c);

% Siplify the above expressions.
TTx_n_c = simplify(collect(expand(TTx_n_c)));
TTx_d_c = simplify(collect(expand(TTx_d_c)));
TTy_n_c = simplify(collect(expand(TTy_n_c)));
TTy_d_c = simplify(collect(expand(TTy_d_c)));