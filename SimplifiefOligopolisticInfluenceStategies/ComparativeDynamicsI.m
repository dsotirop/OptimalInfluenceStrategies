% This script file provides fundamental computational functionality for
% performing comparative statics on the Simplified Oligopolistic Influence
% Model from a Dynamical Systems Perspective. The underlying mathematical
% background is given by Michael R. Caputo in the following paper:
% "The Envelope Theorem and Comparative Statics of Nash Equilibria"

% Clear command window and workspace.
clc
clear

% Define the required symbolic variables.
syms Xx Xy M K Px Py G Go
% Define additional symbolic variables.
syms px py Sx Sy
% Define additional symbolic variables.
syms Lx Ly Tx Ty So

% Define the expressions for the quantities.
Qx = (Xx - K*Xy) - (Px - M*Py);
Qy = (Xy - K*Xx) - (Py - M*Px);

% Define the expressions for the revenues.
Rx = Px * Qx;
Ry = Py * Qy;

% Define the expressions for the cost functions.
Cx = G * Tx^2;
Cy = G * Ty^2;

% Define the expressions for the profits.
Fx = Rx - Cx;
Fy = Ry - Cy;

% Derive the expressions for the first order derivatives of the profits
% with respect to Px and Py.
DFx_Px = diff(Fx,Px);
DFy_Py = diff(Fy,Py);

% Solve the equation DFx_Px == 0 with respect to Px in order to acquire the
% Best Repsponse Function for firm x which determines the corresponding
% value Px_star as a function of Py. Therefore, we have that:
% Px = Px_star(Py).
Px_star = solve(DFx_Px==0,Px);

% Solve the equation DFy_Py == 0 with respect to Py in order to acquire the
% Best Repsponse Function for firm y which determines the corresponding
% value Py_star as a function of Px. Therefore, we have that:
% Py = Py_star(Px).
Py_star = solve(DFy_Py==0,Py);

% At this point we need to acquire the symbolic expressions for the optimal
% values of Px and Py (namely, Px_opt and Py_opt) by performing the following
% substitutions: 
% (i)  Px_opt may be obtained by solving with respect to Px the equation
%      Px_opt = Px(Py(Px_opt)).
% (ii) Py_opt may be obtained by solving with respect to Py the equation
%      Py_opt = Py(Px(Px_opt)).

% Make copies of the original variables in order to perform the necessary
% substitutions.
Px_star_copy = Px_star;
Py_star_copy = Py_star;

% Express variable Py in the Px_star(Py) expression of Px_star as Py_star.
Px_star = subs(Px_star,Py,Py_star_copy);
% Solve the equation Px = Px_star(Py_star(Px) with respect to Px and assign 
% the corresponding value to Px_opt.
Px_opt = solve(Px_star==Px,Px);

% Express variable Px in the Py_star(Px) expression of Py_star as Px_star.
Py_star = subs(Py_star,Px,Px_star_copy);
% Solve the equation Py = Py_star(Px_star(Py) with respect to Py and assign 
% the corresponding value to Py_opt.
Py_opt = solve(Py_star==Py,Py);

% Simplify the expressions for the Px_opt and Py_opt.
Px_opt = simplify(collect(expand(Px_opt)));
Py_opt = simplify(collect(expand(Py_opt)));

% Substitute the optimal price expressions for Pxo and Pyo within the 
% expressions for the profits Fx and Fy.
Fx = subs(Fx,[Px Py],[Px_opt Py_opt]);
Fy = subs(Fy,[Px Py],[Px_opt Py_opt]);

% In this context, we may express quantities Xx and Xy as functions of the
% initial beliefs px and py and the corresponding limiting influences
% according to the following equations:
%
%                  Xx = px + (1 - px)*Sx - px*Sy [1]
%                  Xy = py + (1 - py)*Sy - py*Sx [2]


% Substitute quantities Xx and Xy defined in equations [1] and [2] within the
% expressions for Fxo and Fyo.
Fxo = subs(Fx,[Xx Xy],[px + (1 - px)*Sx - px*Sy,py + (1 - py)*Sy - py*Sx]);
Fyo = subs(Fy,[Xx Xy],[px + (1 - px)*Sx - px*Sy,py + (1 - py)*Sy - py*Sx]);

% At this point we need to define the limiting influece for each firm
% accroding to the following equations:
%
%             Ly*Tx                                Lx*Ty
% Sx = ---------------------- [3] and Sy = ---------------------- [4]
%      Ly*Tx + Lx*Ty + Lx*Ly                Ly*Tx + Lx*Ty + Lx*Ly
%
% Let So = Ly*Tx + Lx*Ty + Lx*Ly [5]
%

% Substitute quantities Sx and Sy within the refined expressions for Fxo 
% and Fyo.
FX = subs(Fxo,[Sx Sy],[(Ly*Tx)/So,(Lx*Ty)/So]);
FY = subs(Fyo,[Sx Sy],[(Ly*Tx)/So,(Lx*Ty)/So]);

% Substitute quantity So into the expressions for FX and FY.
FX = subs(FX,So,Ly*Tx + Lx*Ty + Lx*Ly);
FY = subs(FY,So,Ly*Tx + Lx*Ty + Lx*Ly);

% Compute first order derivatives with respect to Tx and Ty.
DFxTx = diff(FX,Tx);
DFyTy = diff(FY,Ty);

% Compute second order derivatives.
D2FxTxTx = diff(DFxTx,Tx);
D2FxTxTy = diff(DFxTx,Ty);
D2FyTyTy = diff(DFyTy,Ty);
D2FyTxTy = diff(DFyTy,Tx);

% Define additional second order derivatives.
D2FxTxpx = diff(DFxTx,px);
D2FyTypx = diff(DFyTy,px);

% Set the following symbolic variables:
%            1                      1
% Ax = - ---------   and Ay = - ---------
%        D2FxTxTx                D2FyTyTy
%
% Bx = D2FxTxTy      and By = D2FyTxTy
%
% Gx = D2FxTxpx      and Gy = D2FyTypx

% Define the following system of equations with respect to x and y such
% that:
% .     dTx                  .     dTy
% x = -------  = fx(x,y) and y = ------- = fy(x,y)
%       dpx                        dpx
%
% |1     -Ax*Bx| |x|   | Ax*Gx|
% |            | | | = |      |
% |Ay*By   -1  | |y|   |-Ay*Gy|

% Define the necessary symbols.
syms Ax Ay Bx By Gx Gy
S = [1 -Ax*Bx;Ay*By -1];
C = [Ax*Gx;-Ay*Gy];
Z = S\C;
x = Z(1);
y = Z(2);

% Set the quantities Axo, Ayo, Bxo, Byo, Gxo and Gyo.
Axo = -(1/D2FxTxTx);
Ayo = -(1/D2FyTyTy);
Bxo = D2FxTxTy;
Byo = D2FyTxTy;
Gxo = D2FxTxpx;
Gyo = D2FyTypx;
% Substitute the previous quantities within the expresssions for x and y.
fx = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo Gyo]);
fy = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo Gyo]);
% Express the previous quantities as fractions.
[nx,dx] = numden(fx);
[ny,dy] = numden(fy);
% Simplify the numerators.
nx = simplify(collect(expand(nx)));
ny = simplify(collect(expand(ny)));

% Set values for the external parametets.
Lxo = 0.50;
Lyo = 0.50;
pxo = 0.50;
pyo = 0.50;
Mo = 0.4;
Ko = 0.1;
alpha = (Ko*Mo - 2) / (Mo^2 - 4);
beta = (2*Ko - Mo) / (Mo^2 - 4);
Fmax = max(alpha^2,(beta^2+alpha*beta));
Lmin = min(Lxo,Lyo);
Go = Fmax  / (Lmin^2);

% Substitute the above values into the expressions for nx and ny.
nx = subs(nx,[G Lx Ly px py K M],[Go Lxo Lyo pxo pyo Ko Mo]);
ny = subs(ny,[G Lx Ly px py K M],[Go Lxo Lyo pxo pyo Ko Mo]);

% Set the system for deriving the fixed points of the system.
Sys = [nx==0;ny==0];
Sol = solve(Sys,[Tx Ty]);
Sol.Tx = vpa(Sol.Tx)
Sol.Ty = vpa(Sol.Ty)

% Set the right hand sides for the system of non-linear differential
% equations.
DTxpx = fx;
DTypx = fy;
% Replace parameters that remain constant with the corresponding values.
DTxpx = subs(DTxpx,[G Lx Ly py K M],[Go Lxo Lyo pyo Ko Mo]);
DTypx = subs(DTypx,[G Lx Ly py K M],[Go Lxo Lyo pyo Ko Mo]);
% Define the main variables for the system of non-linear differential
% equations.
% syms tx(px) ty(px)
% Use the newly defined symbolic variables to define the final expresssions 
% for DTxpx and DTypx.
% DTxpx = subs(DTxpx,[Tx Ty],[tx ty]);
% DTypx = subs(DTypx,[Tx Ty],[ty ty]);
% Set the equations that define the system.
% plot(t,y,'-o')