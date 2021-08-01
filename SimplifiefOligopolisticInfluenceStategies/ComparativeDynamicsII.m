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

% Define additional second order derivatives with respect to px.
D2FxTxpx = diff(DFxTx,px);
D2FyTypx = diff(DFyTy,px);
% Define additional second order derivatives with respect to py.
D2FxTxpy = diff(DFxTx,py);
D2FyTypy = diff(DFyTy,py);
% Define additional second order derivatives with respect to K.
D2FxTxK = diff(DFxTx,K);
D2FyTyK = diff(DFyTy,K);
% Define additional second order derivatives with respect to M.
D2FxTxM = diff(DFxTx,M);
D2FyTyM = diff(DFyTy,M);
% Define additional second order derivatives with respect to Lx.
D2FxTxLx = diff(DFxTx,Lx);
D2FyTyLx = diff(DFyTy,Lx);
% Define additional second order derivatives with respect to Ly.
D2FxTxLy = diff(DFxTx,Ly);
D2FyTyLy = diff(DFyTy,Ly);

% Set the following symbolic variables that will eventually form the
% system of non-linear differential equations.
%
% Ax = D2FxTxTx
% Ay = D2FyTyTy
% Bx = D2FxTxTy
% By = D2FyTxTy
% Gx = D2FxTxu where u in {px,py,Lx,Ly,K,M}
% Gy = D2FyTyu where u in {px,py,Lx,Ly,K,M}

% Define the following system of equations with respect to x and y such
% that:
% .     dTx                  .     dTy
% x = -------  = fx(x,y) and y = ------- = fy(x,y)
%       dpx                        dpx
%
% |Ax     Bx| |x|   |-Gx|
% |         | | | = |   |
% |By     Ay| |y|   |-Gy|

% Define the necessary symbols.
syms Ax Ay Bx By Gx Gy
S = [Ax Bx;By Ay];
C = [-Gx;-Gy];
Z = S\C;
x = Z(1);
y = Z(2);

% Set the quantities Axo, Ayo, Bxo, Byo, Gxo and Gyo.
Axo = D2FxTxTx;
Ayo = D2FyTyTy;
Bxo = D2FxTxTy;
Byo = D2FyTxTy;

% Set the right hand sides for the system of non-linear differential
% equations.
Gxo_px = D2FxTxpx;
Gyo_px = D2FyTypx;
% Substitute the previous px-related quantities within the expresssions 
% for x and y.
fx_px = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_px Gyo_px]);
fy_px = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_px Gyo_px]);
% Set the right hand sides for the system of non-linear differential
% equations.
DTxpx = fx_px;
DTypx = fy_px;

% Set the quantities Gxo and Gyo.
Gxo_py = D2FxTxpy;
Gyo_py = D2FyTypy;
% Substitute the previous py-related quantities within the expresssions 
% for x and y.
fx_py = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_py Gyo_py]);
fy_py = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_py Gyo_py]);
% Set the right hand sides for the system of non-linear differential
% equations.
DTxpy = fx_py;
DTypy = fy_py;

% Set the quantities Gxo and Gyo.
Gxo_K = D2FxTxK;
Gyo_K = D2FyTyK;
% Substitute the previous K-related quantities within the expresssions 
% for x and y.
fx_K = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_py Gyo_K]);
fy_K = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_py Gyo_K]);
% Set the right hand sides for the system of non-linear differential
% equations.
DTxK = fx_K;
DTyK = fy_K;

% Set the quantities Gxo and Gyo.
Gxo_M = D2FxTxM;
Gyo_M = D2FyTyM;
% Substitute the previous M-related quantities within the expresssions 
% for x and y.
fx_M = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_M Gyo_M]);
fy_M = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_M Gyo_M]);
% Set the right hand sides for the system of non-linear differential
% equations.
DTxM = fx_M;
DTyM = fy_M;

% Set the quantities Gxo and Gyo.
Gxo_Lx = D2FxTxLx;
Gyo_Lx = D2FyTyLx;
% Substitute the previous Lx-related quantities within the expresssions 
% for x and y.
fx_Lx = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_Lx Gyo_Lx]);
fy_Lx = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_Lx Gyo_Lx]);
% Set the right hand sides for the system of non-linear differential
% equations.
DTxLx = fx_Lx;
DTyLx = fy_Lx;

% Set the quantities Gxo and Gyo.
Gxo_Ly = D2FxTxLy;
Gyo_Ly = D2FyTyLy;
% Substitute the previous Ly-related quantities within the expresssions 
% for x and y.
fx_Ly = subs(x,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_Ly Gyo_Ly]);
fy_Ly = subs(y,[Ax Ay Bx By Gx Gy],[Axo Ayo Bxo Byo Gxo_Ly Gyo_Ly]);
% Set the right hand sides for the system of non-linear differential
% equations.
DTxLy = fx_Ly;
DTyLy = fy_Ly;