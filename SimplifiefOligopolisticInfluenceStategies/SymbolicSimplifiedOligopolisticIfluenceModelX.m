% This script file facilitates the conduction of the fundamental symbolic
% operations that underly the eXtended version of the Simplified 
% Oligopolistic Influence Model. The main diversification of the eXtended 
% model lies upon the utilization of an alternative cost function. The
% formulation described here wil assume the same degree of quality and
% price competition.

% For this particular formulation we have alse assummed that the demand
% function is parameterized by a unique competition factor M since we have 
% set that K = M.

% Assuming that Px and Py are the price levels chosen by firms x and y, the
% correpsponding limiting beliefs will be given as Xx and Xy. Additionally,
% M is the degree of competition on the price level and K is the degree of
% competition on the quality level. In this context, the demands for the
% two products will be given by the following equations:
%
% Qx = (Xx - KXy) - (Px - MPy) [1]
% Qy = (Xy - KXx) - (Py - MPx) [2]

% Imposing the simplification assumption K = M, we may write that:
% 
% Qx = (Xx - Px) - M * (Xy - Py) [3]
% Qy = (Xy - Py) - M * (Xx - Px) [4]

% Letting Dx = Xx - Px [5] and Dy = Xy - Py [6], allows the reformulation
% of equations [3] and [4] in the following way:
%
% Qx = Dx - M*Dy [7]
% Qy = Dy - M*Dx [8]

% Mind that according to equations [5] and [6], it is possible to write
% that:
%       Px = Xx - Dx [9] and Py = Xy - Dy [10]

% In this context, the revenue functions for the two firms will be given
% as:
%
%    Rx = Px * Qx = (Xx - Dx) * (Dx - M*Dy) [11]
%    Ry = Py * Qy = (Xy - Dy) * (Dy - M*Dx) [12]
%
% Equations [11] and [12] hold under the assumption that the marginal cost
% c is 0.

% The main diversificarion of this approach deals with the formulation of
% the cost functions which will be defined on the basis of the quantities
% Dx and Dy as:
%               Cx = G * Dx^2 [13] 
%               Cy = G * Dy^2 [14]

% Clear command window and workspace.
clc
clear

% Define the required symbolic variables.
syms Xx Xy M Dx Dy Px Py G Go

% Define the expressions for the quantities.
Qx = Dx - M*Dy;
Qy = Dy - M*Dx;
% Define the expressions for the revenues.
Rx = (Xx - Dx) * Qx;
Ry = (Xy - Dy) * Qy;
% Utilize alternative expressions for the cost that those defined in
% Equations [13] and [14].
% Cx = G * Xx^2;
% Cy = G * Xy^2;
% Cx = G * Dx;
% Cy = G * Dy;
% Cx = (1/2)*Qx^2 + G*Qx + Go;
% Cy = (1/2)*Qy^2 + G*Qy + Go;
% Cx = G * Qx^2;
% Cy = G * Qy^2;
% Cx = G * Dx^2;
% Cy = G * Dy^2;
% Cx = G*Dx*Dy;
% Cy = G*Dx*Dy;
% Cx = M*Px*Dy;
% Cy = M*Py*Dx;
% Cx = Px*(G-Py);
% Cy = (G-Px)*Py;

% Define the expressions for the profits.
Fx = Rx - Cx;
Fy = Ry - Cy;
% Replace quantities Dx and Dy within the expressions for Fx and Fy.
Rx = subs(Rx,[Dx Dy],[Xx - Px,Xy - Py]);
Ry = subs(Ry,[Dx Dy],[Xx - Px,Xy - Py]);
Cx = subs(Cx,[Dx Dy],[Xx - Px,Xy - Py]);
Cy = subs(Cy,[Dx Dy],[Xx - Px,Xy - Py]);
Fx = subs(Fx,[Dx Dy],[Xx - Px,Xy - Py]);
Fy = subs(Fy,[Dx Dy],[Xx - Px,Xy - Py]);

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

% Reformulate the expressions for Px_opt and Py_opt as fractions so that to
% obtain a better analytical formula for the optimal prices of both firms
% for the first stage of the game.
[NPx_opt,DPx_opt] = numden(Px_opt);
[NPy_opt,DPy_opt] = numden(Py_opt);
% Rewrite the numerators of the previous fractional expressions as 
% polynomials of Xx and Xy.
[CNPx_opt,TNPx_opt] = coeffs(NPx_opt,[Xx Xy]);
[CNPy_opt,TNPy_opt] = coeffs(NPy_opt,[Xy Xx]);

% Obtain the optimal expressions for the demand (quantity) functions Qx and
% Qy. Initially, substitute the equivalent expressions for Dx and Dy.
Qx = subs(Qx,[Dx Dy],[Xx-Px Xy-Py]);
Qy = subs(Qy,[Dx Dy],[Xx-Px Xy-Py]);
% Subsequently, substitute the optimal expressions for Px and Py.
Qx_opt = subs(Qx,[Px Py],[Px_opt Py_opt]);
Qy_opt = subs(Qy,[Px Py],[Px_opt,Py_opt]);

% Simplify the expressions for Qx_opt and Qy_opt.
Qx_opt = simplify(collect(expand(Qx_opt)));
Qy_opt = simplify(collect(expand(Qy_opt)));

% Reformulate the expressions for Qx_opt and Qy_opt as fractions so that to
% obtain a better analytical formula for the optimal prices of both firms
% for the first stage of the game.
[NQx_opt,DQx_opt] = numden(Qx_opt);
[NQy_opt,DQy_opt] = numden(Qy_opt);
% Rewrite the numerators of the previous fractional expressions as 
% polynomials of Xx and Xy.
[CNQx_opt,TNQx_opt] = coeffs(NQx_opt,[Xx Xy]);
[CNQy_opt,TNQy_opt] = coeffs(NQy_opt,[Xx Xy]);


% Define additional symbolic variables.
syms px py Sx Sy

% Define the optimal quantities for Pxo and Pyo as functions of Sx and Sy.
Pxo = subs(Px_opt,[Xx Xy],[px + (1 - px)*Sx - px*Sy,py + (1 - py)*Sy - py*Sx]);
Pyo = subs(Py_opt,[Xx Xy],[px + (1 - px)*Sx - px*Sy,py + (1 - py)*Sy - py*Sx]);

% Substitute the optimal price expressions for Pxo and Pyo within the 
% expressions for the profits Fx and Fy.
Fx = subs(Fx,[Px Py],[Px_opt Py_opt]);
Fy = subs(Fy,[Px Py],[Px_opt Py_opt]);

% Decompose quantities Fx and Fy as polynomials of Xx and Xy.
[CFx,TFx] = coeffs(Fx,[Xx Xy]);
[CFy,TFy] = coeffs(Fy,[Xx Xy]);

% In this context, we may express quantities Xx and Xy as functions of the
% initial beliefs px and py and the corresponding limiting influences
% according to the following equations:
%
%                  Xx = px + (1 - px)*Sx - py*Sy [19]
%                  Xy = py + (1 - py)*Sy - px*Sx [20]


% Substitute quantities Xx and Xy defined in equations [19] and [20] within the
% expressions for Fxo and Fyo.
Fxo = subs(Fx,[Xx Xy],[px + (1 - px)*Sx - py*Sy,py + (1 - py)*Sy - px*Sx]);
Fyo = subs(Fy,[Xx Xy],[px + (1 - px)*Sx - py*Sy,py + (1 - py)*Sy - px*Sx]);

% At this point we need to define the limiting influece for each firm
% accroding to the following equations:
%
%             Ly*Tx                                Lx*Ty
% Sx = ---------------------- [17] and Sy = ---------------------- [18]
%      Ly*Tx + Lx*Ty + Lx*Ly                Ly*Tx + Lx*Ty + Lx*Ly
%
% Let So = Ly*Tx + Lx*Ty + Lx*Ly [19]
%

% Define additional symbolic variables.
syms Lx Ly Tx Ty So
% Substitute quantities Sx and Sy within the refined expressions for Fxo 
% and Fyo.
FX = subs(Fxo,[Sx Sy],[(Ly*Tx)/So,(Lx*Ty)/So]);
FY = subs(Fyo,[Sx Sy],[(Ly*Tx)/So,(Lx*Ty)/So]);

% Substitute quantity So into the expressions for FX and FY.
FX = subs(FX,So,Ly*Tx + Lx*Ty + Lx*Ly);
FY = subs(FY,So,Ly*Tx + Lx*Ty + Lx*Ly);

% Compute the first derivatives of the quantities FX and FY with respect to
% Tx and Ty.
DFXTx = diff(FX,Tx);
DFYTy = diff(FY,Ty);

% Express the previous derivatives as fractions.
[NX,DX] = numden(DFXTx);
[NY,DY] = numden(DFYTy);