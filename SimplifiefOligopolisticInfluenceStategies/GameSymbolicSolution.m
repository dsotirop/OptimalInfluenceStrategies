% This script file is an attempt to provide a symbolic solution to a
% particular continuous game where the decision quantities are x and y and
% the corresponding utility functions are given as:
%
%        1          2   1          2
% Fx =  --- * Ax * u + --- * Bx * v  + Gx * u * v + S(x) [1]
%        2              2
%
%        1          2   1          2
% Fy =  --- * Ay * u + --- * By * v  + Gy * u * v + S(y) [2]
%        2              2
%
% where: 
%        u = u(x,y) = Lu * x + Mu * y [3] 
%        v = v(x,y) = Lv * x + Mv * y [4]
% with:
%                  1                        1
%        S(x) = - --- * x [5] and S(y) = - --- * y [6]
%                  2                        2  

% Clear command window and workspace.
clc
clear

% Declare the fundamental symbols for the first stage of the computation.
syms Fx Fy u v Ax Bx Gx Ay By Gy x y Sx Sy
syms Lu Mu Lv Mv

% Declare fundamental quantities for the first stage of the computation.
Fx = 0.5 * Ax * u^2 + 0.5 * Bx * v^2 + Gx * u * v - 0.5 * Sx;
Fy = 0.5 * Ay * u^2 + 0.5 * By * v^2 + Gy * u * v - 0.5 * Sy;

% Compute first order partial derivatives with respect to u and v for both
% quantities.
DFx_u = diff(Fx,u);
DFx_v = diff(Fx,v);
DFy_u = diff(Fy,u);
DFy_v = diff(Fy,v);

% Since Fx and Fy are composite expressions of x and y such that:
%
% Fx = Fx(u(x,y),v(x,y),S(x)) [7] 
% Fy = Fy(u(x,y),v(x,y),S(y)) [8]
%
% derivation of the overall partial derivatives can be done according to:
% DFx_x = DFx_u * Du_x + DFx_v * Dv_x + DSx_x [9] 
% DFy_y = DFy_u * Du_y + DFy_v * Dv_y + DSy_y [10]


% Thus, we need the following additional symbols:
syms Du_x Du_y Dv_x Dv_y DSx_x DSy_y

% Compute the intermediate quantities for DFx_x and DFy_y:
DFx_x = DFx_u * Du_x + DFx_v * Dv_x + DSx_x;
DFy_y = DFy_u * Du_y + DFy_v * Dv_y + DSy_y;

% Collect expressions DFx_x and DFy_y with respect to u and v.
% We know beforehand that DFx_x and DFy_y are linear expresssions with
% respect to u and v.
[Cx,Tx] = coeffs(DFx_x,[u,v]);
[Cy,TY] = coeffs(DFy_y,[u,v]);

% Therefore, we may write that:
%
% DFx_x = Cx(1) * u + Cx(2) * v + Cx(3) [11]
% DFx_x = Cy(1) * u + Cy(2) * v + Cy(3) [12]
% 
% which in turn may be written as:
%
% DFx_x = (Ax*Du_x + Dv_x*Gx) * u + (Bx*Dv_x + Du_x*Gx) * v + DSx_x [13]
% DFy_y = (Ay*Du_y + Dv_y*Gy) * u + (By*Dv_y + Du_y*Gy) * v + DSy_y [14]
%
% or
%
% DFx_x = Au * u + Bu * v + DSx_x [15]
% DFy_y = Av * u + Bv * v + DSy_y [16]
%
% where: Au = Cx(1) Bu = Cx(2) Cu = Cx(3)
%        Av = Cy(1) Bv = Cy(2) Cv = Cy(3)

% Given the additional symbols, the FOCs give the following system of
% linear equations with respect to u and v:
%
% Au * u + Bu * v = -Cu [17]
% Av * u + Bv * v = -Cv [18]

% Define the additional symbols:
syms Au Bu Cu Av Bv Cv
syms Lu Mu Lv Mv

% Derive the solution pair (u_star,v_star) for the system of equations (17)
% and (18).
System = [Au * u + Bu * v + Cu ==0,Av * u + Bv * v + Cv ==0];
Solution = solve(System,[u,v]);
u_star = Solution.u;
v_star = Solution.v;

% Update expressions for u_star and v_star by replacing quantities Au, Bu,
% Cu, Av, Bv and Cv.
u_star = subs(u_star,[Au Bu Cu Av Bv Cv],[Cx(1) Cx(2) Cx(3) Cy(1) Cy(2) Cy(3)]);
v_star = subs(v_star,[Au Bu Cu Av Bv Cv],[Cx(1) Cx(2) Cx(3) Cy(1) Cy(2) Cy(3)]);

% In order to obtain the final expressions for u_star and v_star we need to
% replace the following quantities as:
% DSx_x = -x [19]
% DSy_y = -y [20]
% Du_x = Lu  [21]
% Du_y = Mu  [22]
% Dv_x = Lv  [23]
% Dv_y = Mv  [24]

u_star = subs(u_star,[DSx_x,DSy_y,Du_x,Du_y,Dv_x,Dv_y],[-x -y Lu Mu Lv Mv]);
v_star = subs(v_star,[DSx_x,DSy_y,Du_x,Du_y,Dv_x,Dv_y],[-x -y Lu Mu Lv Mv]);

% Quantities DFx_x and DFy_y will be rewritten by taking into consideration
% equations (3) and (4).
DFx_x = subs(DFx_x,[u,v,DSx_x],[Lu * x + Mu * y,Lv * x + Mv * y,-x]);
DFy_y = subs(DFy_y,[u,v,DSy_y],[Lu * x + Mu * y,Lv * x + Mv * y,-y]);

% Collect expressions DFx_x and DFy_y with respect to x and y.
% We know beforehand that DFx_x and DFy_y are linear expresssions with
% respect to x and y.
[Cxx,Txx] = coeffs(DFx_x,[x y]);
[Cyy,Tyy] = coeffs(DFy_y,[x y]);