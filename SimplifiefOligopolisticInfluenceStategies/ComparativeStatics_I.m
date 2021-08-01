% This script file provides the required functionality in order to perform
% comparative statics analysis on the aquired optimal solutions for the
% game that underlies the Simplified Oligopolistic OPtimal Influence Model.

% We need to determine the behaviour of the Nash Equilibrium point 
% p = (TA,TB) as a function of the exogenous parameters of the model such
% as: PA, PB, LA, LB, K, M, C, and G. Thus, if we let 
% w in {PA, PB, LA, LB, K, M, C,G}, we need to determine the behaviour of 
% p(w) = (TA(w),TB(w)).

clc
clear 
% Declare the necessary symbolic variables.
syms C G LA LB PA PB TA TB alpha beta gamma

% Declare the additional symbolic variables dTA and dTB.
syms dTA dTB

% Set the expressions for the payoff functions fa and fb with respect to
% the above parameters.
fa = (C + gamma)^2 - G*TA^2 + (LB^2*alpha^2*(TA + LA*PA)^2)/(LA*LB + LA*TB ...
     + LB*TA)^2 + (LA^2*beta^2*(TB + LB*PB)^2)/(LA*LB + LA*TB + LB*TA)^2 - ...
     (2*LB*alpha*(C + gamma)*(TA + LA*PA))/(LA*LB + LA*TB + LB*TA) - ...
     (2*LA*beta*(C + gamma)*(TB + LB*PB))/(LA*LB + LA*TB + LB*TA) + ...
     (2*LA*LB*alpha*beta*(TA + LA*PA)*(TB + LB*PB))/(LA*LB + LA*TB + LB*TA)^2;
 
fb = (C + gamma)^2 - G*TB^2 + (LA^2*alpha^2*(TB + LB*PB)^2)/(LA*LB + LA*TB ...
     + LB*TA)^2 + (LB^2*beta^2*(TA + LA*PA)^2)/(LA*LB + LA*TB + LB*TA)^2 - ...
     (2*LA*alpha*(C + gamma)*(TB + LB*PB))/(LA*LB + LA*TB + LB*TA) - ...
     (2*LB*beta*(C + gamma)*(TA + LA*PA))/(LA*LB + LA*TB + LB*TA) + ...
     (2*LA*LB*alpha*beta*(TA + LA*PA)*(TB + LB*PB))/(LA*LB + LA*TB + LB*TA)^2;

% The following line of code ensures the positiveness of the qj and q_j
% quantities.
% fa = subs(fa,[C gamma],[0 0]);
% fb = subs(fb,[C gamma],[0 0]);
 
% Compute the first order derivatives with respect to TA and TB
% accordingly.
Da = diff(fa,TA); %Da = Dfa_TA
Db = diff(fb,TB); %Db = Dfb_TB

% If DTA_w (--> dTA) is the total derivative of TA with respect to parameter 
% w and DTB_w (--> dTB)is the total derivative of TB with respect to parameter w then 
% the following equations hold:
%                               Ax * dTA + Bx * dTB = - Gx [1]
%                               Ay * dTA + By * dTB = - Gy [2]
% where:
%        Ax = DDa_TA = Dfa_TA_TA [3]
%        Bx = DDa_TB = Dfa_TA_TB [4]
%        Gx = DDa_w  = Dfa_TA_w  [5]
%        Ay = DDb_TA = Dfb_TA_TB [6]
%        By = DDb_TB = Dfb_TB_TB [7]
%        Gy = DDb_w  = Dfb_TB_w  [8]

% Compute the quantities Ax and Bx as well as the quantities Ay and By.
Ax = diff(Da,TA);
Bx = diff(Da,TB);
Ay = diff(Db,TA);
By = diff(Db,TB);

% Compute the quantities Gx and Gy where the parameter w is PA. 
Gx = diff(Da,PA);
Gy = diff(Db,PA);

% Re-express the previously defined quantities as fractions.
[NAx,DAx] = numden(Ax);
[NBx,DBx] = numden(Bx);
[NGx,DGx] = numden(Gx);

[NAy,DAy] = numden(Ay);
[NBy,DBy] = numden(By);
[NGy,DGy] = numden(Gy);

% Define the system of linear equations with respect to dTA and dTB.
S_PA = [Ax * dTA + Bx * dTB + Gx==0;Ay * dTA + By * dTB + Gy==0];
% Solve the linear system with respect to dTA and dTB.
S_PA_solution = solve(S_PA,[dTA dTB]);

% Provide an alternative solution by solving the system of linear equations
% with respect to the auxiliary variables Ax, Bx, Gx, Ay, By and Gy.

% Declare the auxiliary symbolic variables.
syms AAx BBx GGx AAy BBy GGy

R = [AAx BBx;AAy BBy];
Q = [-GGx;-GGy];
So = R\Q;

% Keep in mind that the exact value of So is given by the following
% equation:
%        |(BBx*GGy - BBy*GGx)/(AAx*BBy - AAy*BBx) |   |dTA|   |Wa / Ra|
% So  =  |                                        | = |   | = |       | [9]
%        |-(AAx*GGy - AAy*GGx)/(AAx*BBy - AAy*BBx)|   |dTB|   |Wb / Rb|
%
% According to the previous definitions we may write that:
% Wa = Bx * Gy - By * Gx = BGxy - BGyx [10]
% Wb = -(Ax * Gy - Ay * Gx) = -(AGxy - AGyx) [11]
% Ra = Rb = R = Ax * By  - Ay * Bx = ABxy - AByx [12]

% Compute the previously defined intermediate quantities:
ABxy = Ax * By;
AByx = Ay * Bx;
BGxy = Bx * Gy;
BGyx = By * Gx;
AGxy = Ax * Gy;
AGyx = Ay * Gx;
% Decompose the above quantities into numerator - denumerator fractions.
[NABxy,DABxy] = numden(ABxy);
[NAByx,DAByx] = numden(AByx);
[NBGxy,DBGxy] = numden(BGxy);
[NBGyx,DBGyx] = numden(BGyx);
[NAGxy,DAGxy] = numden(AGxy);
[NAGyx,DAGyx] = numden(AGyx);

% Solve the systems of differential equations.
Soo = subs(So,[AAx BBx GGx AAy BBy GGy],[Ax Bx Gx Ay By Gy]);
Soo = simplify(Soo);
syms x(PA) y(PA)
Soo = subs(Soo,[TA TB],[x y]);
eqn1 = diff(x,PA) == Soo(1);
eqn2 = diff(y,PA) == Soo(2);
eqns = [eqn1;eqn2];
% Doo = dsolve(eqns);