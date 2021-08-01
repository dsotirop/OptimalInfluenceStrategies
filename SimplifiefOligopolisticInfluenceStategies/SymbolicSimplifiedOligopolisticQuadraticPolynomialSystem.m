% This script file focuses on extracting the solutions to the system of
% polynomial equations which are derived by applying the KKT conditions on
% the optimization problems associated with the continuous game that underlies
% the Simplified Oligopolistic Optimal Influence Model.

% Given that the first derivatives of the objective functions fa and fb, Da
% and Db are expressed in the following way:
% 
%                   Ua             Ub  
%             Da = ----  and Db = ---- .  
%                   Va             Vb
%
% the derivation of the equilibrium may be reduced to acquring the
% solutions to the following system of bi-variate polynomial equations:
%
%             Ua(TA,TB) = 0
%             Ub(TA,TB) = 0.
%
% It is important to mention that the above system of polynomial equations can
% lead to the Nash Equilibrium points only for the case of internal solutions.

% Declare the necessary symbolic variables. 
% (The relevant quantities are defined elsewhere.)
syms C G LA LB PA PB TA TB alpha beta gamma

% Define the symbolic expressions for the quantities Ua and Ub.
Ua = (-2*G*LB^3)*TA^4 + (-6*G*LA*LB^2)*TA^3*TB + (-6*G*LA*LB^3)*TA^3 + ...
     (-6*G*LA^2*LB)*TA^2*TB^2 + (-12*G*LA^2*LB^2)*TA^2*TB + ...
     (-6*G*LA^2*LB^3)*TA^2 + (-2*G*LA^3)*TA*TB^3 + (-6*G*LA^3*LB)*TA*TB^2 + ...
     (2*LA*LB^2*alpha^2 - 6*G*LA^3*LB^2 - 2*LA*LB^2*alpha*beta - ...
     2*LA*LB^2*alpha*gamma + 2*LA*LB^2*beta*gamma - 2*C*LA*LB^2*alpha + ...
     2*C*LA*LB^2*beta)*TA*TB + (2*LA*LB^3*alpha^2 - 2*G*LA^3*LB^3 - ...
     2*LA*LB^3*alpha*gamma - 2*LA*LB^3*PA*alpha^2 - 2*C*LA*LB^3*alpha + ...
     2*C*LA*LB^3*PA*alpha + 2*C*LA*LB^3*PB*beta - 2*LA*LB^3*PB*alpha*beta + ...
     2*LA*LB^3*PA*alpha*gamma + 2*LA*LB^3*PB*beta*gamma)*TA + ...
     (2*LA^2*LB*alpha*beta - 2*LA^2*LB*beta^2 - 2*LA^2*LB*alpha*gamma + ...
     2*LA^2*LB*beta*gamma - 2*C*LA^2*LB*alpha + 2*C*LA^2*LB*beta)*TB^2 +...
     (2*C*LA^2*LB^2*beta - 4*C*LA^2*LB^2*alpha + 2*LA^2*LB^2*alpha*beta - ...
     4*LA^2*LB^2*alpha*gamma + 2*LA^2*LB^2*beta*gamma + 2*LA^2*LB^2*PA*alpha^2 - ...
     4*LA^2*LB^2*PB*beta^2 + 2*C*LA^2*LB^2*PA*alpha + 2*C*LA^2*LB^2*PB*beta - ...
     4*LA^2*LB^2*PA*alpha*beta + 2*LA^2*LB^2*PB*alpha*beta + 2*LA^2*LB^2*PA*alpha*gamma + ...
     2*LA^2*LB^2*PB*beta*gamma)*TB + 2*LA^2*LB^3*PA*alpha^2 - 2*LA^2*LB^3*PB^2*beta^2 - ...
     2*C*LA^2*LB^3*alpha - 2*LA^2*LB^3*alpha*gamma - 2*LA^2*LB^3*PA^2*alpha^2 + ...
     2*C*LA^2*LB^3*PA*alpha + 2*C*LA^2*LB^3*PB*beta + 2*LA^2*LB^3*PB*alpha*beta + ...
     2*LA^2*LB^3*PA*alpha*gamma + 2*LA^2*LB^3*PB*beta*gamma - 4*LA^2*LB^3*PA*PB*alpha*beta;
 
Ub = (-2*G*LB^3)*TA^3*TB + (-6*G*LA*LB^2)*TA^2*TB^2 + (-6*G*LA*LB^3)*TA^2*TB + ...
     (2*LA*LB^2*alpha*beta - 2*LA*LB^2*beta^2 - 2*LA*LB^2*alpha*gamma + ...
     2*LA*LB^2*beta*gamma - 2*C*LA*LB^2*alpha + 2*C*LA*LB^2*beta)*TA^2 + ...
     (-6*G*LA^2*LB)*TA*TB^3 + (-12*G*LA^2*LB^2)*TA*TB^2 + (2*LA^2*LB*alpha^2 - ...
     6*G*LA^2*LB^3 - 2*LA^2*LB*alpha*beta - 2*LA^2*LB*alpha*gamma + 2*LA^2*LB*beta*gamma - ...
     2*C*LA^2*LB*alpha + 2*C*LA^2*LB*beta)*TA*TB + (2*C*LA^2*LB^2*beta - 4*C*LA^2*LB^2*alpha + ...
     2*LA^2*LB^2*alpha*beta - 4*LA^2*LB^2*alpha*gamma + 2*LA^2*LB^2*beta*gamma + ...
     2*LA^2*LB^2*PB*alpha^2 - 4*LA^2*LB^2*PA*beta^2 + 2*C*LA^2*LB^2*PB*alpha + ...
     2*C*LA^2*LB^2*PA*beta + 2*LA^2*LB^2*PA*alpha*beta - 4*LA^2*LB^2*PB*alpha*beta + ...
     2*LA^2*LB^2*PB*alpha*gamma + 2*LA^2*LB^2*PA*beta*gamma)*TA + (-2*G*LA^3)*TB^4 + ...
     (-6*G*LA^3*LB)*TB^3 + (-6*G*LA^3*LB^2)*TB^2 + (2*LA^3*LB*alpha^2 - 2*G*LA^3*LB^3 - ...
     2*LA^3*LB*alpha*gamma - 2*LA^3*LB*PB*alpha^2 - 2*C*LA^3*LB*alpha + 2*C*LA^3*LB*PB*alpha + ...
     2*C*LA^3*LB*PA*beta - 2*LA^3*LB*PA*alpha*beta + 2*LA^3*LB*PB*alpha*gamma + ...
     2*LA^3*LB*PA*beta*gamma)*TB + 2*LA^3*LB^2*PB*alpha^2 - 2*LA^3*LB^2*PA^2*beta^2 - ...
     2*C*LA^3*LB^2*alpha - 2*LA^3*LB^2*alpha*gamma - 2*LA^3*LB^2*PB^2*alpha^2 + ...
     2*C*LA^3*LB^2*PB*alpha + 2*C*LA^3*LB^2*PA*beta + 2*LA^3*LB^2*PA*alpha*beta + ...
     2*LA^3*LB^2*PB*alpha*gamma + 2*LA^3*LB^2*PA*beta*gamma - 4*LA^3*LB^2*PA*PB*alpha*beta;
 

 
% The following sequence of simplifications may be omitted by commenting out
% the corresponding lines of code.
% syms L P
% Ua = subs(Ua,[LA LB PA PB],[L L P P]);
% Ub = subs(Ub,[LA LB PA PB],[L L P P]);

% In order to ensure that quantities Qa and Qb of the simplified
% oligopolistic model are positive, we have to assume that C = gamma = 0.
% Ua = subs(Ua,[C gamma],[0 0]);
% Ub = subs(Ub,[C gamma],[0 0]);
 
 
% Extract coefficients and corresponding multivariate monomial terms for
% the expressions Ua and Ub with respect to the terms TA and TB.
% Ua(TA,TB) ---> <CUa,TUa> with respect to {TA,TB}.
% Ub(TA,TB) ---> <CUb,TUb> with respect to {TA,TB}.
[CUa,TUa] = coeffs(Ua,[TA,TB]); 
[CUb,TUb] = coeffs(Ub,[TA,TB]);

% Extract coefficients and corresponding multivariate monomial terms for
% the expressions Ua and Ub with respect to the term TA.
% Ua(TA,TB) ---> <CUaa,TUaa> with respect to {TA}.
% Ub(TA,TB) ---> <CUba,TUba> with respect to {TA}.
[CUaa,TUaa] = coeffs(Ua,TA);
[CUba,TUba] = coeffs(Ub,TA);

% Extract coefficients and corresponding multivariate monomial terms for
% the expressions Ua and Ub with respect to the term TB.
% Ua(TA,TB) ---> <CUab,TUab> with respect to {TB}.
% Ub(TA,TB) ---> <CUbb,TUbb> with respect to {TB}.
[CUab,TUab] = coeffs(Ua,TB);
[CUbb,TUbb] = coeffs(Ub,TB);

% Define the system of polynomials for which the Groebner basis will be 
% determined.
Pab = [Ua,Ub];
% 
% % Define the fundamental symbolic variables upon which the Groebner basis
% % will be computed.
% Vars = [TA,TB];
% 
% % Compute the Groebner basis.
% G = gbasis(Pab,Vars);
% 
% % Get the number of polynomials.
% N = length(G);
% 
% % Initialize cell array containers for the coefficients and terms.
% C = cell(1,N);
% T = cell(1,N);
% 
% Ca = cell(1,N);
% Ta = cell(1,N);
% 
% Cb = cell(1,N);
% Tb = cell(1,N);
% 
% % Extract the multivariate polynomial coefficients and corrsponding terms 
% % for each compoment of G with respect to :
% % (i):   TA and TB concurrently ---> extracted terms stored in (C,T)
% % (ii):  TA                     ---> extracted terms stored in (Ca,Ta)
% % (iii): TB                     ---> extracted terms stored in (Cb,Tb)
% for t = 1:1:N
%     Gt = G(t);
%     [C{t},T{t}] = coeffs(Gt,[TA,TB]);
%     [Ca{t},Ta{t}] = coeffs(Gt,TA);
%     [Cb{t},Tb{t}] = coeffs(Gt,TB);
% end

% Try to solve the original system.
% cond1 = alpha>0;
% cond2 = beta>0;
% cond3 = (0<=TA) & (TA<=1);
% cond4 = (0<=TB) & (TB<=1);
% cond5 = (0<=PA) & (PA<=1);
% cond6 = (0<=PB) & (PB<=1);
% cond7 = (0<LA) & (LA<1);
% cond8 = (0<LB) & (LB<1);
% cond9 = G>0;
% conds = [cond1,cond2,cond3,cond4,cond5,cond6,cond7,cond8,cond9];
% assume(conds);
% eqns = [Ua==0,Ub==0];
% S = solve(eqns,[TA TB]);