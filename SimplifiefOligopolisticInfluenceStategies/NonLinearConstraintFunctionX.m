function [C,Ceq] = NonLinearConstraintFunctionX(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function encodes the non-linear constraints which are imposed on the
% combined optimization problem which arises within the context of the
% Simplified Oligopolistic Model of Optimal Influences in the interaction
% network between two firms (Firm A and Firm B) and the single consumer C.

% The non-linear constraint function enforces the requirement that the
% underlying combined optimization objective accepts non-negative values
% the fact that its actual minimum is zero. This in turn, is an immediate
% consequence of the fact that the combined optimization objective is
% derived by the KKT conditions of the individual optimization problems.

% The main diversification of this constraint function is that it takes
% into consideration the existance of the upper-limit-related Lagrangian
% multipliers.

% Get the influence-related optimization variables.
TA = T(1);
TB = T(2);

% Get the Lagrangian multipliers-related optimization variables.
bA = T(3);
bB = T(4);

% Compute the values of the first order derivatives Da and Db.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Evaluate the combined optimization objective.
F = -(Da * TA + Db * TB) + bA * (1-TB) + bB * (1-TA);

% Since the non-linear inequality constraints should be of the form:
% C <= 0, we need to set C = -F such that C<=0 implies that -F <= 0 which
% results in F >=0.
C = -F;

% Mind that there exist no linear equality constraints for the problem at
% hand.
Ceq = [];

end

