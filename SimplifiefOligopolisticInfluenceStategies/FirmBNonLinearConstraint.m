function [Cineq,Ceq] = FirmBNonLinearConstraint(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB)

% This function implements the non-linear constraints which are associated
% with the optimization problem solved by Firm A. Specifically, the signle
% non-linear constraint aims at reassuring the non-negative of the optimal
% profit for the first firm. Mind that there exist no non-linear equality
% constraints.

% Set the non-linear inequlity constraint function. Mind that inequality
% constraints with respect to vector x are expressed in the form:
% Cineq(x) <= 0.
Cineq = (-FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB));

% Set the non-linear equality constraints.
Ceq = [];

end