function [Q_A,Q_B] = OligopolisticQOptimal(XA,XB,C,K,M)

% This function computes the optimal quantities (Q_A and Q_B) for products 
% A and B with respect to limiting agents' beliefs for each product and the
% parameters C, K and M.

% Compute the optimal value for the quantity of product A.
Q_A = -(2*(2*XA - 2*C - 2*K*XB + M*XB + C*M^2 + C*M - K*M*XA))/(M^2 - 4);

% Compute the optimal value for the quantity of product B.
Q_B = -(2*(2*XB - 2*C - 2*K*XA + M*XA + C*M^2 + C*M - K*M*XB))/(M^2 - 4);


end

