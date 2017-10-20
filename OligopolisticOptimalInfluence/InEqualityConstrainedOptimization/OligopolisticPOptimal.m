function [pA,pB] = OligopolisticPOptimal(XA,XB,C,K,M)

% This function computes the optimal prices (pA and pB) for products A and
% B with respect to limiting agents' beliefs for each product and the
% parameters C, K and M.

% Compute the optimal value for the price of product A.
pA = -(2*C + 2*XA - 2*K*XB + M*XB + C*M - K*M*XA)/(M^2 - 4);

% Compute the optimal value for the price of product B.
pB = -(2*C + 2*XB - 2*K*XA + M*XA + C*M - K*M*XB)/(M^2 - 4);


end

