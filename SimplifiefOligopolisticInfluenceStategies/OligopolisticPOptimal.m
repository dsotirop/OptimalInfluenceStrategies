function [pA,pB] = OligopolisticPOptimal(XA,XB,alpha,beta,gamma)

% This function computes the optimal prices (pA and pB) for the products of
% firms A and B respectively given the agents' limiting beliefs for each
% product (XA and XB) and the external optimization parameters alpha, beta
% and gamma.

% Compute the optimal price for the product of Firm A.
pA = alpha * XA + beta * XB - gamma;

% Compute the optimal price for the product of Firm B.
pB = beta * XA + alpha * XB - gamma;


end