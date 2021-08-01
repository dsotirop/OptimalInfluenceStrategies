function [QA,QB] = OligopolisticQOptimal(XA,XB,alpha,beta,gamma_prime)

% This function computes the optimal quantities (QA and QB) for the products 
% of firms A and B respectively given the agents' limiting beliefs for each
% product (XA and XB) and the external optimization parameters alpha, beta
% and gamma'.

% Compute the optimal quantity for the product of Firm A.
QA = alpha * XA + beta * XB - gamma_prime;

% Compute the optimal quantity for the product of Firm B.
QB = beta * XA + alpha * XB - gamma_prime;

end

