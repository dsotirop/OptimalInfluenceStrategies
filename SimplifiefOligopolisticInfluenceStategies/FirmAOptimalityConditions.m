function [SD] = FirmAOptimalityConditions(Solutions,TB,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function evaluates the First and Second Order optimality conditions   
% for each pair of solutions of the form (TAopt(TB),TB), indicating that
% TAopt(TB) is Firm's A best response to a given investement level of
% Firm B. In this setting, Solutions is column vector storing the best 
% response solutions for the investement level of Firm A as acquired by the 
% execution of the function "FirmAProfitMaximizer" whereas TB is just a 
% scalar value.
%
%               SD(k) = [DDa(TAopt,TB)].
% Mind that:
%
% DDa: is the second order derivative of fa with respect to TA
% Solutions: is a column vector storing best response solutuions for Firm A 

% Set the optimal solutions vector for the investment level of Firm A.
TAopt = Solutions;
% Construct a vector of the same dimensionality as TAopt storing the
% investement level of Firm B.
TBo = repmat(TB,size(TAopt,1),1);

% Evaluate Second Order Derivatives.
SD = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBo);

end