function [SD] = FirmBOptimalityConditions(Solutions,TA,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function evaluates the First and Second Order optimality conditions   
% for each pair of solutions of the form (TA,TBopt(TA)), indicating that
% TBopt(TA) is Firm's B best response to a given investement level of
% Firm A. In this setting, Solutions is column vector storing the best 
% response solutions for the investement level of Firm B as acquired by the 
% execution of the function "FirmBProfitMaximizer". In this context, TA
% stores just a scalar value.
%
%            D(k) = [DDb(TA,TBopt)].
% Mind that: 
% 
% DDb: is the second order derivative of fb with respect to TB
% Solutions: is a column vector storing best response solutuions for Firm B 

% Set the optimal solutions vector for the investment level of Firm B.
TBopt = Solutions;
% Construct a vector of the same dimensionality as TBopt storing the
% investement level of Firm A.
TAo = repmat(TA,size(TBopt,1),1);

% Evaluate Second Order Derivatives.
SD = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAo,TBopt);

end