function [F] = NonLinearObjectiveFunctionX(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This code file defines the non-linear objective function for the combined
% optimization problem which is derived for the two-player continuous game
% with coupled constraints according to the methodology presented in
% "Algoritms of Informatics" and "Mathematics and Its Applications".
% -------------------------------------------------------------------------
% The main differentation of this function lies upon the fact that there 
% exist two additional optimization variables (bA,bB) which correspond to 
% the upper bound-related Lagrangian multipliers of the original constrained
% optimization problems associated with each firm.

% The input vector T is assumed to be of the form:
% T = [TA TB bA bB].

% Get the influence-related optimization variables.
TA = T(1);
TB = T(2);

% Get the Lagrangian multipliers-related optimization variables.
bA = T(3);
bB = T(4);

% The objective function of the combined minimization problem may be
% formulated as:
%
% F = -[Da * TA + Db * TB] + bA * (1-TB) + bB * (1-TA) [1] where:
% 
%               d fa(TA,TB)                       d fb(TA,TB)
% Da(TA,TB) = -------------- [2] and Db(TA,TB) = -------------- [3]
%                 d TA                               d TB 
%
% According to the derivations conducted within the scope of
% the script file "SymbolicSimplifiedOligopolisticInfluenceModelQuadratic.m"
% it is known that the first-order derivatives Da and Db with respect to TA
% and TB are fractional function of the of form:
%
%              Ua(TA,TB)                       Ub(TA,TB)
% Da(TA,TB) = ----------- [4] and Db(TA,TB) = ----------- [5]
%              V(TA,TB)                         V(TA,TB) 
%
% Ua and Ub are 4th degree polynomials with respect to TA and TB.
% V is a 3rd degree polynomial with respect to TA or TB.


% Compute the values of the first order derivatives Da and Db.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Evaluate the combined optimization objective.
F = -(Da * TA + Db * TB) + bA * (1-TB) + bB * (1-TA);

end

