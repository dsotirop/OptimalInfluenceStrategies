function [F] = NonLinearObjectiveFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function defines the most general form of the system of nonlinear
% equations that has to be solved in order to obtain the optimal influence
% strategies within the simplified oligopolostic environment of the two
% firms {Firm A, Firm B} and one consumer {C}. However, the system of
% nonlinear equations will be reformulated as a nonlinear function
% optimization problem.

% The input vector T is assumed to be of the form:
% T = [TA TB]

% Get the influence related variables (the actual optimization varibles).
TA = T(1);
TB = T(2);

% The system of nonlinear equations corresponds to the first order
% conditions on the profit functions fa and fb for products A and B. In 
% fact, we have a nonlinear system of the following form:
%          [F1 F2] = [F(1) F(2)] =  [0 0] [I]
% where
%               d fa                  d Fb
%       F(1) = --------  and F(2) = -------- [II] 
%               d TA                  d Tb

F(1) = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
F(2) = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Finally, the nonlinear objective function to be minimized will be given
% as:
%                                    2         2     
%                      F(T) = (F1 - 0)+ (F2 - 0) 
%
% where each Fi with i in {1,2} corresponds to a particular equation of
% the original system.

F = sum(F.^2);

end

