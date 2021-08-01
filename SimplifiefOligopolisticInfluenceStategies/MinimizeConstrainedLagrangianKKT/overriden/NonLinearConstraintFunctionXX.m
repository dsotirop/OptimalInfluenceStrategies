function [C,Ceq] = NonLinearConstraintFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function codifies the non-linear equality and inequality constraints
% which are imposed on the combined optimization problem which arises in
% the context of the Simplified Oligopolistic Model of Optimal Influences
% within the interaction network between two firms (Firm A and Firm B) and 
% the single consumer C.

% The input vector T is assumed to be of the form:
% T = [TA TB LambdaA_1 LambdaA_2 LambdaB_1 LambdaB_2].

% Get the influence-related variables (the actual optimization variables).
TA = T(1);
TB = T(2);

% Get the lagrangian-related variables with respect to TA.
LambdaA_1 = T(3);
LambdaA_2 = T(4);

% Get the lagrangian-related variables with respect to TB.
LambdaB_1 = T(5);
LambdaB_2 = T(6);

% Compute the values of the first order derivatives Da and Db.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Evaluate the combined optimization objective.
DLa = Da + LambdaA_1 - LambdaA_2;
DLb = Db + LambdaB_1 - LambdaB_2;

% This function implements the [C.6] and [C.7] non-linear inequality 
% constraints defined in the memo file MinimizeConstrainedLagrangianKKT.
% Mind that the non-linear inequality constraints with respect to some 
% variable x are expressed in a row-wise manner as C(x) <= 0. This is why 
% quantities DLa and DLb should be negated.
C = [-DLa;-DLb];

% This function also implements the [C.2], [C.3], [C.4] and [C.5]
% non-linear equality constraints defined in the memo file 
% MinimizeConstrainedLagrangianKKT. Mind that the non-linear equality 
% constraints with respect to some variable x are expressed in a row-wise
% manner as Ceq(x) = 0.
Ceq = [LambdaA_1 * TA;LambdaA_2 * (1-TB-TA);LambdaB_1 * TB;LambdaB_2 * (1-TA-TB)];


end

