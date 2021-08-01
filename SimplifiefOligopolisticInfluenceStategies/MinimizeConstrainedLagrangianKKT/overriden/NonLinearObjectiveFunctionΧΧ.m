function [F] = NonLinearObjectiveFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This code file defines the non-linear objective functio for the combined
% optimization problem which is derived for the two-player continuous game
% with coupled constraints involving two firms (Firm A and Firm B) and a
% sigle consumer (Consumer C).

% This function implements the [C.1a] or [C.1b] alternatives defined in the
% memo file MinimizeConstrainedLagrangianKKT.

% The combined objective function codifies the KKT requirements which 
% correspond to the following equations: 
%
%  d La(TA,TB)              d Lb(TA,TB)   
% ------------ = 0 [I] and ------------ = 0 [II]
%    d TA                     d TB 
%
% The optimization procedure attempts to minimize the sum of the quantities
% defined in Eqs.(I) and (II) which is equivalent to:
% [Da(TA,TB) + LamdaA_1 - LambdaA_2] + [Db(TA,TB) + LamdaB_1 - LambdaB_2] [III]

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
% Uncomment the following code line in case the [C.1a] is to be evaluated.
%F = DLa + DLb;
% Uncomment the following code line in case the [C.1b] is to be evaluated.
% Evaluate the alternative for the combined optimization objective.
F = DLa^2 + DLb^2;
end