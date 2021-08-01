function [F] = NonLinearObjectiveFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This code file defines the non-linear objective functio for the combined
% optimization problem which is derived for the two-player continuous game
% with coupled constraints involving two firms (Firm A and Firm B) and a
% sigle consumer (Consumer C).

% This function implements the [D.1] combined augmemnted objective function
% defined in the memo file MinimizeConstrainedLagrangianKKT.

% The combined objective function codifies the KKT requirements which 
% correspond to the following equations: 
%
%  d La(TA,TB)              d Lb(TA,TB)   
% ------------ = 0 [I] and ------------ = 0 [II]
%    d TA                     d TB 
% -------------------------------------------------------------------------
% CASE I: beta >= 0
% -------------------------------------------------------------------------
% The optimization procedure attempts to minimize the sum of the quantities
% defined in Eqs.(I) and (II) which is equivalent to:
%                         2                          2
% [Da(TA,TB) + MUa - MUab] + [Db(TA,TB) + MUb - MUab] [III]
%
% The input vector T is assumed to be of the form:
% T = [TA TB MUa MUb MUab].
% -------------------------------------------------------------------------
% CASE II: beta < 0
% -------------------------------------------------------------------------
% The optimization procedure attempts to minimize the sum of the quantities
% defined in Eqs.(I) and (II) which is equivalent to:
%                                                                   2                          
% [Da(TA,TB) + MUa - MUab1 + alpha * LB * MUab2 + beta * LB * MUab3] +                     2 
%                                                                   2
% [Db(TA,TB) + MUb - MUab1 + beta * LA * MUab2 + alpha * LA * MUab3] [III]
%
% The input vector T is assumed to be of the form:
% T = [TA TB MUa MUb MUab1 MUab2 MUab3].


% Get the influence-related variables (the actual optimization variables).
TA = T(1);
TB = T(2);

% Get the lagrangian-related variables that are common for both cases.
% MUa  = T(3);
% MUb  = T(4);
% Get the lagrangian-related variables for CASE I:
if(beta >=0)
    MUab = T(3);
% Get the lagrangian-related variables for CASE II:
else
    MUab1 = T(3);
    MUab2 = T(4);
    MUab3 = T(5);
end
    
% Compute the values of the first order derivatives Da and Db.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Evaluate the combined optimization objective for CASE I:
if(beta >= 0)
    % DLa = Da + MUa - MUab;
    DLa = Da - MUab;
    % DLb = Db + MUb - MUab;
    DLb = Db - MUab;
else
% Evaluate the combined optimization objective for CASE II:
    % DLa = Da + MUa - MUab1 + alpha * LB * MUab2 + beta * LB * MUab3;
    % DLb = Db + MUb - MUab1 + beta * LA * MUab2 + alpha * LA * MUab3;
    DLa = Da - MUab1 + alpha * LB * MUab2 + beta * LB * MUab3;
    DLb = Db - MUab1 + beta * LA * MUab2 + alpha * LA * MUab3;
end
% Evaluate the combined optimization objective.
F = DLa^2 + DLb^2;
end