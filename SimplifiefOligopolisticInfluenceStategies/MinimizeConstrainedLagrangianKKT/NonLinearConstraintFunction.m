function [Cineq,Ceq] = NonLinearConstraintFunction(T,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function codifies the non-linear equality and inequality constraints
% which are imposed on the combined optimization problem which arises in
% the context of the Simplified Oligopolistic Model of Optimal Influences
% within the interaction network between two firms (Firm A and Firm B) and 
% the single consumer C.

% The input vector T is assumed to be of the form:
% T = [TA TB MUa MUb MUab].

% Get the influence-related variables (the actual optimization variables).
TA = T(1);
TB = T(2);

% Get the lagrangian-related variables that are common for both cases.
% MUa  = T(3);
% MUb  = T(4);
% Get the lagrangian-related variables for CASE I:
if(beta >=0)
    % MUab = T(5);
    MUab = T(3);
% Get the lagrangian-related variables for CASE II:
else
    % MUab1 = T(5);
    % MUab2 = T(6);
    % MUab3 = T(7);
    MUab1 = T(3);
    MUab2 = T(4);
    MUab3 = T(5);
end

% Compute the values of the profit functions for the two firms Fa and Fb.
Fa = FirmAProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Fb = FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Compute the values of the first order derivatives Da and Db.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TB);

% Evaluate the combined optimization objective for CASE I:
if(beta >= 0)
    % DLa = Da + MUa - MUab;
    % DLb = Db + MUb - MUab;
    DLa = Da - MUab;
    DLb = Db - MUab;
else
% Evaluate the combined optimization objective for CASE II:
    % DLa = Da + MUa - MUab1 + alpha * LB * MUab2 + beta * LB * MUab3;
    % DLb = Db + MUb - MUab1 + beta * LA * MUab2 + alpha * LA * MUab3;
    DLa = Da - MUab1 + alpha * LB * MUab2 + beta * LB * MUab3;
    DLb = Db - MUab1 + beta * LA * MUab2 + alpha * LA * MUab3;
end

% This function implements the [D.5] and [D.6] non-linear inequality 
% constraints for the case of beta >=0 or the [Bo.7] and [Bo.8] for the 
% case of beta < 0 defined in  the memo file MinimizeConstrainedLagrangianKKT.
% Two additional non-linear inquality constraints may be incorporated in
% order to ensure that the optimal profit values will be non-negative.
% Mind that the non-linear inequality constraints with respect to some 
% variable x are expressed in a row-wise manner as C(x) <= 0. This is why 
% quantities DLa and DLb should be negated.
Cineq = [];
% Cineq = [-Fa;-Fb];
% Cineq = [-Fa;-Fb;-DLa;-DLb];
% Cineq = [-DLa;-DLb];

% CASE I: beta >= 0
% This function also implements the [D.2], [D.3] and [D.4]
% non-linear equality constraints defined in the memo file 
% MinimizeConstrainedLagrangianKKT. 
if (beta >= 0)
    % Ceq = [MUa * TA;MUb * TB;MUab * (1-TA-TB)];
    Ceq = MUab * (1-TA-TB);
else
% CASE II: beta < 0
% This function also implements the [Bo.2], [Bo.3], [Bo.4], [Bo.5] and [Bo.6]
% non-linear equality constraints defined in the memo file 
% MinimizeConstrainedLagrangianKKT.    
    % Ceq = [MUa * TA;MUb * TB;MUab1 * (1-TA-TB);...
           % MUab2 * (alpha * LB * TA + beta * LA * TB + LA * LB * (alpha * PA + beta * PB));...
           % MUab3 * (beta * LB * TA + alpha * LA * TB + LA * LB * (alpha * PB + beta * PA))];
    Ceq = [MUab1 * (1-TA-TB);...
           MUab2 * (alpha * LB * TA + beta * LA * TB + LA * LB * (alpha * PA + beta * PB));...
           MUab3 * (beta * LB * TA + alpha * LA * TB + LA * LB * (alpha * PB + beta * PA))];
end
% Mind that the non-linear equality constraints with respect to some 
% variable x are expressed in a row-wise manner as Ceq(x) = 0.

end