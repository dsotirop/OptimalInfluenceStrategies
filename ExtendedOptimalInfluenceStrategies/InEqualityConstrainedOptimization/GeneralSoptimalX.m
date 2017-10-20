function [S0,S1,S2] = GeneralSoptimalX(T1opt,T2opt,Lambda1opt,Lambda2opt,Theta1,Theta2)

% This function returns the optimal limiting influences with respect to 
% the optimal T1,T2,Lambda1 and Lambda2 with respect to Theta1 and Theta2.

Denom = Lambda1opt*Theta1 + Lambda1opt*Theta2 + Lambda2opt*Theta1 + Lambda2opt*Theta2 + T1opt*T2opt + T1opt*Theta1 + T2opt*Theta2 + Lambda1opt*T2opt + Lambda2opt*T1opt;
Nom0 = T1opt*T2opt + T1opt*Theta1 + T2opt*Theta2;
Nom1 = Lambda1opt*Theta1 + Lambda2opt*Theta1 + Lambda1opt*T2opt;
Nom2 = Lambda1opt*Theta2 + Lambda2opt*Theta2 + Lambda2opt*T1opt;
S0 = Nom0 / Denom;
S1 = Nom1 / Denom;
S2 = Nom2 / Denom;
end

