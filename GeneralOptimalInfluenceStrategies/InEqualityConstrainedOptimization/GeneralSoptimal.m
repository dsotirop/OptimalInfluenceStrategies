function [S0,S1,S2] = GeneralSoptimal(T1opt,T2opt,Lambda1,Lambda2,Theta1,Theta2)

% This function returns the optimal limiting influences with respect to 
% the optimal T1 and T2 given Lambda1, Lambda2, Theta1 and Theta2.

Denom = Lambda1*Theta1 + Lambda1*Theta2 + Lambda2*Theta1 + Lambda2*Theta2 + T1opt*T2opt + T1opt*Theta1 + T2opt*Theta2 + Lambda1*T2opt + Lambda2*T1opt;
Nom0 = T1opt*T2opt + T1opt*Theta1 + T2opt*Theta2;
Nom1 = Lambda1*Theta1 + Lambda2*Theta1 + Lambda1*T2opt;
Nom2 = Lambda1*Theta2 + Lambda2*Theta2 + Lambda2*T1opt;
S0 = Nom0 / Denom;
S1 = Nom1 / Denom;
S2 = Nom2 / Denom;
end

