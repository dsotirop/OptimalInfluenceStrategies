function [S0,S1,S2] = Soptimal(T1opt,T2opt,Theta)

% This function returns the optimal limiting influences with respect to 
% the optimal T1 and T2 for a given Theta.

Denom = T1opt + 3*T2opt + Theta + 4*T1opt*T2opt + 2*Theta*T1opt + 1;
Nom1 = 2*(T2opt+2*T1opt*T2opt+Theta*T1opt);
Nom2 = Theta + T2opt;
Nom3 = 1 + T1opt;
S0 = Nom1 / Denom;
S1 = Nom2 / Denom;
S2 = Nom3 / Denom;

end

