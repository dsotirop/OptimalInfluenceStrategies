function [F] = NonLinearFunction(T,P1_const,P2_const,theta_const,A_const,C_const,G_const)

% This is the implementation of the non-linear function whose minimum
% with respect to T corrsponds to the objective of optimal influence
% problem.

% Setup the objective function.
T1 = T(1);
T2 = T(2);
R = (P1_const+2*T1)*theta_const + P2_const + 2*T2 + P1_const*T2 + P2_const*T1 + 4*T1*T2;
Q = (2*T1+1)*theta_const + T1 + 3*T2 + 4*T1*T2 + 1;
F = G_const*(T1+T2)- ((A_const-C_const)^2/2)* (R/Q);

end

