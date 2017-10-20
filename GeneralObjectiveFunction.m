function [F] = GeneralObjectiveFunction(T,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma)

% This function provides the most general implementation of the objective 
% function of the minimization problem that corresponds to problem of determining the
% determining the optimal influences.
T1 = T(1);
T2 = T(2);
F = Gamma*(T1 + T2) - (Delta^2*((P0*T2 + P0*Theta1 + Lambda2*P2)*T1 +...
    + Lambda1*P1*T2 + Lambda1*P1*Theta1 + Lambda2*P1*Theta1 + Lambda1*P2*Theta2 + ...
    + Lambda2*P2*Theta2 + P0*T2*Theta2))/(2*((Lambda2 + T2 + Theta1)*T1 + ...
    + Lambda1*Theta1 + Lambda1*Theta2 + Lambda2*Theta1 + Lambda2*Theta2 +...
    + T2*Theta2 + Lambda1*T2));

end

