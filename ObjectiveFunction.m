function [F] = ObjectiveFunction(T,P1,P2,Theta,Delta,Gamma)

% This function provides the implementation of the objective function of
% the minimization problem that corresponds to problem of determining the
% optimal influences.

% Mind that the Delta parameters corresponds to the difference (A - C).

% Setup the objective function.
T1 = T(1);
T2 = T(2);
R = (P1+2*T1)*Theta + P2 + 2*T2 + P1*T2 + P2*T1 + 4*T1*T2;
Q = (2*T1+1)*Theta + T1 + 3*T2 + 4*T1*T2 + 1;
F = Gamma*(T1+T2)- (Delta^2/2)* (R/Q);

end

