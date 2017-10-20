function [XA,XB] = OligopolisticXOptimal(P_A_A,P_A_1,P_A_2,P_A_B,P_B_A,P_B_1,P_B_2,P_B_B,SA,S1,S2,SB)

% This function computes the consensus beliefs for both products A and B with
% respect to the corresponding initial beliefs (PA and PB) and limiting 
% influences vectors (SA, S1, S2 and SB).
% Quantities (SA, S1, S2 and SB) and (XA and XB) take their optimal values.

% Form the vector of initial beliefs for product A.
PA = [P_A_A,P_A_1,P_A_2,P_A_B];

% Form the vector of initial beliefs for prodict B.
PB = [P_B_A,P_B_1,P_B_2,P_B_B];

% Form the vector of limiting influences for each participant in the
% network.
S  = [SA,S1,S2,SB];

% Compute the limiting belief for product A.
XA = PA * S';

% Compute the limiting belief for product B.
XB = PB * S';

end

