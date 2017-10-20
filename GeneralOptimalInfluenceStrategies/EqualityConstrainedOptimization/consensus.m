function X = consensus(P0,P1,P2,S0,S1,S2)

% This function computes the consensus belief with respect to P and S.
P = [P0,P1,P2];
S = [S0,S1,S2];
X = P*S';


end

