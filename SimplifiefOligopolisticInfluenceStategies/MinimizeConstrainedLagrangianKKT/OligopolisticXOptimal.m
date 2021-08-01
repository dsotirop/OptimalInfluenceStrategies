function [XA,XB] = OligopolisticXOptimal(PA,PB,SA,SC,SB)

% This function computes the consensus beliefs XA and XB for both products 
% A and B with respect to the corresponding initial beliefs PA and PB of 
% consumer C and the optimal limiting influences of each agent in the network
% SA, SC and SB.

% Form the vector PPA = [PAA PA PBA] of initial agents' beliefs for the
% product of Firm A given that PAA = 1 and PBA = 0.
PPA = [1 PA 0];
% Form the vector PPB = [PAB PB PBB] of initial agents' beliefs for the 
% product of Firm B given that PBB = 1 and PAB = 0.
PPB = [0 PB 1];
% Form the vector of optimal limiting influences for each participant in the
% network.
S = [SA;SC;SB];
% Compute the optimal limiting beliefs of consumer C for the products of
% firms A and B.
XA = PPA * S;
XB = PPB * S;


end

