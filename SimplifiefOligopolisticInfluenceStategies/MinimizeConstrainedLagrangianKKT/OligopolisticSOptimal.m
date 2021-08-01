function [SA,SC,SB] = OligopolisticSOptimal(LA,LB,TA,TB)

% This function computes the optimal limiting influences SA, SC and SB for
% the Simplified Oligopolistic Network model involving two firms (Firm A
% and Firm B) and a single consumer (Consumer C). The optimal quantities
% for the limiting influences of the two firms (SA and SB) as well as the
% optimal quantity for the limiting influence of the consumer (SC) are
% obtained as a function of the optimal limiting investment levels TA and TB
% as well as the external optimization parameters LA and LB by the combined 
% optimization process which extracts the Nash equilibrium for the 
% continuous game between the two firms.

% Compute the common denominator for the three quantities SA, SC and SB.
D = (LA*LB + LA*TB + LB*TA);
% Compute the three different nominatores NA, NC and NB for the
% corresponding quantities SA, SC and SB.
NA = LB*TA;
NC = LA*LB;
NB = LA*TB;
% Compute the actual values for the three quantities SA, SC and SB.
SA = NA / D;
SC = NC / D;
SB = NB / D;

end

