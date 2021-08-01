function [FA,FB,FA_rev,FB_rev,FA_cost,FB_cost] = OligopolisticFOptimal(pA,pB,C,G,TA,TB)

% This function computes the optimal profits (FA and FB) for boths firms 
% (Firm A and Firm B) with respect to the corresponding optimal prices 
% (pA and pB), the optimal investment levels (TA and TB) and the external 
% optimization parameters C and G (Gamma). Moreover, the optimal profit 
% quantities are decomposed into the corresponding constituent parts which 
% are, namely, the revenues (FA_rev and FB_rev) and the investment costs
% (FA_cost and FB_cost).

% Compute the optimal profit the for the Firm A.
FA = (pA - C)^2 - G * TA^2;
% Compute the revenue and cost constituents for the profit of Firm A.
FA_rev = (pA - C)^2;
FA_cost = G * TA^2;
% Compute the optimal profit the for the Firm B.
FB = (pB - C)^2 - G * TB^2;
% Compute the revenue and cost constituents for the profit of Firm B.
FB_rev = (pB - C)^2;
FB_cost = G * TB^2;

end

