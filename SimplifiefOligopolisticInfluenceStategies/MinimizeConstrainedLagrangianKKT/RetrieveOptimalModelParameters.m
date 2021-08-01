function [Sopt,Xopt,Popt,Qopt,Fopt] = RetrieveOptimalModelParameters(Topt,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime)

% This function computes the model parameter values at the optimal solution 
% points stored in matrix Topt.

% Topt stores in a row-wise manner optimal [TAopt TBopt] pairs
% corresponding to the optimal investment levels.
% Sopt stores in a row-wise manner optimal [SAopt SCopt SBopt] triplets
% corresponding to the optimal limiting influence levels for the three
% agents of the networl {FirmA, Consumer C, FirmB}.
% Xopt stors in a row-wise manner optimal [XAopt,XBopt] pairs corresponding
% to the optimal limiting beliefs for the two products A and B.
% Popt stores in a row-wise manner optimal [pAopt pBopt] pairs corrsponding
% to the optimal prices for products A and B.
% Qopt stores in a row-wise manner optimal [QAopt QBopt] pairs corrsponding
% to the optimal quantities for products A and B.
% Fopt stores in a  row-wise manner optimal
% [FAopt,FBopt,FA_rev_opt,FB_rev_opt,FA_cost_opt,FB_cost_opt] vectors
% corresponding to:
% (i): the overall profit FAopt and FBopt for the two firms.
% (ii): the net profit (revenue) FA_rev_opt and FB_rev_opt for the two firms.
% (iii): the cost FA_cost_opt and FB_cost_opt for the two firms.

% Get the number of solution points.
N = size(Topt,1);

% Initialize matrix containers Sopt, Xopt, Popt, Qopt and Fopt.
Sopt = zeros(N,3);
Xopt = zeros(N,2);
Popt = zeros(N,2);
Qopt = zeros(N,2);
Fopt = zeros(N,6);

% Loop through the various solution points.
for k = 1:1:N
    % Retrieve the current optimal solution points TA and TB.
    TA = Topt(k,1);
    TB = Topt(k,2);
    % Compute the optimal influence levels for the current optimal solution
    % point.
    [SA,SC,SB] = OligopolisticSOptimal(LA,LB,TA,TB);
    % Compute the limiting beliefs for the current optimal solution point.
    [XA,XB] = OligopolisticXOptimal(PA,PB,SA,SC,SB);
    % Compute the optimal quantities for the current optimal solution
    % point.
    [pA,pB] = OligopolisticPOptimal(XA,XB,alpha,beta,gamma);
    % Compute the optimal prices for the current optimal solution point.
    [QA,QB] = OligopolisticQOptimal(XA,XB,alpha,beta,gamma_prime);
    % Compute the optimal prices for the current optimal solution point.
    [FA,FB,FA_rev,FB_rev,FA_cost,FB_cost] = OligopolisticFOptimal(pA,pB,C,G,TA,TB);
    % Form the Sopt vector.
    Sopt(k,:)= [SA,SC,SB];
    % Form the Xopt vector.
    Xopt(k,:) = [XA,XB];
    % Form the Popt vector.
    Popt(k,:) = [pA,pB];
    % Form the Qopt vector.
    Qopt(k,:) = [QA,QB];
    % Form the Fopt vector.
    Fopt(k,:) = [FA,FB,FA_rev,FB_rev,FA_cost,FB_cost];
end

end