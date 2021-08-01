% This script file tests the efficiency of the 
% ExtractClusterRepresentativeSolutions function in identifying the distinct 
% solutions within a given set of solutions which are actually slight
% modifucations of a restricted set of different solutions.

% Set the maximum number of iterations.
MaxIterations = 100;
% Set the minimum number of different solutions.
Kmin = 1;
% Set the maximum number of different solutions.
Kmax = 10;
% Set the minimum number of dimensions per point.
Dmin = 2;
% Set the maximum number of dimensions per point.
Dmax = 10;
% Set the minimum number of solutions.
Nmin = 10;
% Set the maximum number of solutions.
Nmax = 20;
% Set the minimum value for the DiversificationFactor.
DiversificationFactorMin = 2;
% Set the maximum value for the DiversificationFactor.
DiversificationFactorMax = 10;
% Mind that the DiversificationFactor controls the number of digits that
% will be the same within each group of slightly different solutions.

% Perform main testing iteration:
for i = 1:1:MaxIterations
    % Draw the actual number of different solutions.
    K = randi([Kmin Kmax],1);
    % Draw the actual numbers of slightly different solutions per solution.
    N = randi
    % Draw the actual number of dimensions per point.
    D = randi([Dmin,Dmax],1);
    % Draw the actual diversification factor.
    DiversificationFactor = randi([DiversificationFactorMin,DiversificationFactorMax],1);
    % Compute the actual Diversification.
    Diversification = 10^(-DiversificationFactor);
    % Generate the initial data points.
    Po = rand(K,D);
end;