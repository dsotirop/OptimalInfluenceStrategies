% This script file provides an algorithmic estimation of the Best Response
% Curves for the firms (A) and (B) that pertain to the Simplified Oligopolistic 
% Optimal Influece model where both firms compete in order to determine their
% optimal influence level on the single consumer (C) of the considered network.

% Estimation of the Best Reaction Curves for each firm in the network will
% be conducted by solving an instance of the associated maximization tasks
% given the influence level selected by the opponent firm. Specifically,
% the reaction sets for both firms are defined according to the following
% equations:
% For all TBo in [0,1], Ra(TBo) = {TAopt in [0,1]: min fa(TA,TBo) s.t. 0<= TA <=1-TB} [1] 
%                                                  TA
% For all TAo in [0,1], Rb(TAo) = {TBopt in [0,1]: min fb(TAo,TB) s.t. 0<= TB <=1-TA} [2] 
%                                                  TB

clc
clear

% Firstly, we need to define the fundamental parameters of the particular
% instance of the Simplified Oligopolistic game.
LA = 0.50;
LB = 0.50;
PA = 0.50;
PB = 0.50;
M = 0.30;
K = 0.20;
C = 0;

% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);
gamma_prime = gamma * (M - 1); % gamma_prime is the gamma' parameter.
% Determine the minimum value for G.
% Mind that G is the Gamma parameter defined elsewhere.
Fmax = max(alpha^2,(3*beta^2+2*alpha*beta));
Lmin = min(LA,LB);
% Mind that G is the Gamma parameter defined elsewhere.
G = Fmax  / (Lmin^2);
% Set the fraction of the minimum G to be used.
Gfraction = 0.50; % This value should be set to 1.
G = Gfraction * G;

% Construct the pairs of points that define the best reaction curves Ra, Rb
% for firms A and B according to equations [1] and [2]. The reaction curve
% for the first firm will be stored in a matrix Ra which will be holding
% pairs of (TAopt,TB) values for any given value of TB. Accordingly, the
% reaction curve for the second firm will be stored in a matrix Rb which
% will be holding pairs of (TA,TBopt) values for any given value of TA.
% For both curves the indepented variables TA and TB will be taking values
% within the TRange interval [0:dt:1] whose density will be parameterized
% by dt.

% The reaction curve for the first firm (Firm A) may be 
% formulated as:
%
%               Ra = {(TA*(t),t), for all t in [0,1]} [3].
%
% Mind that t in equation [3] represents any given value for TB (TB = t).
% Special attention is alse required by the term TA*(t) which might not be
% a singleton value since it is derived by the representative solution
% points of the optimization problem defined in [1].

% The reaction curve for the second firm (Firm B) may be 
% formulated as:
%
%               Rb = {(t,TB*(t)), for all t in [0,1]} [4].
%
% Mind that t in equation [4] represents any given value for TA (TA = t).
% Special attention is alse required by the term TB*(t) which might not be
% a singleton value since it is derived by the representative solution
% points of the optimization problem defined in [2].

% Set the number of solution to be used by the optimizer.
SolutionsNumber = 2000;
% Define the dt parameter controlling the density of the TRange interval.
dt = 0.005;
% Define the corresponding TRange interval.
TRange = 0.00:dt:1.00;
% Initialize curve containers Ra and Rb.
Ra = zeros(length(TRange),2);
Rb = zeros(length(TRange),2);
% Initialize FilterFlag container FFa and FFb.
FFa = zeros(length(TRange),1);
FFb = zeros(length(TRange),1);

% Loop through the various values for both independent parameters TA and TB
parfor k = 1:length(TRange)
    % Set the current value for the direct influence levels of each 
    % adverserial firm.
    t = TRange(k);
    TA = t;
    TB = t;
    % Get the vector of optimal influence level TAopt by solving the
    % optimization problem defined in [1] for the given value of TB = t.
    fprintf('Estimating Best Response TAopt when TB = %d\n',TB);
    [TAopt,~,FilterFlagA] = FirmABestResponse(SolutionsNumber,TB,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
    % Update the corresponding entry within the matrix storing the reaction
    % curve for Firm A.
    Ra(k,:) = [TAopt,TB];
    % Update the corresponding FilterFlag container for the optimization
    % problem solved for the first firm.
    FFa(k) = FilterFlagA;
    % Get the vector of optimal influence level TBopt by solving the
    % optimization problem defined in [2] for the given value of TA = t.
    fprintf('Estimating Best Response TBopt when TA = %d\n',TA);
    [TBopt,~,FilterFlagB] = FirmBBestResponse(SolutionsNumber,TA,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
    % Update the corresponding entry within the matrix storing the reaction
    % curve for Firm B.
    Rb(k,:) = [TA,TBopt];
    % Update the corresponding FilterFlag container for the optimization
    % problem solved for the second firm.
    FFb(k) = FilterFlagB;
end

% Plot the best response curves for Ra and Rb.
PlotBestResponseCurves(Ra,Rb);
% Plot the profit functions and corresponding first and second derivatives
% on the meshgrid generated by value pairs (TA,TB) in TRange x TRange.
PlotFundamentalFunctions3D(C,G,LA,LB,PA,PB,alpha,beta,gamma,TRange);