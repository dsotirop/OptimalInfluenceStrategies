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

tic
% Firstly, we need to define the fundamental parameters of the particular
% instance of the Simplified Oligopolistic game.
LA = 0.25;
LB = 0.25;
PA = 0.1;
PB = 0.1;
M = 0.4;
K = 0.1;
C = 0;

% Found two equilibrium points for LA=LB=0.25, PA=0.1, PB=0.3, M=K=0.3, 
% C = 0.0001 and G = 0.2

% Found three equilibrium points for LA=LB=0.25, PA=PB=0.1, M=K=0.3, 
% C = 0.0001 and G = 0.2

% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);
% Determine the minimum value for G.
% Mind that G is the Gamma parameter defined elsewhere.
Fmax = max(alpha^2,(beta^2+alpha*beta));
Lmin = min(LA,LB);
% Mind that G is the Gamma parameter defined elsewhere.
G = Fmax  / (Lmin^2);
% Set the fraction of the minimum G to be used.
Gfraction = 1.00; % This value should be set to 1.
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

% In order to characterize each point of the reaction sets Ra and Rb it is
% required to evaluate the first and second order derivatives. The first 
% order derivatives for each point within the reaction curves Ra and Rb will 
% be stored by matrices FDa and FDb in the following way:
%
% FDa(p) = {(Da(p),Db(p)): for all p = (TA*(t),t) in Ra} [5].
% FDb(p) = {(Da(p),Db(p)): for all p = (t,TB*(t)) in Rb} [6].
% The second order derivatives for each point within the reaction curves Ra  
% and Rb will be stored by matrices SDa and SDb in the following way:
%
% SDa(p) = {(DDa(p),DDb(p)): for all p = (TA*(t),t) in Ra} [7].
% SDb(p) = {(DDa(p),DDb(p)): for all p = (t,TB*(t)) in Rb} [8].

% Set the number of solution to be used by the optimizer.
SolutionsNumber = 200;

% Define the dt parameter controlling the density of the TRange interval.
dt = 0.01;
% Define the corresponding TRange interval.
TRange = 0.0:dt:1.0;
% Initialize curve containers Ra and Rb.
Ra = [];
Rb = [];
% Initialize the objective functions values containers Fa and Fb
% for the points pertaining to the reaction curves Ra and Rb respectively.
Fa = [];
Fb = [];
% Initialize the first order derivative containers FDa and FDb for the  
% points pertaining to the reaction curves Ra and Rb respectively.
FDa = [];
FDb = [];
% Initialize the second order derivative containers SDa and SDb for the  
% points pertaining to the reaction curves Ra and Rb respectively.
SDa = [];
SDb = [];
% Initialize the revenue functions containers Rev_a and Rev_b for the
% points pertaining to the reaction curves Ra and Rb respectively.
Rev_a = [];
Rev_b = [];
% Initialize the first order derivative containers D_Rev_a and D_Rev_b for
% the points pertaining to the reaction curves Ra and Rb respectively.
D_Rev_a = [];
D_Rev_b = [];
% Initialize the second order derivative containers D2_Rev_a and D2_Rev_b for
% the points pertaining to the reaction curves Ra and Rb respectively.
D2_Rev_a = [];
D2_Rev_b = [];


% Loop through the various values for both independent parameters TA and
% TB.
for t = TRange
    % Set the values of TA and TB for the reaction sets Ra and Rb
    % respectively.
    TA = t;
    TB = t;
    % Get the vector of optimal influence levels TAopt by solving the
    % optimization problem defined in [1] for the given value of TB = t.
    fprintf('Estimating Best Response TAopt when TB = %d\n',TB);
    [Raa,~,~]  = FirmABestResponse(SolutionsNumber,TB,C,G,LA,LB,PA,PB,alpha,beta,gamma);
    % Get the vector of optimal influence levels TBopt by solving the
    % optimization problem defined in [2] for the given value of TA = t.
    fprintf('Estimating Best Response TBopt when TA = %d\n',TA);
    [Rbb,~,~] = FirmBBestResponse(SolutionsNumber,TA,C,G,LA,LB,PA,PB,alpha,beta,gamma);
    
    % Loop through the various optimal influence levels stored in Raa:
    for ka = 1:1:length(Raa)
        % For each solution of the optimization problem [1] append the 
        % point (TAopt,TB) to the Ra curve.
        % For each point (TAopt,TB) append the values of the objective
        % functions fa(TAopt,TB) and fb(TAopt,TB) to the set Fa.
        % For each point (TAopt,TB) append the first order derivatives
        % Da(TAopt,TB) and Db(TAopt,TB) to the set FDa.
        % For each point (TAopt,TB) append the second order derivatives
        % DDa(TAopt,TB) and DDb(TAopt,TB) to the set SDa. 
        % For each point (TAopt,TB) append the values of the revenue
        % functions ra(TAopt,TB) and rb(TAopt,TB) to the set Rev_a.
        % For each point (TAopt,TB) append the values of the first order
        % derivatives dra(TAopt,TB) and drb(TAopt,TB) to the set D_Rev_a.
        % For each point (TAopt,TB) append the value of the second order
        % derivatives d2ra(TAopt,TB) and d2rb(TAopt,TB) to the set
        % D2_Rev_a.
        TAopt = Raa(ka);
        Ra = [Ra;[TAopt,TB]];
        fa = FirmAProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        fb = FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        ra = FirmARevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        rb = FirmBRevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        dra = FirmARevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        drb = FirmBRevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        d2ra = FirmARevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        d2rb = FirmBRevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        Fa  = [Fa;[fa,fb]];
        FDa = [FDa;[Da,Db]];
        SDa = [SDa;[DDa,DDb]];
        Rev_a = [Rev_a;[ra rb]];
        D_Rev_a = [D_Rev_a;[dra drb]];
        D2_Rev_a = [D2_Rev_a;[d2ra d2rb]];
    end
    
    % Loop through the various optimal influence levels stored in Rbb:
    for kb = 1:1:length(Rbb)
        % For each solution of the optimization problem [2] append the 
        % point (TA,TBopt) to the Rb curve.
        % For each point (TA,TBopt) append the values of the objective 
        % functions fa(TA,TBopt) and fb(TA,TBopt) to the set Fb
        % For each point (TA,TBopt) append the first order derivatives
        % Da(TA,TBopt) and Db(TA,TBopt) to the set FDb.
        % For each point (TA,TBopt) append the second order derivatives
        % DDa(TA,TBopt) and DDb(TA,TBopt) to the set SDb.
        % For each point (TA,TBopt) append the values of the revenue
        % functions ra(TA,TBopt) and rb(TA,TBopt) to the set Rev_b.
        % For each point (TA,TBopt) append the values of the first order
        % derivatives dra(TA,TBopt) and drb(TA,TBopt) to the set D_Rev_b.
        % For each point (TA,TBopt) append the value of the second order
        % derivatives d2ra(TA,TBopt) and d2rb(TA,TBopt) to the set
        % D2_Rev_b.
        TBopt = Rbb(kb);
        Rb = [Rb;[TA,TBopt]];
        fa = FirmAProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        fb = FirmBProfit(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        ra = FirmARevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        rb = FirmBRevenue(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        dra = FirmARevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        drb = FirmBRevenueFirstDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        d2ra = FirmARevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        d2rb = FirmBRevenueSecondDerivative(C,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        Fb = [Fb;[fa fb]];
        FDb = [FDb;[Da,Db]];
        SDb = [SDb;[DDa,DDb]];
        Rev_b = [Rev_b;[ra rb]];
        D_Rev_b = [D_Rev_b;[dra drb]];
        D2_Rev_b = [D2_Rev_b;[d2ra d2rb]];
    end
end

% Plot the best response curves for Ra and Rb.
PlotBestResponseCurves(Ra,Rb);

% Plot the profit functions and corresponding first and second derivatives
% on the meshgrid generated by value pairs (TA,TB) in TRange x TRange.
PlotFundamentalFunctions3D(C,G,LA,LB,PA,PB,alpha,beta,gamma,TRange);
toc