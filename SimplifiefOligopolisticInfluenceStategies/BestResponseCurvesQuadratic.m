% This script file provides a graphical representation of the reaction
% curves for firms A and B according to the Simplified Oligopolistic
% Optimal Influence Strategies model.

% This .m file is actually an extension to the original BestResponseCurves
% version with the only differantion being the incorporation of a quadratic
% cost term on the profit functions of both firms. This, in turn, leads to 
% higher degree (4th degree) polynomials Ua(TA) and Ub(TB). 

% The reaction curves for both firms A and B are initially given by the 3rd 
% degree bivariate polynomials Ua(TA,TB) and Ub(TA,TB). These polynomials 
% may be utilized in order to obtain the reaction sets for any given TAo in
% [0,1] and TBo in [0,1]. Specifically, the reaction sets for both firms
% are defined according to the following equations:
% For all TBo in [0,1], Ra(TBo) = {TA in [0,1]: Ua(TA,TBo) = 0} [1]
% For all TAo in [0,1], Rb(TAo) = {TB in [0,1]: Ub(TAo,TB) = 0} [2]

% However, solving equation Ua(TA,TBo) = 0 with respect to TA reduces to 
% finding the roots of a third degree polynomial with respect to TA for a
% given value of TBo. Accordingly, solving equation Ub(TAo,TB) = 0 with 
% respect to TB reduces to finding the roots of a third degree polynomial 
% with respect to TB for a given value of TAo.

% Mind that polynomials will be represented through their coeffients which
% will be stored within variables Caa and Cbb. The corresponding monomial
% terms for each coefficient stored in vectors Caa and Cbb will be:
% (i):  {TA^4,TA^3,TA^2,TA,1}
% (ii): {TB^4,TB^3,TB^2,TB,1}


clc
clear all

% Firstly, we need to define the fundamental parameters of the particular
% instance of the Simplified Oligopolistic game.
LA = 0.25;
LB = 0.25;
PA = 0.4;
PB = 0.4;
M = 0.3;
K = 0.1;
C = 0.0001;
G = 0.2; % Mind that G is the Gamma parameter defined elsewhere.

% Found two equilibrium points for LA=LB=0.25, PA=0.1, PB=0.3, M=K=0.3, 
% C = 0.0001 and G = 0.2

% Found three equilibrium points for LA=LB=0.25, PA=PB=0.1, M=K=0.3, 
% C = 0.0001 and G = 0.2

% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);

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
% a singleton value since it is derived as the root of the polynomial
% equation Ua(TA*,t) = 0.

% The reaction curve for the second firm (Firm B) may be 
% formulated as:
%
%               Rb = {(t,TB*(t)), for all t in [0,1]} [4].
%
% Mind that t in equation [4] represents any given value for TA (TA = t).
% Special attention is alse required by the term TB*(t) which might not be
% a singleton value since it is derived as the root of the polynomial
% equation Ub(t,TB*) = 0.

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

% Mind that the utilized extended containers are formed by grouping
% together solution points without taking into consideration the fact that
% the constituents points are not real numbers.

% Define the dt parameter controlling the density of the TRange interval.
dt = 0.001;
% Define the corresponding TRange interval.
TRange = [0:dt:1.0];
% Initialize curve containers Ra and Rb.
Ra = [];
Rb = [];
% Initialize extended curve containers for Ra and Ra.
RaX = [];
RbX = [];
% Initialize the first order derivative containers FDa and FDb for the  
% points pertaining to the reaction curves Ra and Rb respectively.
FDa = [];
FDb = [];
% Initialize the extended first order derivative containers FDa and FDb   
% for the points pertaining to the reaction curves Ra and Rb respectively.
FDaX = [];
FDbX = [];
% Initialize the second order derivative containers SDa and SDb for the  
% points pertaining to the reaction curves Ra and Rb respectively.
SDa = [];
SDb = [];
% Initialize the extended second order derivative containers SDa and SDb   
% for the points pertaining to the reaction curves Ra and Rb respectively.
SDaX = [];
SDbX = [];

% Loop through the various values for both independent parameters TA
% and TB.
for t = TRange
    % Set the values of TA and TB for the reaction curves Rb and Ra
    % respectively.
    TA = t;
    TB = t;
    % Define the corresponding polynomials Ua and Ub through their
    % associated coefficients stored in vectors Caa and Cbb for the given
    % values of TB and TA respectively.
    Caa = CaaPolyQuadratic(C,G,LA,LB,PA,PB,TB,alpha,beta,gamma);
    Cbb = CbbPolyQuadratic(C,G,LA,LB,PA,PB,TA,alpha,beta,gamma);
    % Since Caa and Cbb define the particular forms of the polynomials 
    % TAopt = {TAo: Ua(TB)=0} and TBopt = {TBo: Ub(TA)=0}, we need to find
    % the roots for the particular forms of the polymials.
    Raa = roots(Caa);
    Rbb = roots(Cbb);
    % At this point we will take into consideration all roots that pertain
    % to the set of the real numbers irrespectively of the actual
    % contraints imposed to the underlying optimization problem.
    
    % Loop through the various roots of the polynomial Raa:
    for ka = 1:1:length(Raa)
        % For each real root append the point (TAopt,TB) to the Ra curve.
        % For each point (TAopt,TB) append the first order derivatives
        % Da(TAopt,TB) and Db(TAopt,TB) to the set FDa.
        % For each point (TAopt,TB) append the second order derivatives
        % DDa(TAopt,TB) and DDb(TAopt,TB) to the set SDa.
        if(isreal(Raa(ka)))
            TAopt = Raa(ka);
            Ra = [Ra;[TAopt,TB]];
            Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
            Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
            DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
            DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
            FDa = [FDa;[Da,Db]];
            SDa = [SDa;[DDa,DDb]];
        end
        % Populate the corresponding extended containers irrepsective of
        % whether the associated roots are real numbers.
        TAopt = Raa(ka);
        RaX = [RaX;[TAopt,TB]];
        Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TB);
        FDaX = [FDaX;[Da,Db]];
        SDaX = [SDaX;[DDa,DDb]];
    end
    
    % Loop through the various roots of the polynomial Rbb:
    for kb = 1:1:length(Rbb)
        % For each real root append the point (TA,TBopt) to the Rb curve.
        % For each point (TA,TBopt) append the first order derivatives
        % Da(TA,TBopt) and Db(TA,TBopt) to the set FDb.
        % For each point (TA,TBopt) append the second order derivatives
        % DDa(TA,TBopt) and DDb(TA,TBopt) to the set SDb.
        if(isreal(Rbb(kb)))
            TBopt = Rbb(kb);
            Rb = [Rb;[TA,TBopt]];
            Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
            Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
            DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
            DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
            FDb = [FDb;[Da,Db]];
            SDb = [SDb;[DDa,DDb]];
        end
        % Populate the corresponding extended containers irrepsective of
        % whether the associated roots are real numbers.
        TBopt = Rbb(kb);
        RbX = [RbX;[TA,TBopt]];
        Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TA,TBopt);
        FDbX = [FDbX;[Da,Db]];
        SDbX = [SDbX;[DDa,DDb]];
    end
end

% Plot the best response curves for Ra and Rb.
PlotBestResponseCurves(Ra,Rb);

% Plot the profit functions and corresponding first and second derivatives
% on the meshgrid generated by value pairs (TA,TB) in TRange x TRange.
PlotFundamentalFunctions3D(C,G,LA,LB,PA,PB,alpha,beta,gamma,TRange);
