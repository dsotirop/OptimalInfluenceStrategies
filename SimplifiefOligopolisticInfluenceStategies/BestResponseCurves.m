% This script file provides a graphical representation of the reaction
% curves for firms A and B according to the Simplified Oligopolistic
% Optimal Influence Strategies model.

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
% (i):  {TA^3,TA^2,TA,1}
% (ii): {TB^3,TB^2,TB,1}

clc
clear all

% Firstly, we need to define the fundamental parameters of the particular
% instance of the Simplified Oligopolistic game.
LA = 0.5;
LB = 0.5;
PA = 0.2;
PB = 0.2;
M = 0.1;
K = 0.1;
C = 0.0;
G = 0.3; % Mind that G is the Gamma parameter defined elsewhere.

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

% Define the dt parameter controlling the density of the TRange interval.
dt = 0.01;
% Define the corresponding TRange interval.
TRange = [0:dt:1.0];
% Initialize curve containers Ra and Rb.
Ra = [];
Rb = [];
% Loop through the various values for both the independent parameters TA
% and TB.
for t = TRange
    % Set the values of TA and TB for the reaction curves Rb and Ra
    % respectively.
    TA = t;
    TB = t;
    % Define the corresponding polynomials Ua and Ub through their
    % associated coefficients stored in vectors Caa and Cbb for the given
    % values of TB and TA respectively.
    Caa = CaaPolyX(C,G,LA,LB,PA,PB,TB,alpha,beta,gamma);
    Cbb = CbbPolyX(C,G,LA,LB,PA,PB,TA,alpha,beta,gamma);
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
        if(isreal(Raa(ka)))
            TAopt = Raa(ka);
            Ra = [Ra;[TAopt,TB]];
        end
    end
    
    % Loop through the various roots of the polynomial Rbb:
    for kb = 1:1:length(Rbb)
        % For each real root append the point (TA,TBopt) to the Rb curve.
        if(isreal(Rbb(kb)))
            TBopt = Rbb(kb);
            Rb = [Rb;[TA,TBopt]];
        end
    end
end

% Plot the best response curves for Ra and Rb.
PlotBestResponseCurves(Ra,Rb);