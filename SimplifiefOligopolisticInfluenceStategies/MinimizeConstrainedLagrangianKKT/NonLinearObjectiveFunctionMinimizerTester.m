% This script file provides the fundamental testing functionality for
% evaluating the process of solving the underlying constrained optimization 
% problem. 

% Clear workspace and command window
clc
clear

% Set the number of firms.
FirmsNumber = 2;

% Set of parameters for which no solution could be found.
% (1): [LA = 0.50,LB = 0.50,PA = 0.80,PB = 0.15,M = 0.10,K = 0.30,C = 0]
% (2): [LA = 0.50,LB = 0.50,PA = 0.80,PB = 0.20,M = 0.10,K = 0.30,C = 0]

% Initialize the external optimization variables.
LA = 0.50;
LB = 0.50;
PA = 0.50;
PB = 0.50;
M = 0.30;
K = 0.40;
C = 0;
% Additional parameters definition.
alpha = (K*M - 2) / (M^2 - 4);
beta = (2*K - M) / (M^2 - 4);
gamma = C / (M - 2);
gamma_prime = gamma * (M - 1); % gamma_prime is the gamma' parameter.
% Determine the minimum value for G.
Fmax = max(alpha^2,(3*beta^2+2*alpha*beta));
Lmin = min(LA,LB);
% G = 0.2; % Mind that G is the Gamma parameter defined elsewhere.
G = Fmax  / (Lmin^2);
% Set the fraction of the minimum G to be used.
Gfraction = 0.20; % This value should be set to 1.
G = Gfraction * G;
% G = 0.230816594160215; 

% Set the dimensionality of the search space for the case of beta>=0:
if(beta>=0)
    Dimensionality = FirmsNumber + 1;
else
% Set the dimensionality of the search space for the case of beta<0:
    Dimensionality = 2*FirmsNumber + 1;
end

% Set the number of different initial points to be considered.
N = 2000;
% Set lower and upper bounds for optimization variables.
lb = zeros(1,Dimensionality);
% Set upper bounds for the optimization variables for the case beta < 0:
if(beta>=0)
    ub = [1 1 Inf];
% Set upper bounds for the optimization variables for the case beta > 0:
else
    ub = [1 1 Inf Inf Inf];
end
% Set the display flag parameter to 'iter'.
DisplayFlag = 'iter';
% Set the tolerance value for the minimizer. (Preferable value = 1e-10)
Tolerance = 1e-15;
% Set the Fvals tolerance value for filtering the obtained solutions. (Preferable value = 1e-15)
FvalsTolerance = 1e-15;
% Set the derivative tolerance value for filering the obtained solutions. (Preferable value = 1e-08)
DerivativeTolerance = 1e-08;
% Set the maximum number of iterations to be conducted by the optimizer.
MaxIterations = 2000;
% Set the maximum number of function evaluations to be conducted by the
% optimizer.
MaxFunctionEvaluations = 20000;
% Set the minimum digits accuracy parameter controlling the extraction of
% the representative solution. (Preferable value = 8)
MinimumDigitsAccuracy = 7;

% Run the solver.
[Solutions,Fvals,ExitFlags,Cineq,Ceq] = NonLinearObjectiveFunctionMinimizer(DisplayFlag,Tolerance,MaxIterations,MaxFunctionEvaluations,Dimensionality,N,lb,ub,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Filter Solutions.
[Status,Filtered,Solutions,Fvals,FilterFlag,FD,SD,Fopt,Sopt,Xopt,Popt,Qopt,ExitFlags,Ceq,Cineq] = FilterSolutions(Solutions,ExitFlags,Fvals,Ceq,Cineq,FvalsTolerance,DerivativeTolerance,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);

if(FilterFlag==0)
    % Isolate Solution Points.
    SOL = Solutions(:,1:2);
    % Plot filtered solutions.
    figure('Name','Filtered Solutions');
    plot(SOL(:,1),SOL(:,2),'*m','LineWidth',2.0);
    grid on
    xlabel('T^{\ast}_A');
    ylabel('T^{\ast}_B');
    % Extract Representative Solutions.
    % [RSOL,DigitsAccuracy] = ExtractRepresentativeSolution(SOL,MinimumDigitsAccuracy);
    [RSOL,DigitsAccuracy,~] = ExtractClusterRepresentativeSolutions(Solutions,Fvals,MinimumDigitsAccuracy,FirmsNumber);
    [Sopt,Xopt,Popt,Qopt,Fopt] = RetrieveOptimalModelParameters(RSOL,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
else
    warning('Optimization process failed: %d',FilterFlag);
    RSOL = [];
    DigitsAccuracy = [];
end

% Plot the constraint regions.
PlotConstraintRegions(alpha,beta,LA,LB,PA,PB,RSOL);
% Report final solutions.
ReportOptimalSolutions(RSOL,DigitsAccuracy);

if(FilterFlag==0)
    
    % Perform comparative statics analysis with respect to PA(==px). 
    % This analysis can only be performed when the Nash equilibrium is obtained
    % for the zeroth value of the exogenous parameter, that is PA = 0.
    if(PA==0)
        pxspan = 0:0.01:1;
        T0 = [RSOL(1);RSOL(2)];
        opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
        [px,T] = ode45(@(px,T) OligopolyODEpx(px,T,G,K,LA,LB,M,PB),pxspan, T0);
        figure('Name','Comparative Statics wrt Px(0)');
        plot(px,T(:,1),'-r',px,T(:,2),'-g','LineWidth',2.5)
        grid on
        legend({'Tx(px)','Ty(px)'});
        xlabel('px');
        ylabel('T^{*}(px)');
    end

    % Perform comparative statics analysis with respect to PB(==py).
    % This analysis can only be performed when the Nash equilibrium is obtained
    % for the zeroth value of the exogenous parameter, that is PB = 0.
    if(PB==0)
        pyspan = 0:0.01:1;
        T0 = [RSOL(1);RSOL(2)];
        opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
        [py,T] = ode45(@(py,T) OligopolyODEpy(py,T,G,K,LA,LB,M,PA),pyspan, T0);
        figure('Name','Comparative Statics wrt Py(0)');
        plot(py,T(:,1),'-r',py,T(:,2),'-g','LineWidth',2.5)
        grid on
        legend({'Tx(py)','Ty(py)'});
        xlabel('py');
        ylabel('T^{*}(py)');
    end

    % Perform comparative statics analysis with respect to K.
    % This analysis can only be performed when the Nash equilibrium is obtained
    % for the zeroth value of the exogenous parameter, that is K = 0.
    if(K==0)
        kspan = 0:0.01:(M/2);
        T0 = [RSOL(1);RSOL(2)];
        opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
        [k,T] = ode45(@(K,T) OligopolyODEK(K,T,G,LA,LB,M,PA,PB),kspan, T0);
        figure('Name','Comparative Statics wrt K');
        plot(k,T(:,1),'-r',k,T(:,2),'-g','LineWidth',2.5)
        grid on
        legend({'Tx(K)','Ty(K)'});
        xlabel('K');
        ylabel('T^{*}(K)');
    end

    % Perform comparative statics analysis with respect to M.
    % This analysis can only be performed when the Nash equilibrium is obtained
    % for the zeroth value of the exogenous parameter, that is M = 0.
    if(M==0)
        mspan = (2*K):0.01:1;
        T0 = [RSOL(1);RSOL(2)];
        opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
        [m,T] = ode45(@(M,T) OligopolyODEK(M,T,G,K,LA,LB,PA,PB),mspan, T0);
        figure('Name','Comparative Statics wrt M');
        plot(m,T(:,1),'-r',m,T(:,2),'-g','LineWidth',2.5)
        grid on
        legend({'Tx(M)','Ty(M)'});
        xlabel('M');
        ylabel('T^{*}(M)');
    end

    % Perform comparative statics analysis with respect to Lx(=LA).
    % This analysis can only be performed when the Nash equilibrium is obtained
    % for the zeroth value of the exogenous parameter, that is LA = 0.01.
    if(LA==0.01)
        lxspan = 0.01:0.01:1;
        T0 = [RSOL(1);RSOL(2)];
        opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
        [lx,T] = ode45(@(Lx,T) OligopolyODELx(Lx,T,G,K,LB,M,PA,PB),lxspan, T0);
        figure('Name','Comparative Statics wrt Lx');
        plot(lx,T(:,1),'-r',lx,T(:,2),'-g','LineWidth',2.5)
        grid on
        legend({'Tx(Lx)','Ty(Lx)'});
        xlabel('Lx');
        ylabel('T^{*}(Lx)');
    end
 
    % Perform comparative statics analysis with respect to Lx(=LA).
    % This analysis can only be performed when the Nash equilibrium is obtained
    % for the zeroth value of the exogenous parameter, that is LA = 0.01.
    if(LB==0.01)
        lyspan = 0.01:0.01:1;
        T0 = [RSOL(1);RSOL(2)];
        opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
        [ly,T] = ode45(@(Ly,T) OligopolyODELy(Ly,T,G,K,LA,M,PA,PB),lyspan, T0);
        figure('Name','Comparative Statics wrt Lx');
        plot(ly,T(:,1),'-r',ly,T(:,2),'-g','LineWidth',2.5)
        grid on
        legend({'Tx(Ly)','Ty(Ly)'});
        xlabel('Ly');
        ylabel('T^{*}(Ly)');
    end
end