% This script estimates the optimal influences T1_opt, T2_opt for a set of
% given ranges for the optimal influence model parameters defined for the 
% worst case within a 5 dimmensional grid space.

% Define ranges and corresponding increment step for each parameter of the
% optimal influence model.

P0 = 1;

P1_MIN = 0.1;
P1_MAX = 0.1;
P1_STEP = 0.01;

P2_MIN = 0.6;
P2_MAX = 0.6;
P2_STEP = 0.0005;

THETA_MIN = 0.0;
THETA_MAX = 2.0;
THETA_STEP = 0.01;

DELTA_MIN = 1.5;
DELTA_MAX = 1.5;
DELTA_STEP = 0.05;

GAMMA_MIN = 0.5;
GAMMA_MAX = 0.5;
GAMMA_STEP = 0.5;

P1_RANGE = [P1_MIN:P1_STEP:P1_MAX];
P2_RANGE = [P2_MIN:P2_STEP:P2_MAX];
THETA_RANGE = [THETA_MIN:THETA_STEP:THETA_MAX];
DELTA_RANGE = [DELTA_MIN:DELTA_STEP:DELTA_MAX];
GAMMA_RANGE = [GAMMA_MIN:GAMMA_STEP:GAMMA_MAX];

% Initialize T_opt and Fval vectors for storing best values for each
% parameter tuple configuration.
Sopts = [];
Topts = [];
Fvals = [];
Params = [];
X = [];
% Perform the actual grid searching.
for P1 = P1_RANGE
    %P2 = P1;
    for P2 = P2_RANGE
        for Theta = THETA_RANGE
            for Delta = DELTA_RANGE
                for Gamma = GAMMA_RANGE
                    params = [P1 P2 Theta Delta Gamma];
                    Params = [Params;params];
                    [T1_opt,T2_opt,Fval] = OptimalInfluences(P1,P2,Theta,Delta,Gamma);
                    [S0_opt,S1_opt,S2_opt] = Soptimal(T1_opt,T2_opt,Theta); 
                    x = consensus(P0,P1,P2,S0_opt,S1_opt,S2_opt);
                    T = [T1_opt T2_opt]
                    S = [S0_opt,S1_opt,S2_opt];
                    X = [X;x];
                    Sopts = [Sopts;S];
                    Topts = [Topts;T];
                    Fvals = [Fvals;Fval];
                end;
            end;
        end;
    end;
end;

% Compute corresponding maximum profit by performing sign reversion.
Fvals = -Fvals;

% Set the parameters' indices that are allowed to vary.
ParamIndex = 3;
% Set parameters' indices and corresponding values for the parameters that
% remain constant.
ConstIndices = [1 2 4 5];
ConstValues = [0.1 0.6 1.5 0.5];
% Perform the actual plotting.
plot_parameter_quadruples(Topts,Fvals,Params,ParamIndex,ConstIndices,ConstValues)
% Plot optimal influences with respect to the parameter that is allowed to
% vary.
figure('Name','Optimal Influences');
S0optimals = Sopts(:,1);
S1optimals = Sopts(:,2);
S2optimals = Sopts(:,3);
Thetas = Params(:,3);
hold on
plot(Thetas,S0optimals,'*b');
plot(Thetas,S1optimals,'*r');
plot(Thetas,S2optimals,'*g');
xlabel('Theta');
ylabel('SOopt/S1opt/S2opt');
grid on
figure('Name','Optimal Influences Differences');
hold on
plot(Thetas,S0optimals,'*b');
plot(Thetas,S1optimals+S1optimals,'*k');
xlabel('Theta');
ylabel('SOopt/S1opt+S2opt');
hold off
grid on
figure('Name','Limit Beliefs');
plot(Thetas,X,'*m');
xlabel('Theta');
ylabel('Limit Belief');
grid on