% This script estimates the optimal influences T1_opt, T2_opt for a set of
% given ranges for the optimal influence model parameters defined for the 
% worst case within a 5 dimmensional grid space.

% Define ranges and corresponding increment step for each parameter of the
% optimal influence model.

P1_MIN = 0.8;
P1_MAX = 0.8;
P1_STEP = 0.1;

P2_MIN = 0.8;
P2_MAX = 0.8;
P2_STEP = 0.01;

THETA_MIN = 0.0;
THETA_MAX = 2.0;
THETA_STEP = 0.1;

DELTA_MIN = 2.0;
DELTA_MAX = 2.5;
DELTA_STEP = 0.001;

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
Topts = [];
Fvals = [];
Params = [];
% Perform the actual grid searching.
for P1 = P1_RANGE
    for P2 = P2_RANGE
        for Theta = THETA_RANGE
            for Delta = DELTA_RANGE
                for Gamma = GAMMA_RANGE
                    params = [P1 P2 Theta Delta Gamma];
                    Params = [Params;params];
                    [T1_opt,T2_opt,Fval] = OptimalInfluences(P1,P2,Theta,Delta,Gamma);
                    T = [T1_opt T2_opt]
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
ParamIndices = [3 4];
% Set parameters' indices and corresponding values for the parameters that
% remain constant.
ConstIndices = [1 2 5];
ConstValues = [0.8 0.8 0.5];
% Perform the actual plotting.
plot_parameter_triples(Topts,Fvals,Params,ParamIndices,ConstIndices,ConstValues);