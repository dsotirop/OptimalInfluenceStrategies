% This script estimates the optimal influences T1_opt, T2_opt for a set of
% given ranges for the optimal influence model parameters defined for the 
% worst case within a 9 dimmensional grid space.

% Define ranges and corresponding increment step for each parameter of the
% optimal influence model.

P0_MIN = 1.0;
P0_MAX = 1.0;
P0_STEP = 0.01;

P1_MIN = 0.1;
P1_MAX = 0.1;
P1_STEP = 0.01;

P2_MIN = 0.0;
P2_MAX = 1.0;
P2_STEP = 0.01;

LAMBDA1_MIN = 0.25;
LAMBDA1_MAX = 0.25;
LAMBDA1_STEP = 0.01;

LAMBDA2_MIN = 0.25;
LAMBDA2_MAX = 0.25;
LAMBDA2_STEP = 0.01;

THETA1_MIN = 0.2;
THETA1_MAX = 0.2;
THETA1_STEP = 0.01;

THETA2_MIN = 0.0;
THETA2_MAX = 1.0;
THETA2_STEP = 0.01;

DELTA_MIN = 2.1;
DELTA_MAX = 2.1;
DELTA_STEP = 0.01;

GAMMA_MIN = 0.59;
GAMMA_MAX = 0.59;
GAMMA_STEP = 0.01;

% Set the corresponding ranges for each one of the model parameters.
P0_RANGE = [P0_MIN:P0_STEP:P0_MAX];
P1_RANGE = [P1_MIN:P1_STEP:P1_MAX];
P2_RANGE = [P2_MIN:P2_STEP:P2_MAX];
LAMBDA1_RANGE = [LAMBDA1_MIN:LAMBDA1_STEP:LAMBDA1_MAX];
LAMBDA2_RANGE = [LAMBDA2_MIN:LAMBDA2_STEP:LAMBDA2_MAX];
THETA1_RANGE = [THETA1_MIN:THETA1_STEP:THETA1_MAX];
THETA2_RANGE = [THETA2_MIN:THETA2_STEP:THETA2_MAX];
DELTA_RANGE = [DELTA_MIN:DELTA_STEP:DELTA_MAX];
GAMMA_RANGE = [GAMMA_MIN:GAMMA_STEP:GAMMA_MAX];

% Initialize T_opt and Fval vectors for storing best values for each
% parameter tuple configuration.
Sopts = [];
Topts = [];
Fvals = [];
Params = [];
X = [];
Frevenue = [];
Fcost = [];
% Perform the actual grid searching.
for P0 = P0_RANGE
    for P1 = P1_RANGE
        %P2 = P1;
        for P2 = P2_RANGE
            for Lambda1 = LAMBDA1_RANGE
                for Lambda2 = LAMBDA2_RANGE
                %Lambda2 = Lambda1;
                    for Theta2 = THETA2_RANGE
                        for Theta1 = THETA1_RANGE
                            for Delta = DELTA_RANGE
                                for Gamma = GAMMA_RANGE
                                    % if block should be commented out when
                                    % enforcing the T1 + T2 = 1 equality
                                    % constraint.
                                    %if(Theta1+Theta2<1)
                                        params = [P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma];
                                        Params = [Params;params];
                                        %[Theta1,Theta2]
                                        [T1_opt,T2_opt,Fval] = GeneralOptimalInfluences(P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma);
                                        [S0_opt,S1_opt,S2_opt] = GeneralSoptimal(T1_opt,T2_opt,Lambda1,Lambda2,Theta1,Theta2);
                                        x = consensus(P0,P1,P2,S0_opt,S1_opt,S2_opt);
                                        T = [T1_opt T2_opt]
                                        [Fr,Fc] = GeneralPartialObjectiveFunction(T,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma);
                                        S = [S0_opt,S1_opt,S2_opt];
                                        Sopts = [Sopts;S];
                                        X = [X;x];
                                        Topts = [Topts;T];
                                        Fvals = [Fvals;Fval];
                                        Frevenue = [Frevenue;Fr];
                                        Fcost = [Fcost;Fc];
                                    %end;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

% Compute corresponding maximum profit by performing sign reversion.
Fvals = -Fvals;

% Set the parameters' indices that are allowed to vary.
ParamIndices = [3 7];
% Set parameters' indices and corresponding values for the parameters that
% remain constant.
ConstIndices = setdiff([1:9],ParamIndices);
ConstValues = [1.0 0.1 0.25 0.25 0.2 2.1 0.59];
% Do the actual plotting.
plot_parameters_tuples(Topts,Fvals,Frevenue,Fcost,Sopts,X,Params,ParamIndices,ConstIndices,ConstValues);

ParamSubRanges = [0.18 0.26 0.34;0.15 0.30 0.40];
slice_plot_parameters_tuples(Topts,Fvals,Sopts,X,Params,ParamIndices,ParamSubRanges,ConstIndices,ConstValues);