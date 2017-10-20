% This script estimates the optimal influences T1_opt, T2_opt and associated
% optimal lambda parameters Lambda1_opt,Lambda2_opt for a set of
% given ranges for the optimal influence model parameters defined for the 
% worst case within a 9 dimmensional grid space.

% Define ranges and corresponding increment step for each parameter of the
% optimal influence model.

P0_MIN = 1.0;
P0_MAX = 1.0;
P0_STEP = 0.01;

P1_MIN = 0.2;
P1_MAX = 0.2;
P1_STEP = 0.01;

P2_MIN = 0.2;
P2_MAX = 0.2;
P2_STEP = 0.01;

THETA1_MIN = 0.00;
THETA1_MAX = 1.00;
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

K_MIN = 1.0;
K_MAX = 1.0;
K_STEP = 0.01;

BETA_MIN = 0.45;
BETA_MAX = 0.45;
BETA_STEP = 0.01;

% Set the corresponding ranges for each one of the model parameters.
P0_RANGE = [P0_MIN:P0_STEP:P0_MAX];
P1_RANGE = [P1_MIN:P1_STEP:P1_MAX];
P2_RANGE = [P2_MIN:P2_STEP:P2_MAX];
THETA1_RANGE = [THETA1_MIN:THETA1_STEP:THETA1_MAX];
THETA2_RANGE = [THETA2_MIN:THETA2_STEP:THETA2_MAX];
DELTA_RANGE = [DELTA_MIN:DELTA_STEP:DELTA_MAX];
GAMMA_RANGE = [GAMMA_MIN:GAMMA_STEP:GAMMA_MAX];
K_RANGE = [K_MIN:K_STEP:K_MAX];
BETA_RANGE = [BETA_MIN:BETA_STEP:BETA_MAX];

% Initialize TLopt and Fval vectors for storing best values for each
% parameter tuple configuration. Mind that TL_opt is now a four-dimensional
% optimization vector that takes into consideration parameters
% T1,T2,Lambda1 and Lambda2.
Sopts = [];
TLopts = [];
Fvals = [];
Params = [];
X = [];
Flags1 = [];
Flags2 = [];
% Perform the actual grid searching. Mind that now there exist only 9
% dimensions to be traversed.
for K = K_RANGE
    for P0 = P0_RANGE
        for P1 = P1_RANGE
            for P2 = P2_RANGE
                for Theta1 = THETA1_RANGE
                    for Theta2 = THETA2_RANGE
                        for Delta = DELTA_RANGE
                            for Gamma = GAMMA_RANGE
                                for Beta = BETA_RANGE
                                    T1_MIN = max(0,Theta1+K-1);
                                    T1_MAX = min(1-Theta2,K);
                                    T2_MIN = max(0,Theta2+K-1);
                                    T2_MAX = min(1-Theta1,K);
                                    if(and((T1_MIN<T1_MAX),(T2_MIN<T2_MAX)))
                                        params = [P0,P1,P2,Theta1,Theta2,Delta,Gamma,K,Beta];
                                        Params = [Params;params];
                                        [T1_opt,T2_opt,Lambda1_opt,Lambda2_opt,Fval,Flag1,Flag2] = GeneralOptimalInfluencesX(P0,P1,P2,Theta1,Theta2,Delta,Gamma,K,Beta);
                                        [S0_opt,S1_opt,S2_opt] = GeneralSoptimalX(T1_opt,T2_opt,Lambda1_opt,Lambda2_opt,Theta1,Theta2);
                                        x = consensus(P0,P1,P2,S0_opt,S1_opt,S2_opt);
                                        TL = [T1_opt T2_opt,Lambda1_opt,Lambda2_opt]
                                        S = [S0_opt,S1_opt,S2_opt];
                                        Sopts = [Sopts;S];
                                        X = [X;x];
                                        TLopts = [TLopts;TL];
                                        Fvals = [Fvals;Fval];
                                        Flags1 = [Flags1;Flag1];
                                        Flags2 = [Flags2;Flag2];
                                    end;
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
ParamIndices = [4 5];
% Set parameters' indices and corresponding values for the parameters that
% remain constant.
ConstIndices = setdiff([1:9],ParamIndices);
ConstValues = [1.0 0.2 0.2 2.1 0.59 1.0 0.45];
% Do the actual plotting.
plot_parameters_tuplesX(TLopts,Fvals,Sopts,X,Params,ParamIndices,ConstIndices,ConstValues);

ParamSubRanges = [0.1 0.25 0.4;0.1 0.25 0.4];
slice_plot_parameters_tuplesX(TLopts,Fvals,Sopts,X,Params,ParamIndices,ParamSubRanges,ConstIndices,ConstValues);