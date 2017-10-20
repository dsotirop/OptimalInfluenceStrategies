% This script file investigates the range of the counter intuitive
% influence regions DeltaTheta = (Theta_hat - Theta_star) as a function of 
% the DeltaP = (P2 - P1) difference.

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

THETA1_MIN = 0.5;
THETA1_MAX = 0.5;
THETA1_STEP = 0.01;

THETA2_MIN = 0.0;
THETA2_MAX = 1.0;
THETA2_STEP = 0.01;

DELTA_MIN = 1.94;
DELTA_MAX = 1.94;
DELTA_STEP = 0.01;

GAMMA_MIN = 1.0;
GAMMA_MAX = 1.0;
GAMMA_STEP = 0.5;

K_MIN = 1.0;
K_MAX = 1.0;
K_STEP = 0.01;

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
K_RANGE = [K_MIN:K_STEP:K_MAX];

% Initialize T_opt and Fval vectors for storing best values for each
% parameter tuple configuration.
Sopts = [];
Topts = [];
Fvals = [];
Params = [];
X = [];
Flags1 = [];
Flags2 = [];
% Perform the actual grid searching.
for K = K_RANGE
    for P0 = P0_RANGE
        for P1 = P1_RANGE
            for P2 = P2_RANGE
                for Lambda1 = LAMBDA1_RANGE
                    for Lambda2 = LAMBDA2_RANGE
                        for Theta2 = THETA2_RANGE
                            for Theta1 = THETA1_RANGE
                                for Delta = DELTA_RANGE
                                    for Gamma = GAMMA_RANGE
                                    % if block should be commented out when
                                    % enforcing the T1 + T2 = K equality
                                    % constraint.
                                    T1_MIN = max(0,Theta1+K-1);
                                    T1_MAX = min(1-Theta2,K);
                                    T2_MIN = max(0,Theta2+K-1);
                                    T2_MAX = min(1-Theta1,K);
                                    if(and((T1_MIN<T1_MAX),(T2_MIN<T2_MAX)))
                                        params = [P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma,K];
                                        Params = [Params;params];
                                        [T1_opt,T2_opt,Fval,Flag1,Flag2] = GeneralOptimalInfluences(P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma,K);
                                        [S0_opt,S1_opt,S2_opt] = GeneralSoptimal(T1_opt,T2_opt,Lambda1,Lambda2,Theta1,Theta2);
                                        x = consensus(P0,P1,P2,S0_opt,S1_opt,S2_opt);
                                        T = [T1_opt T2_opt]
                                        S = [S0_opt,S1_opt,S2_opt];
                                        Sopts = [Sopts;S];
                                        X = [X;x];
                                        Topts = [Topts;T];
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
end;
           
% Compute corresponding maximum profit by performing sign reversion.
Fvals = -Fvals;

% Set the varying parameters indices storing the theta parameter last.
ParamIndices = [3 7];
% Set the constant paramters indices.
ConstIndices = setdiff([1:10],ParamIndices);
% Set the constant parameters values.
ConstValues = [1 0.1 0.25 0.25 0.5 1.94 1 1];

% Set the name of the common constant parameter name.
% For the moment the number of varying parameters must be exactly 2 since 
% there is no implementation for a number of 3 of varying parameters. So
% the value of the variable CommonParameterName must be set to the empty
% string as it appears below.
CommonParameterName = '';

% Plot influence regions.
plot_influence_regions(Topts,Sopts,Flags1,Flags2,Params,ParamIndices,ConstIndices,ConstValues,CommonParameterName);