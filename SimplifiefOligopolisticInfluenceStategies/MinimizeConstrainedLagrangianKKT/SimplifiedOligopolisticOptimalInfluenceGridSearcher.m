% This script file computes the optimal investment levels TAopt and TBopt
% for the two firms (Firm A and Firm B) of the Simplified Oligopolistic 
% Optimal Influence model which also involves a single consumer C. The rest
% of the model parameters (Sopt,Xopt,Popt,Qopt,Fopt) are also computed.

% Specifically, the internal model parameters whose optimal values are to
% be determined during the grid searching optimization process are the
% follwing:
% (i):   [TAopt,TBopt]
% (ii):  [SAopt,SCopt,SBopt]
% (iii): [XAopt,XBopt]
% (iv):  [PAopt,PBopt] (These are the optimal prices!!!)
% (v):   [QAopt,QBopt]
% (vi):  [FAopt FBopt FA_rev_opt FB_rev_opt FA_cost_opt FB_rev_opt]

% The grid searching process will be conducted within the 8-dimensional
% space defined by the Cartesian product of the following external model
% parameters:
% (i):    LA (direct influnce exerted by consumer C on Firm A).
% (ii):   LB (direct influnce exerted by consumer C on Firm B).
% (iii):  PA (initial belief consumer C holds for product A).
% (iv):   PB (initial belief consumer C holds for product B).
% (v):    M  (sensitivity coefficient).
% (vi):   K  (sensitivity coefficient).
% (vii):  C  (marginal cost).
% (viii): G  (marginal influence cost or Gamma).

% IMPORTANT NOTE!!!
% Mind that the underlying continuous game may not have an equlibrium point
% for any given configuration of the external parameters. Since all the 
% internal parameters to be determined accept positive values, the value of 
% (-1) will be utilized in order to indicate the absence of an equilibrium
% point for a particular configuration of the external parameters.

% Define ranges and corresponding increment step for each external parameter 
% of the Simplified Oligopolistic Optimal Influence model.

clc
clear

% Construct a cell array storing the names of model parameters.
ParamsNames = {'LA','LB','PA','PB','M','K','C','G'};
% Get the number of different parameters.
ParamsNum = length(ParamsNames);


% Parameter #1
LA_MIN = 0.50;
LA_MAX = 0.50;
LA_STEP = 0.01;
% Parameter #2
LB_MIN = 0.50;
LB_MAX = 0.50;
LB_STEP = 0.01;
% Parameter #3
PA_MIN = 0.45;
PA_MAX = 0.45;
PA_STEP = 0.01;
% Parameter #4
PB_MIN = 0.40;
PB_MAX = 0.40;
PB_STEP = 0.01;
% Parameter #5
M_MIN = 0.50;
M_MAX = 0.50;
M_STEP = 0.01;
% Parameter #6
K_MIN = 0.00;
K_MAX = 1.00;
K_STEP = 0.01;
% Parameter #7
C_MIN = 0.0;
C_MAX = 0.0;
C_STEP = 0.01;
% Parameter #8 (G or Gamma!!!)
% Mind that G is the Gamma parameter defined elsewhere.
% Determine the minimum value for G.
alpha0 = (K_MAX*M_MAX - 2) / (M_MIN^2 - 4);
beta0 = (2*K_MAX - M_MIN) / (M_MIN^2 - 4);
Fmax = max(alpha0^2,(3*beta0^2+2*alpha0*beta0));
Lmin = min(LA_MIN,LB_MIN);
G0 = Fmax  / (Lmin^2);
% When Gfraction is set to 1.00 it represents an overestimated minimum
% value for the Gamma parameter controlling the convexity of the profit
% functions for both firms. Thus, it is reasonable to experiment with 
% significanlty lower values for the Gamma parameter. 
Gfraction = 1.50;  
G0 = Gfraction * G0;
G_MIN = G0;
G_MAX = G0;
G_STEP = 0.01;

% Set the corresponding ranges for each parameter of the model.
LA_RANGE = LA_MIN:LA_STEP:LA_MAX;
LB_RANGE = LB_MIN:LB_STEP:LB_MAX;
PA_RANGE = PA_MIN:PA_STEP:PA_MAX;
PB_RANGE = PB_MIN:PB_STEP:PB_MAX;
M_RANGE = M_MIN:M_STEP:M_MAX;
K_RANGE = K_MIN:K_STEP:K_MAX;
C_RANGE = C_MIN:C_STEP:C_MAX;
G_RANGE = G_MIN:G_STEP:G_MAX;

% Initialize a cell array storing the range values for each exogenous
% parameter of the model.
PARAMETERS = cell(1,ParamsNum);
% Populate the cell array with the corresponding parameter ranges.
PARAMETERS{1} = LA_RANGE;
PARAMETERS{2} = LB_RANGE;
PARAMETERS{3} = PA_RANGE;
PARAMETERS{4} = PB_RANGE;
PARAMETERS{5} = M_RANGE;
PARAMETERS{6} = K_RANGE;
PARAMETERS{7} = C_RANGE;
PARAMETERS{8} = G_RANGE;
% Get the range length for each parameter stored in the PARAMETERS cell 
% array.
RANGE_SIZES = cellfun(@length,PARAMETERS);
% Check the number of parameters that were allowed to vary during the grid
% searching process.
if((sum(RANGE_SIZES>1))>1)
    % If the exists more than one parameter with range length greater than 1
    % then throw an error.
    error('More than one parameters where allowed to vary during the grid searching process!')
else
    % Determine the unique parameter index which is allowed to vary. 
    ParamIndex = vec2ind(double(RANGE_SIZES>1)');
end

% Get the number of points to be evaluated during the grid searching
% process.
RangeLength = RANGE_SIZES(ParamIndex);

% Initialize matrices storing the fundamental intermediate quantities of
% the underlying optimization problem within the simplified oligopolistic 
% enviroment as well as the values of the external optimization parameters 
% for each step of the multi-dimensional grid searching process.
Topts = zeros(RangeLength,2); % Optimal investment levels vectors [TA,TB].
Sopts = zeros(RangeLength,3); % Optimal limiting influences vectors [SA,SC,SB].
Xopts = zeros(RangeLength,2); % Optimal limiting beliefs vectors [XA,XB].
Popts = zeros(RangeLength,2); % Optimal prices vectors [PA,PB].
Qopts = zeros(RangeLength,2); % Optimal quantities vectors [QA,QB].
Fopts = zeros(RangeLength,6); % Optimal profits vectos [FAopt,FBopt,FA_rev_opt,FB_rev_opt,FA_cost_opt,FB_cost_opt].
FilterFlags = zeros(RangeLength,1); % Solution filtering flag which may be indicative of a non-existing solution. 
DigitAccuracies = zeros(RangeLength,1); % Length of maximal sequence of identical digits within the obtained optimal solutions.

% Initialize the matrix storing the 8-tuples of the varying optimization
% parameters.
Params = zeros(RangeLength,ParamsNum);
% Populate matrix Params.
for param_index = 1:ParamsNum
    if(param_index==ParamIndex)
        Params(:,param_index) = PARAMETERS{param_index}';
    else
        Params(:,param_index) = repmat(PARAMETERS{param_index},RangeLength,1);
    end
end

% Initialize internal solver parameters.
% Set the number of initial solution points.
N = 2000;
% Set the tolerance value for the minimizer. (Preferable value = 1e-10)
Tolerance = 1e-15;
% Set the Fvals tolerance value for filtering the obtained solutions. (Preferable value = 1e-15)
FvalsTolerance = 1e-15;
% Set the derivative tolerance value for filering the obtained solutions. (Preferable value = 1e-08).
DerivativeTolerance = 1e-08;
% Set the maximum number of iterations to be conducted by the optimizer.
MaxIterations = 2000;
% Set the maximum number of function evaluations to be conducted by the
% optimizer.
MaxFunctionEvaluations = 20000;
% Set the display flag parameter to 'off' or to 'iter'.
DisplayFlag = 'off';
% Set the minimum digits accuracy.
MinimumDigitsAccuracy = 7;

fprintf('Grid Evaluation Process in Progess...\n');
% Perform the actual grid searching in parallel by processing matrix Params
% in a row-wise manner.
for grid_index = 1:1:RangeLength
   % Acquire the current value for each parameter pertaining to the grid
   % searching process.
   LA = Params(grid_index,1);
   LB = Params(grid_index,2);
   PA = Params(grid_index,3);
   PB = Params(grid_index,4);
   M = Params(grid_index,5);
   K = Params(grid_index,6);
   C = Params(grid_index,7);
   G = Params(grid_index,8);
   alpha = (K*M - 2) / (M^2 - 4);
   beta = (2*K - M) / (M^2 - 4);
   gamma = C / (M - 2);
   gamma_prime = gamma * (M - 1); % gamma_prime is the gamma' parameter.
   [TAopt,TBopt,FilterFlag,DigitsAccuracy] = SimplifiedOligopolisticOptimalInfluences(N,DisplayFlag,Tolerance,FvalsTolerance,DerivativeTolerance,MaxIterations,MaxFunctionEvaluations,MinimumDigitsAccuracy,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
   DigitsAccuracy = min(DigitsAccuracy);
   T = [TAopt TBopt];
   if(FilterFlag==0)
       [S,X,P,Q,F] = RetrieveOptimalModelParameters(T,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
   else
       S = -1*ones(1,3);
       X = -1*ones(1,2);
       P = -1*ones(1,2);
       Q = -1*ones(1,2);
       F = -1*ones(1,6);
   end
   % Print current solution.
   param_value_string = strcat([ParamsNames{ParamIndex} ' = ' num2str(Params(grid_index,ParamIndex),10)]);
   fprintf('%s TAopt = %f TBopt = %f FilterFlag = %d DigitsAccuracy = %d\n',param_value_string,TAopt,TBopt,FilterFlag,DigitsAccuracy);
   Topts(grid_index,:) = T;
   FilterFlags(grid_index) = FilterFlag;
   DigitAccuracies(grid_index) = DigitsAccuracy;
   Sopts(grid_index,:) = S;
   Xopts(grid_index,:) = X;
   Popts(grid_index,:) = P;
   Qopts(grid_index,:) = Q;
   Fopts(grid_index,:) = F;
end

% Set parameters' indices and corresponding values for the parameters that
% remain constant.
ConstIndices = setdiff(1:length(params),ParamIndex);
ConstValues = params(ConstIndices);
% Perform plotting operations.
PlotParametersTuples(Topts,Sopts,Xopts,Popts,Qopts,Fopts,Params,ParamIndex,ConstIndices,ConstValues)