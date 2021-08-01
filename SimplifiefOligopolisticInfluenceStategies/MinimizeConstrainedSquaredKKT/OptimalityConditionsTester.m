function [FD,SD] = OptimalityConditionsTester(Solutions,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function evaluates the First and Second Order optimality conditions   
% for each pair of optimal solutions (TAopt,TBopt). Thus, for each solution  
% pair (TAopt,TBopt) which is stored in a row-wise manner within matrix Solutions 
% Solutions there exists a corresponding row-vector entry in matrices FD and 
% SD such that:
%               FD(k,:) = [Da(TAopt,TBopt),Db(TAopt,TBopt)] and  
%               SD(k,:) = [DDa(TAopt,TBopt),Db(TAopt,TBopt)].
% Mind that Da: is the first order derivative of fa with respect to TA, 
%           Db: is the first order derivative of fb with respect to TB,
%          DDa: is the second order derivative of fa with respect to TA and
%          DDb: is the second order derivative of fb with respect to TB.

% Get the number of solutions to be tested for first order optimality.
N = size(Solutions,1);

% Initialize the output variables FD and SD which will be storing the first 
% and second order derivatives for both objective functions for each solution 
% pair. Therefore, FD ans SD should be a [N x 2] matrices.
FD = zeros(N,2);
SD = zeros(N,2);

% Loop through the various solutions in order to compute FOCs ans SOCs for 
% both firms for each solution pair.
for k = 1:1:N
    % Get the obtained optimal solutions per variable.
    TAopt = Solutions(k,1);
    TBopt = Solutions(k,2);
    % Evaluate first-order derivatives.
    Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
    Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
    % Set the first-order derivatives vector.
    FD(k,:) = [Da Db];
    % Evaluate second-order derivatives.
    DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
    DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
    % Set the second-order derivatives vector.
    SD(k,:) = [DDa DDb];
end;

end

