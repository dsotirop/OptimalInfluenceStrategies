function [FD,SD] = OptimalityConditionsTester(Solutions,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function evaluates the First and Second Order optimality conditions   
% for each pair of optimal solutions (TAopt,TBopt). Thus, for each solution  
% pair (TAopt,TBopt) which is stored in a row-wise manner within matrix Solutions 
% Solutions there exists a corresponding row-vector entry in matrices FD and 
% SD such that:
%               FD(k,:) = [Da(TAopt,TBopt),Db(TAopt,TBopt)] and  
%               SD(k,:) = [DDa(TAopt,TBopt),Db(TAopt,TBopt)].
% Mind that: 
%           Da: is the first order derivative of fa with respect to TA
%           Db: is the first order derivative of fb with respect to TB
%          DDa: is the second order derivative of fa with respect to TA 
%          DDb: is the second order derivative of fb with respect to TB

% Isolate the optimal investment levels for each firm: TAopt and TBopt.
TAopt = Solutions(:,1);
TBopt = Solutions(:,2);
% Evaluate First Order Derivatives.
Da = FirmAProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
Db = FirmBProfitFirstDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
% Evaluate Second Order Derivatives.
DDa = FirmAProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);
DDb = FirmBProfitSecondDerivative(C,G,LA,LB,PA,PB,alpha,beta,gamma,TAopt,TBopt);

% Construct the output variables FD and SD which will be storing the first 
% and second order derivatives for both objective functions for each solution 
% pair. Therefore, FD ans SD should be a [N x 2] matrices.
FD = [Da,Db];
SD = [DDa,DDb];

end