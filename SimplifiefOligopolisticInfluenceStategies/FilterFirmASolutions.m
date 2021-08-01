function [Solutions,FilterFlag,SD,Fopt,Sopt,Xopt,Popt,Qopt] = FilterFirmASolutions(Solutions,TB,ExitFlags,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime)

% This function filters out the true best response solutions for Firm A
% given the investement level selected by Firm B.
% -------------------------------------------------------------------------
% Solutions: is a column vector storing the optimal direct influence level
%            exerted by Firm A (TAopt(TB)) as a response to the direct
%            influence exerted by Firm B. Thus, vector Solutions stores the
%            outcome of the non-linear optimization process for Firm A.
% TB: is a scalar value corresponding to the direct influence exerted by
%     Firm B.
% -------------------------------------------------------------------------
% The obtained solutions will be filtered according to the following
% criteria:
% (i):   ExitFlag >= 0, indicating that the optimization process was
%                       terminated normally.
% (ii):  DDa < 0, indicating that the equilibrium point for the first firm
%                 satisfies the Second Order Conditions.
% (iii): FAopt >= 0, indicating that the equilibrium point for the first
%                    firm incurs a non-negative profit value.
% (iv):  PAopt >= 0, indicating that the equilibrium point for the first 
%                    firm incurs a non-negative price value.
% (v):   QAopt >= 0, indicating that the equilibrium point for the first 
%                    firm incurs a non-negative quantity value.
% -------------------------------------------------------------------------
% FilterFlag is the output variable indicating the exact stage during the 
% filtering process where no valid solutions could be attained.
% -------------------------------------------------------------------------
% FilterFlag =   0 ==>  Optimization process terminated normally.
% FilterFlag =  -1 ==>  Optimization process did not reach any acceptable 
%                       solutions since the associated ExitFlags were all 
%                       negative: [(i): ExitFlags > 0]
% FilterFlag =  -2 ==>  Optimization process yields solutions that do not
%                       satisfy the Second Order Convexity Conditions for
%                       the first firm: [(ii): DDa > 0]
% FilterFlag =  -3 ==>  Optimization process yields solutions that are
%                       associated with negative profit values for the first
%                       firm: [(iii): FAopt < 0]
% FilterFlag =  -4 ==>  Optimization process yields solutions that are
%                       associated with negative price values for the first
%                       firm: [(iv): PAopt < 0]
% FilterFlag = -5 ==>   Optimization process yields solutions that are
%                       associated with negative quantity values for the 
%                       first firm: [(v): QAopt < 0]
% -------------------------------------------------------------------------
% IMPORTANT NOTE:
% 
% Notice that there is no need in testing the validity of the First Order
% Conditions for the obtained solutions for Firm A since they are
% internally handled by Matlab. Also take into consideration that
% optimality conditions are tested only for the solution points obtained
% for Firm A since they are the result of an underlying optimization
% process.
% -------------------------------------------------------------------------


% FILTERING PROCESS

% Initialize the output variable.
FilterFlag = 0;
% Initialize a variable storing the number of criterions that need to be
% examined.
CriterionsNumber = 5;
% Get the number of obtained solutions.
SolutionsNumber = size(Solutions,1);
% Initialize a logical matrix storing the evaluation status of each
% solution regarding a particular condition.
I = zeros(SolutionsNumber,CriterionsNumber);

% Set the optimal solutions vector for the investment level of Firm A.
TAopt = Solutions;
% Construct a vector of the same dimensionality as TAopt storing the
% investement level of Firm B.
TBo = repmat(TB,size(TAopt,1),1);
% Merge the optimal solution points for Firm A with the ad hoc value for
% the direct influnce for Firm B.
Topts = [TAopt,TBo];
% Obtain second order optimality conditions for the optimal direct influence
% levels of Firm A.
SD = FirmAOptimalityConditions(TAopt,TB,C,G,LA,LB,PA,PB,alpha,beta,gamma);
% Compute fundamental model quantities for each solution point.
[Sopt,Xopt,Popt,Qopt,Fopt] = RetrieveOptimalModelParameters(Topts,...
                             C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
                         
% Evaluate the status of the obtained solutions regarding condition (i):
I(:,1) = (ExitFlags >= 0);
% Evalute the status of the obtained solutions regarding condition (ii):
I(:,2) = (SD < 0);
% Evalute the status of the obtained solutions regarding condition (iii):
I(:,3) = (Fopt(:,1) >= 0);
% Evalute the status of the obtained solutions regarding condition (iv):
I(:,4) = (Qopt(:,1) >= 0);
% Evalute the status of the obtained solutions regarding condition (v):
I(:,5) = (Popt(:,1) >= 0);

% Perform the column-wise conjunction of the boolean values stored in
% matrix I and store the outcome in the boolean vector Io.
Io = ones(SolutionsNumber,1);
for k = 1:CriterionsNumber
    Io = and(Io,I(:,k));
    if(sum(Io)==0)
        FilterFlag = -k;
        break
    end
end

% Mind that the logical ones of the boolean vector Io identify the
% solutions that satisfy all the required conditions. Therefore, the next
% step is to acquire all the output variables according to the valid
% solution indices stored in the logical vector Io.
Solutions = Solutions(Io);
SD = SD(Io,:);
Fopt = Fopt(Io,:);
Sopt = Sopt(Io,:);
Xopt = Xopt(Io,:);
Popt = Popt(Io,:);
Qopt = Qopt(Io,:);
end