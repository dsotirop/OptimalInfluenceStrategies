function [Status,Filtered,Solutions,Fvals,FilterFlag,FD,SD,Fopt,Sopt,Xopt,Popt,Qopt,ExitFlags,Ceq,Cineq] = FilterSolutions(Solutions,ExitFlags,Fvals,Ceq,Cineq,FvalsTolerance,DerivativeTolerance,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime)

% This function filters out the true equibrium points for a given instance 
% of the two-player continuous game that underlies the Simplified 
% Olipolistic Influence model between the two firms (FirmA and FirmB) and 
% the single consumer (Consumer C).

% FvalsTolerance is the maximum acceptable value for the optimization
% parameter Fvals. For this particular problem, the minimization of the
% squared expresssions of the first order KKT conditions, Fvals corresponds
% to the acquired minimum value. Therefore, it should be as small as
% possible. [FvalsTolerance <= 1e-15]

% DerivativeTolerance is the parameter value indicating whether the first
% derivative is equal to zero. Specifically, we assume that:
% |FD| <= DerivativeTolerance ==> FD = 0.
% [DerivativeTolerance <= 1e-10]
% ----------------------------------------------------------------------------------
% CASE I: (beta >= 0)
% ----------------------------------------------------------------------------------
% In fact, First Order Conditions will be evaluated within the Lagrangian
% context according to which we have that:
%
% La = fa + MUa * TA + MUab * (1-TA-TB)
% Lb = fb + MUb * TB + MUab * (1-TB-TA)
% 
% Therefore, the actual formulation for the First Order Conditions will be
% given by the following equations:
%
% DLa = Da + MUa - MUab = 0
% DLb = Db + MUb - MUab = 0
% ----------------------------------------------------------------------------------
% CASE II: (beta < 0)
% ----------------------------------------------------------------------------------
% When First Order Conditions are evaluated within the Lagrangian context
% for the case of negative beta we have that:
% 
% La = fa + MUa * TA + MUab1 * (1 - TA - TB) = 0  
% + MUab2 * (alpha * LB * TA + beta * LA * TB + LA * LB * (alpha * PA + beta * PB))
% + MUab3 * (beta * LB * TA + alpha * LA * TB + LA * LB * (alpha * PB + beta * PA))
%
% Therefore, the actual formulation for the First Order Conditions will be 
% given by the following equations:
% 
% DLa = Da + MUa - Mab1 + alpha * LB * MUab2 + beta * LB * MUab3 = 0
% DLb = Db + MUb - Mab1 + beta * LA * MUab2 + alpha * LA * MUab3 = 0
% ----------------------------------------------------------------------------------
% CASE I: (beta >= 0)
% ----------------------------------------------------------------------------------
% By letting: 
%            DeltaMUa = MUab - MUa
%            DeltaMUb = MUab - Mub
%
% First Order Conditions may be re-expressed as:
%
% Da - DeltaMUa = 0 ==> |Da - DeltaMUa| <= DerivativeTolerance.
% Db - DeltaMUb = 0 ==> |Db - DeltaMUb| <= DerivativeTolerance.
% ----------------------------------------------------------------------------------
% CASE II : (beta < 0)
% ----------------------------------------------------------------------------------
% By letting:
%            DeltaMUa = (MUab1 - alpha*LB*MUab2 - beta*LB*MUab3) - MUa
%            DeltaMUb = (MUab1 - beta*LA*MUab2 - alpha*LA*MUab3) - MUb
% ----------------------------------------------------------------------------------


% The obtained solutions will be filtered according to the following
% criteria:
% (i):  ExitFlag >= 0,indicating that the optimization process was
%                     terminated normally.
% (ii): Fval <= FvalsTolerance, indicating that the optimization process
%                               actually achieved a minimum since the
%                               underlying optimization problem is a
%                               minimization problem.
% (iii):|Da - DeltaMUa| <= DerivativeTolerance, 
%                                    indicating that the equlibrium point
%                                    for the first firm satisfies the  
%                                    First Order Conditions.  
% (iv): |Db - DeltaMUb| <= DerivativeTolerance,
%                                    indicating that the equilibrium point
%                                    for the second firm satisfies the
%                                    First Order Conditions.
% (v): DDa < 0, indicating that the equilibrium point for the first firm
%               satisfies the Second Order Conditions.
% (vi):DDb < 0, indicating that the equilibrium point for the second firm 
%               satisfies the Second Order Conditions.
% (vii):  FAopt >= 0, indicating that the equilibrium point for the first
%                     firm incurs a non-negative profit value.
% (viii): FBopt >= 0, indicating that the equilibrium point for the second 
%                     incurs a non-negative profit value.
% (ix):   PAopt >= 0, indicating that the equilibrium point for the first 
%                     firm incurs a non-negative price value.
% (x):    PBopt >= 0, indicating that the equilibrium point for the second
%                     firm incurs a non-negative price value.
% (xi):   QAopt >= 0, indicating that the equilibrium point for the first 
%                     firm incurs a non-negative quantity value.
% (xii):  QBopt >= 0, indicating that the equlibrium point for the second 
%                     firm incurs a non-negative quantity value.

% FilterFlag is the output variable indicating the exact stage during the 
% filtering process where no valid solutions could be attained.  

% FilterFlag =   0 ==>  Optimization process terminated normally.
% FilterFlag =  -1 ==>  Optimization process did not reach any acceptable 
%                       solutions since the associated ExitFlags were all 
%                       negative: [(i): ExitFlags > 0]
% FilterFlag =  -2 ==>  Optimizaiton process did not reach any acceptable
%                       solutions satisfying the FvalsTolerance threshold:
%                       [(ii): Fvals > FvalsTolerance]
% FilterFlag =  -3 ==>  Optimization process yields solutions that do not
%                       satisfy the First Order Lagrangian Conditions for
%                       the first firm:
%                       [(iii): |Da - DeltaMUa| > DerivativeTolerarnce]
% FilterFlag =  -4 ==>  Optimization process yields solutions that do not 
%                       satisfy the First Order Lagrangian Conditions for  
%                       the second firm:
%                       [(iv): |Db - DeltaMUb| > DerivativeTolerance]
% FilterFlag =  -5 ==>  Optimization process yields solutions that do not
%                       satisfy the Second Order Convexity Conditions for
%                       the first firm: [(v): DDa > 0]
% FilterFlag =  -6 ==>  Optimization process yields solutions that do not
%                       satisfy the Second Order Convexity Conditions for
%                       the second firm: [(vi): DDb > 0]
% FilterFlag =  -7 ==>  Optimization process yields solutions that are
%                       associated with negative profit values for the first
%                       firm: [(vii): FAopt < 0]
% FilterFlag =  -8 ==>  Optimization process yields solutions that are
%                       associated with negative profit values for the
%                       second firm: [(viii): FBopt < 0]
% FilterFlag =  -9 ==>  Optimization process yields solutions that are
%                       associated with negative price values for the first
%                       firm: [(ix): PAopt < 0]
% FilterFlag =  -10 ==> Optimization process yields solutions that are
%                       associated with negative price values for the second
%                       firm: [(x): PBopt < 0]
% FilterFlag = -11 ==>  Optimization process yields solutions that are
%                       associated with negative quantity values for the 
%                       first firm: [(xi): QAopt < 0]
% FilterFlag = -12 ==>  Optimization process yields solutions that are
%                       associated with negative quantity values for the 
%                       second firm: [(xii): QBopt < 0]

% FILTERING PROCESS

% Initialize the output variable.
FilterFlag = 0;
% Initialize a variable storing the number of criterions that need to be
% examined.
CriterionsNumber = 12;
% Get the number of obtained solutions.
SolutionsNumber = size(Solutions,1);
% Initialize a logical matrix storing the evaluation status of each
% solution regarding a particular condition.
I = zeros(SolutionsNumber,CriterionsNumber);
% Initialize a logical matrix storing the cumulative evaluation status for
% each criterion per solution in a repetitive manner.
Status = zeros(SolutionsNumber,CriterionsNumber);

% Isolate the optimal vector of solutions.
Topt = Solutions(:,1:2);
% Isolate the optimal MUa.
% MUa_opt = Solutions(:,3);
% Isolate the optimal MUb.
% MUb_opt = Solutions(:,4);
% Case I: beta > 0
if(beta >= 0)
    % Isolate the optimal MUab.
    MUab_opt = Solutions(:,3);
% Case II: beta < )
else
    % Isoloate the optimal MUab1, MUab2 and MUab3.
    MUab1_opt = Solutions(:,3);
    MUab2_opt = Solutions(:,4);
    MUab3_opt = Solutions(:,5);
end
% Case I: (beta >= 0)
if(beta >= 0)
    % Evaluate the lagrange multipliers differences.
    % DeltaMUa = MUab_opt - MUa_opt;
    % DeltaMUb = MUab_opt - MUb_opt;
    DeltaMUa = MUab_opt;
    DeltaMUb = MUab_opt;
% Case II: (beta < 0)
else
    % Evaluate the lagrange multipliers differences.
    % DeltaMUa = (MUab1_opt - alpha*LB*MUab2_opt - beta*LB*MUab3_opt) - MUa_opt;
    % DeltaMUb = (MUab1_opt - beta*LA*MUab2_opt - alpha*LA*MUab3_opt) - MUb_opt;
    DeltaMUa = (MUab1_opt - alpha*LB*MUab2_opt - beta*LB*MUab3_opt);
    DeltaMUb = (MUab1_opt - beta*LA*MUab2_opt - alpha*LA*MUab3_opt);
end    
% Obtain first and second order derivatives at the retained solution 
% points.
[FD,SD] = OptimalityConditionsTester(Topt,C,G,LA,LB,PA,PB,alpha,beta,gamma);
% Compute fundamental model quantities for each solution point.
[Sopt,Xopt,Popt,Qopt,Fopt] = RetrieveOptimalModelParameters(Solutions(:,1:2),...
                             C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
                         
% Evaluate the status of the obtained solutions regarding condition (i):
I(:,1) = (ExitFlags >= 0);
% Evaluate the status of the obtained solutions regarding condition (ii):
I(:,2) = (Fvals <= FvalsTolerance);
% Evaluate the status of the obtained solutions regarding condition (iii): 
I(:,3) = (abs(FD(:,1) - DeltaMUa) <= DerivativeTolerance);
% Evaluate the status of the obtained solutions regarding condition (iv):
I(:,4) = (abs(FD(:,2) - DeltaMUb) <= DerivativeTolerance);
% Evaluate the status of the obtained solutions regarding condition (v):
I(:,5) = (SD(:,1) < 0);
% Evaluate the status of the obtained solutions regarding condition (vi):
I(:,6) = (SD(:,2) < 0);
% Evaluate the status of the obtained solutions regarding condition (vii):
I(:,7) = (Fopt(:,1) >= 0);
% Evaluate the status of the obtained solutions regarding conditions
% (viii):
I(:,8) = (Fopt(:,2) >= 0);
% Evaluate the status of the obtained solutions regarding condition (ix):
I(:,9) = (Qopt(:,1) >= 0);
% Evaluate the status of the obtained solutions regarding condition (x):
I(:,10) = (Qopt(:,2) >= 0);
% Evaluate the status of the obtained solutions regarding condition (xi):
I(:,11) = (Popt(:,1) >= 0);
% Evaluate the status of the obtained solutions regarding condition (xii):
I(:,12) = (Popt(:,2) >= 0);

% Perform the column-wise conjunction of the boolean values stored in
% matrix I and store the outcome in the boolean vector Io.
Io = ones(SolutionsNumber,1);
for k = 1:CriterionsNumber
    Io = and(Io,I(:,k));
    Status(:,k) = Io;
    if(sum(Io)==0)
        FilterFlag = -k;
        break
    end
end

% Initialize the Filtered output variable.
Filtered = struct;

% Compute optimization criteria on the set of the filtered solutions for
% the case where the resulting filter flag is not equal to 0. Mind that the
% last column index within the Status logical matrix that corresponds to
% valid solutions is -(FilterFlag+1) for FilterFlag <= -2 and (+1) for 
% FilterFlag = -1.
if(FilterFlag~=0)
    % Get the column index pointing to the last no empty set of valid
    % solutions as indicated by the boolean values stored in the Status
    % matrix.
    if(FilterFlag==-1)
        LastIndex = -FilterFlag;
    else
        LastIndex = -(FilterFlag+1);
    end
    % Get the last logical vector of valid indices.
    Ivalid = Status(:,LastIndex);
    % Re-order unfiltered solutions according to the acquired minimum values.
    Filtered.Solutions = Solutions(Ivalid==1,:);
    Filtered.ExitFlags = ExitFlags(Ivalid==1);
    Filtered.Fvals = Fvals(Ivalid==1);
    Filtered.Cineq = Cineq(Ivalid==1,:);
    Filtered.Ceq = Ceq(Ivalid==1,:);
    Topt = Filtered.Solutions(:,1:2);
    % Compute the first and second order derivatives for the acquired solution
    % points.
    [Filtered.FD,Filtered.SD] = OptimalityConditionsTester(Topt,C,G,LA,LB,PA,PB,alpha,beta,gamma);
    [~,Isorted] = sort(Filtered.Fvals,'ascend');
    Filtered.Solutions = Filtered.Solutions(Isorted,:);
    Filtered.ExitFlags = Filtered.ExitFlags(Isorted);
    Filtered.Fvals = Filtered.Fvals(Isorted,:);
    Filtered.Cineq = Filtered.Cineq(Isorted,:);
    Filtered.Ceq = Filtered.Ceq(Isorted,:);
    Filtered.FD = Filtered.FD(Isorted,:);
    Filtered.SD = Filtered.SD(Isorted,:);
    % Retrieve unfiltered optimal model parameters.
    [Filtered.Sopt,Filtered.Xopt,Filtered.Popt,Filtered.Qopt,Filtered.Fopt] = ...
    RetrieveOptimalModelParameters(Filtered.Solutions(:,1:2),C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
    % In case of a negative beta compute the DeltaMUa and DeltaMUb values for
    % each one of the unfiltered solutions.
    if(beta < 0)
        Filtered.MUab1 = Filtered.Solutions(:,3);
        Filtered.MUab2 = Filtered.Solutions(:,4);
        Filtered.MUab3 = Filtered.Solutions(:,5);
        % Filtered.MUa = Filtered.Solutions(:,3);
        % Filtered.MUb = Filtered.Solutions(:,4);
        % Filtered.DeltaMUa = (Filtered.MUab1 - alpha*LB*Filtered.MUab2 - beta*LB*Filtered.MUab3) - Filtered.MUa;
        % Filtered.DeltaMUb = (Filtered.MUab1 - beta*LA*Filtered.MUab2 - alpha*LA*Filtered.MUab3) - Filtered.MUb;
        Filtered.DeltaMUa = (Filtered.MUab1 - alpha*LB*Filtered.MUab2 - beta*LB*Filtered.MUab3);
        Filtered.DeltaMUb = (Filtered.MUab1 - beta*LA*Filtered.MUab2 - alpha*LA*Filtered.MUab3);
        Filtered.DLa = Filtered.FD(:,1) - Filtered.DeltaMUa;
        Filtered.DLb = Filtered.FD(:,2) - Filtered.DeltaMUb;
        Filtered.DL = [Filtered.DLa,Filtered.DLb];
    else
        Filtered.MUab = Filtered.Solutions(:,3);
        Filtered.DeltaMUa = Filtered.MUab;
        Filtered.DeltaMUb = Filtered.MUab;
        Filtered.DLa = Filtered.FD(:,1) - Filtered.DeltaMUa;
        Filtered.DLb = Filtered.FD(:,2) - Filtered.DeltaMUb;
        Filtered.DL = [Filtered.DLa,Filtered.DLb];
    end
end

% Mind that the logical ones of the boolean vector Io identify the
% solutions that satisfy all the required conditions. Therefore, the next
% step is to acquire all the output variables according to the valid
% solution indices stored in the logical vector Io.
Solutions = Solutions(Io,:);
Fvals = Fvals(Io,:);
FD = FD(Io,:);
SD = SD(Io,:);
Fopt = Fopt(Io,:);
Sopt = Sopt(Io,:);
Xopt = Xopt(Io,:);
Popt = Popt(Io,:);
Qopt = Qopt(Io,:);
ExitFlags = ExitFlags(Io);
Ceq = Ceq(Io,:);
Cineq = Cineq(Io,:);