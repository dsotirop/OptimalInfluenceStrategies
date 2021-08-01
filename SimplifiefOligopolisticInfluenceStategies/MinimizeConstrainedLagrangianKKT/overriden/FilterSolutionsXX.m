function [Solutions,Fvals,FilterFlag,FD,SD,Fopt,Sopt,Xopt,Popt,Qopt] = FilterSolutions(Solutions,ExitFlags,Fvals,FvalsTolerance,DerivativeTolerance,C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime)

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

% The obtained solutions will be filtered according to the following
% criteria:
% (i):   ExitFlag >= 0, indicating that the optimization process was
%                       terminated normally.
% (ii):  Fval <= FvalsTolerance, indicating that the optimization process
%                               actually achieved a minimum since the
%                               underlying optimization problem is a
%                               minimization problem.
% (iii):  Da = LambdaA_2 - LambdaA_1, indicating that the equlibrium point
%                                    for the first firm satisfies the  
%                                    First Order Conditions. [Da = FD(1)].
% By letting DLa = LambdaA_2 - LambdaA_1, condition (iii) takes the 
% following form: |Da - DLa| <= DerivativeTolerance.  
% (iv): Db = LambdaB_2 - LambdaB_1, indicating that the equilibrium point
%                                    for the second firm satisfies the
%                                    First Order Conditions. [Da = FD(2)].
% By letting DLb = LambdaB_2 - LambdaB_1, condition (iv) takes the 
% following form: |Db - DLb| <= DerivativeTolerance. 
% (v): IF Da = 0 <=> |Da| <= DerivativeTolerance THEN DDa < 0, indicating
%                                    that the equilibrium point for the
%                                    first firm satisfies the Second Order
%                                    Conditions.
% (vi): IF Db= 0 <=> |Db| <= DerivativeTolerance THEN DDb < 0, indicating 
%                                    that the equilibrium point for the 
%                                    second firm satisfies the Second Order 
%                                    Conditions.
% (vii):   Fopt >= 0,  indicating that the equilibrium point for both firms
%                      incurs a non-negative profit.
% (viii):  Popt >= 0,  indicating that the equilibrium point for both firms 
%                      incurs a non-negative price.
% (ix):    Qopt >= 0,  indicating that the equilibrium point for both firms
%                      incurs a non-negative quantity.

% FilterFlag values indicate the possible reasons for not determining an
% acceptable equilibrium point:
% FilterFlag = 0  ==>  Optimization process terminated with no problems.
% FilterFlag = -1 ==>  Optimization process did not reach any acceptable 
%                      solutions since the associated ExitFlags were all 
%                      negative.
% FilterFlag = -2 ==>  Optimization process yields solutions that do not
%                      satisfy the First Order Lagrangian conditions.
% FilterFlag = -3 ==>  Optimization process yields internal solutions that
%                      do not satisfy the second order conditions.
% FilterFlag = -4 ==>  Optimization process yields solutions that are
%                      associated with negative profit values for at least
%                      one of the firms.
% FilterFlag = -5 ==>  Optimization process yields solutions that are
%                      associated with negative price values for at least 
%                      one of the firms.
% FilterFlag = -6 ==>  Optimization process yields solutions that are
%                      associated with negative quantity values for at
%                      least one of the firms.

% FILTERING PROCESS

% Initialize output variables.
FD = [];
SD = [];
Fopt = [];
Sopt = [];
Xopt = [];
Popt = [];
Qopt = [];
FilterFlag = 0;

% STEP 1: Retain Solutions for which ExitFlag is greater than or equal to 
% zero.
Solutions = Solutions(ExitFlags>=0,:);
Fvals = Fvals(ExitFlags>=0);
% Check whether the Solutions' matrix is not empty.
if(isempty(Solutions))
    FilterFlag = -1;
    % Optimization process failed!
else
    % STEP 2: Retain Solutions for which Fval is less than or equal to 
    % FvalsTolerance.
    Solutions = Solutions(Fvals<=FvalsTolerance,:);
    Fvals = Fvals(Fvals<=FvalsTolerance);
    % Sort Solutions and Fvals in the ascending order of Fvals.
    [Fvals,sorting_order] = sort(Fvals,'ascend');
    Solutions = Solutions(sorting_order,:);
    % Isolate the optimal vector of solutions.
    Topt = Solutions(:,1:2);
    % Isolate the optimal LambdaA.
    LambdaAopt = Solutions(:,3:4);
    % Isolate the optimal LambdaB.
    LambdaBopt = Solutions(:,5:6);
    % Obtain first and second order derivatives at the retained solution 
    % points.
    [FD,SD] = OptimalityConditionsTester(Topt,C,G,LA,LB,PA,PB,alpha,beta,gamma);
    % STEPS 3 and 4: Retain Solutions and auxiliary variables for which the 
    % Lagrangian First Order Conditions are met for the first and second
    % firm.
    DeltaLambdaA = LambdaAopt(:,2) - LambdaAopt(:,1);
    DeltaLambdaB = LambdaBopt(:,2) - LambdaBopt(:,1);
    I1 = (abs(FD(:,1) - DeltaLambdaA) <= DerivativeTolerance);
    I2 = (abs(FD(:,2) - DeltaLambdaB) <= DerivativeTolerance);
    I = and(I1,I2);
    Solutions = Solutions(I,:);
    Fvals = Fvals(I);
    FD = FD(I,:);
    SD = SD(I,:);
    % Check whether First Order Conditions are met for both firms.
    if(isempty(Solutions))
        FilterFlag = -2;
        % First Order Conditions were not met for both firms.
    else
        % STEPS 5 and 6: Retain Solutions and auxiliary variables for which
        % the Second Order Conditions are met.
        % Obtain an index vector to all currently available solutions.
        Io = 1:1:size(Solutions,1);
        % Obtain indices to solution points for which the first order
        % derivatives are zero for the first firm.
        IF1 = (abs(FD(:,1)) <= DerivativeTolerance);
        % Obtain indices to solution points for which the second order
        % derivatives are negative for the first firm.
        IS1 = (SD(:,1) < 0);
        % Obtain indices to solution points that meet the first criterion 
        % and fail to meet to meet the second criterion for the first firm.
        % These solutions are to be excluded in a future step.
        I1 = and(IF1,not(IS1));
        % Obtain indices to solution points for which the first order
        % derivatives are zero for the second firm.
        IF2 = (abs(FD(:,2)) <= DerivativeTolerance);
        % Obtain indices to solution points for which the second order 
        % derivatives are negative.
        IS2 = (SD(:,2) < 0);
        % Obtain indices to solution points that meet the first critering
        % and fail to meet the second criterion for the second firm.
        % These solutions are to be excluded in a future step.
        I2 = and(IF2,not(IS2));
        % Therefore, the remaining solution indices at the completion of
        % elimination STEPS 5 and 6 will be given by removing I1 and I2
        % from Io such that I = Io \ (I1 U I2).
        I = setdiff(Io,union(I1,I2));
        Solutions = Solutions(I,:);
        Fvals = Fvals(I);
        FD = FD(I,:);
        SD = SD(I,:);
        % Check for the existence of remaining solution points after testing
        % the second order solutions for the internal solutions.
        if(isempty(Solutions))
            FilterFlag = -3;
            % Second Order Conditions were not met for both firms at the
            % internal solution points.
        else
            % STEP 7: Retain Solutions and auxiliary variables for which 
            % the profit values for both firms are non-negative.
            % Compute all the remaining model parameters for each solution
            % point.
            [Sopt,Xopt,Popt,Qopt,Fopt] = ...
            RetrieveOptimalModelParameters(Solutions(:,1:2),C,G,LA,LB,PA,PB,alpha,beta,gamma,gamma_prime);
            % Obtain indices to solutions points for which the acquired 
            % profit for the first firm is non-negative.
            I1 = (Fopt(:,1) >= 0);
            % Obtain indices to solution points for which the acquired
            % profit for the second firm is non-negative.
            I2 = (Fopt(:,2) >= 0);
            % Obtain indices to solution points that meet both criteria.
            I = and(I1,I2);
            % Retain corresponding model parameters.
            Solutions = Solutions(I,:);
            Fvals = Fvals(I);
            FD = FD(I,:);
            SD = SD(I,:);
            Sopt = Sopt(I,:);
            Xopt = Xopt(I,:);
            Popt = Popt(I,:);
            Qopt = Qopt(I,:);
            Fopt = Fopt(I,:);
            % Check for the existence of remaining solutions after ensuring
            % the non-negativeness of the overall profit for the two firms.
            if(isempty(Solutions))
                FilterFlag = -4;
                % Retained Solutions involve negative profit values for at
                % least one of the firms.
            else
                % STEP 8: Retain Solutions and auxiliary variables for
                % which the price values for both firms are non-negative.
                % Obtain indices to solution points for which the price
                % value for the first firm is non-negative.
                I1 = (Popt(:,1) >= 0);
                % Obtain indices to solution points for which the price
                % value for the second firm is non-negative.
                I2 = (Popt(:,2) >= 0);
                % Obtain indices to solution points that meet both
                % criteria.
                I = and(I1,I2);
                % Retain corresponding model parameters.
                Solutions = Solutions(I,:);
                Fvals = Fvals(I);
                FD = FD(I,:);
                SD = SD(I,:);
                Sopt = Sopt(I,:);
                Xopt = Xopt(I,:);
                Popt = Popt(I,:);
                Qopt = Qopt(I,:);
                Fopt = Fopt(I,:);
                % Check for the existence of remaining solutions after
                % eliminating solution points associated with negative
                % price values.
                if(isempty(Solutions))
                    FilterFlag = -5;
                    % Retained Solutions involve negative price values for at
                    % least one of the firms.
                else
                    % STEP 9: Retain Solutions and auxiliary variables for 
                    % which the quantity values for both firms are
                    % non-negative.
                    % Obtain indices to solution points for which the
                    % acquired quantity value for the first firm is non-negative.
                    I1 = (Qopt(:,1) >= 0);
                    % Obtain indices to solution points for which the
                    % acquired quantity value for the second firm is non-negative.
                    I2 = (Qopt(:,2) >= 0);
                    % Obtain indices to solution points the meet both
                    % criteria.
                    I = and(I1,I2);
                    % Retain corresponding model parameters.
                    Solutions = Solutions(I,:);
                    Fvals = Fvals(I);
                    FD = FD(I,:);
                    SD = SD(I,:);
                    Sopt = Sopt(I,:);
                    Xopt = Xopt(I,:);
                    Popt = Popt(I,:);
                    Qopt = Qopt(I,:);
                    Fopt = Fopt(I,:);
                    % Check for the existence of remaining solutions after
                    % eliminating solution points associated with negative
                    % quantity values.
                    if(isempty(Solutions))
                        FilterFlag = -6;
                        % Retained Solutions involve negative quantity values
                        % for at least one firm.
                    end
                end
            end
        end
    end

end

