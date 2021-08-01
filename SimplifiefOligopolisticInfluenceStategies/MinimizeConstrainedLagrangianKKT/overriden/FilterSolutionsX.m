function [Solutions,Fvals,Cineq,Ceq,FD,SD] = FilterSolutions(Solutions,ExitFlags,Fvals,FvalsTolerance,Cineq,Ceq,C,G,LA,LB,PA,PB,alpha,beta,gamma)

% This function filters out the true equibrium points for a given instance 
% of the two-player continuous game that underlies the Simplified 
% Olipolistic Influence model between the two firms (FirmA and FirmB) and 
% the single consumer (Consumer C).

% FvalsTolerance is the maximum acceptable value for the optimization
% parameter Fvals. For this particular problem, the minimization of the
% squared expresssions of the first order KKT conditions, Fvals corresponds
% to the acquired minimum value. Therefore, it should be as small as
% possible. [FvalsTolerance <= 1e-15]

% Initially, we need to retain all solutions for which the associated
% ExitFlags are greater than or equal to zero.

Solutions = Solutions(ExitFlags>=0,:);
Fvals = Fvals(ExitFlags>=0);
Cineq = Cineq(ExitFlags>=0,:);
Ceq = Ceq(ExitFlags>=0,:);

% Extract the optimal solutions points Topt with respect to TA and TB.
Topt = Solutions(:,1:2);

% Compute the first and second order derivatives for the acquired solution
% points.
[FD,SD] = OptimalityConditionsTester(Topt,C,G,LA,LB,PA,PB,alpha,beta,gamma);

% Find the sorted indices of Fvals in ascending order.
[~,Isort_ascend] = sort(Fvals,'ascend');

% Sort all solution-related parameters according to the aorder of Fvals 
% parameter.
Fvals = Fvals(Isort_ascend);
Solutions = Solutions(Isort_ascend,:);
Cineq = Cineq(Isort_ascend,:);
Ceq = Ceq(Isort_ascend,:);
FD = FD(Isort_ascend,:);
SD = SD(Isort_ascend,:);

% Secondly, we need to retain all solution points (TA,TB) for which the
% following second order conditions hold:
% (i):  SDa = DDa < 0, since the underlying optimizaiton problem is a 
%                      maximization problem.
% (ii): SDb = DDb < 0, since the underlying optimization problem is a
%                      maximization problem.
Isd_neg = and((SD(:,1)<0),(SD(:,2)<0));
Fvals = Fvals(Isd_neg);
Solutions = Solutions(Isd_neg,:);
Cineq = Cineq(Isd_neg,:);
Ceq = Ceq(Isd_neg,:);
FD = FD(Isd_neg,:);
SD = SD(Isd_neg,:);

% Finally, we need to reserve all solution points (TA,TB) for which the
% assocociated Fvals values are less than or equal to FvalsTolerance.
Ifval_tol = Fvals<=FvalsTolerance;
Fvals = Fvals(Ifval_tol);
Solutions = Solutions(Ifval_tol,:);
Cineq = Cineq(Ifval_tol,:);
Ceq = Ceq(Ifval_tol,:);
FD = FD(Ifval_tol,:);
SD = SD(Ifval_tol,:);

if(isempty(Fvals))
    fprintf('Equilibrium point could not be found!\n');
end


end

