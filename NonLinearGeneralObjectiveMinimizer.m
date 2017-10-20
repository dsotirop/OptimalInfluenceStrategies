function [Solution,Fval] = NonLinearGeneralObjectiveMinimizer(Dimensionality,N,lb,ub,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma)

% This function performs numerical optimization of the nonlinear
% general objective function defined in the function file "GeneralObectiveFunction"
% by utilizing the fmincon routine.

rng default 
pts = 100*rand(N,Dimensionality);
soln = zeros(N,Dimensionality);
fvals = zeros(N,Dimensionality);
exflags = zeros(N,Dimensionality);
%opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
%opts = optimoptions(@fmincon,'Algorithm','active-set','Display','off');
opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
% Set the function handle to the function to be minimized.
Fobj = @(T)GeneralObjectiveFunction(T,P0,P1,P2,Lambda1,Lambda2,Theta1,Theta2,Delta,Gamma);
Aeq = [];
Beq = [];
%Aeq = [1 1];
%Beq = 1;
A = [];
B = [];
for k = 1:N
    %[a,b,c] = fmincon(Fobj,pts(k,:),A,B,Aeq,Beq,lb,ub,[],opts);
    %soln(k,:) = a;
    %if(isempty(b))
        %fvals(k,:) = 0;
    %else
        %fvals(k,:) = b;
    %end;
    %exflags(k,:) = c;
    [soln(k,:),fvals(k,:),exflags(k,:)] = fmincon(Fobj,pts(k,:),A,B,Aeq,Beq,lb,ub,[],opts);
end;
Solutions = soln;
Fvals = fvals;
ExitFlags = exflags;

SolutionsDiffs = diff(Solutions);
ExitFlagsDiffs = diff(ExitFlags);
SolutionsDiffsSum = sum(sum(SolutionsDiffs));
ExitFlagsDiffsSum = sum(sum(ExitFlagsDiffs));
% The following code segment (commented at the moment) provides additional
% control over the minimization process in cases where non-feasible
% solutions cannot be obtained.
%if(sum(ExitFlags)==-2*length(ExitFlags))
    %error('No Feasible Solutions');
%end;
if(SolutionsDiffsSum~=0)
    error('Minimization returned different optimal solutions');
end;
if(ExitFlagsDiffsSum~=0)
    error('Minimization did not succeed');
end;
Solution = Solutions(1,:);
Fval = Fvals(1,:);
end