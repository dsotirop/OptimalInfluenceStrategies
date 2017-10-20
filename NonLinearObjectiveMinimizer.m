function [Solution,Fval] = NonLinearObjectiveMinimizer(Dimensionality,N,lb,ub,P1,P2,Theta,Delta,Gamma)

% This function performs numerical optimization of the nonlinear
% objective function defined in the function file "ObectiveFunction"
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
Fobj = @(T)ObjectiveFunction(T,P1,P2,Theta,Delta,Gamma);
Aeq = [];
Beq = [];
%Aeq = [1 1];
%Beq = 1;
%Aeq = [];
%Beq = [];
A = [];
B = [];
for k = 1:N
    [soln(k,:),fvals(k,:),exflags(k,:)] = fmincon(Fobj,pts(k,:),A,B,Aeq,Beq,lb,ub,[],opts);
end;
Solutions = soln;
Fvals = fvals;
ExitFlags = exflags;

SolutionsDiffs = diff(Solutions);
ExitFlagsDiffs = diff(ExitFlags);
SolutionsDiffsSum = sum(sum(SolutionsDiffs));
ExitFlagsDiffsSum = sum(sum(ExitFlagsDiffs));
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

