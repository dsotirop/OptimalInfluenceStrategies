function [Solutions,Fvals,ExitFlags] = SolveNonLinearSystemFmincon(EquationsNum,N,lb,ub)

% This function performs numerical optimization on the nonlinear
% optimization problem defined in the function file "NonLinearSystem"
% by utilizing the fmincon routine.

rng default 
pts = 100*rand(N,EquationsNum);
soln = zeros(N,EquationsNum);
fvals = zeros(N,EquationsNum);
exflags = zeros(N,EquationsNum);
%opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
%opts = optimoptions(@fmincon,'Algorithm','active-set','Display','off');
opts = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
for k = 1:N
    [soln(k,:),fvals(k,:),exflags(k,:)] = fmincon(@(x)0,pts(k,:),[],[],[],[],lb,ub,@fminconstr,opts);
end;
Solutions = soln;
Fvals = fvals;
ExitFlags = exflags;
end