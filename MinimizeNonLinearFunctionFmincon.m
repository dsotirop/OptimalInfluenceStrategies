function [Solutions,Fvals,ExitFlags] = MinimizeNonLinearFunctionFmincon(Dimensionality,N,lb,ub,P1_const,P2_const,theta_const,A_const,C_const,G_const)

% This function performs numerical optimization of the nonlinear
% objective function defined in the function file "NonLinearFunction"
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
Fobj = @(T)NonLinearFunction(T,P1_const,P2_const,theta_const,A_const,C_const,G_const);
for k = 1:N
    [soln(k,:),fvals(k,:),exflags(k,:)] = fmincon(Fobj,pts(k,:),[],[],[],[],lb,ub,[],opts);
end;
Solutions = soln;
Fvals = fvals;
ExitFlags = exflags;
end

